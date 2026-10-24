#!/usr/bin/env python
"""Count salt bridges, hydrogen bonds and hydrophobic contacts across the interface.

Definitions, stated because each is a distance-only proxy:
  salt bridge   charged N of K/R/H within 4.0 A of carboxylate O of D/E
  hydrogen bond N/O donor-acceptor pair within 3.5 A. NO ANGLE TERM is applied,
                so this over-counts relative to a real H-bond detector; use it to
                compare designs with each other, not as an absolute count.
  hydrophobic   side-chain carbon of A/V/L/I/M/F/W/P within 4.5 A of the same on
                the other side
Only binder-to-receptor pairs are counted; intramolecular contacts are ignored.
Reported separately for the whole binder, the designed region (11-140) and the
redesigned anchor (1-10), since the anchor contacts the TMD and the designed
region the ECD.
"""
import argparse, csv, glob, json, os, statistics as st, sys
from collections import Counter

ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
CATION = {"LYS": {"NZ"}, "ARG": {"NE", "NH1", "NH2"}, "HIS": {"ND1", "NE2"}}
ANION = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}
APOLAR_RES = set("AVLIMFWP")
BACKBONE = {"N", "CA", "C", "O"}
ANCHOR_HI = 10


def read_cif(path):
    names, rows = [], []
    for l in open(path):
        if l.startswith("_atom_site."):
            names.append(l.strip().split(".")[1]); continue
        if not (l.startswith("ATOM") or l.startswith("HETATM")): continue
        t = l.split()
        if len(t) != len(names): continue
        d = dict(zip(names, t))
        el = (d.get("type_symbol") or "").strip()
        if el == "H": continue
        try:
            rn = int(d.get("auth_seq_id") or d.get("label_seq_id"))
            xyz = (float(d["Cartn_x"]), float(d["Cartn_y"]), float(d["Cartn_z"]))
        except (TypeError, ValueError, KeyError):
            continue
        rows.append({"ch": d.get("auth_asym_id") or d.get("label_asym_id"), "rn": rn,
                     "res": (d.get("label_comp_id") or "").upper(),
                     "at": (d.get("label_atom_id") or "").strip().strip('"'),
                     "el": el, "xyz": xyz})
    return rows


def d2(p, q):
    return (p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2


def analyse(path):
    at = read_cif(path)
    by = {}
    for r in at: by.setdefault(r["ch"], []).append(r)
    if len(by) < 2: return None
    chs = sorted(by, key=lambda c: -len({r["rn"] for r in by[c]}))
    tgt, bnd = chs[0], chs[1]
    B, T = by[bnd], by[tgt]
    # coarse prefilter: only receptor atoms near the binder matter
    bmin = [min(r["xyz"][i] for r in B) - 6 for i in range(3)]
    bmax = [max(r["xyz"][i] for r in B) + 6 for i in range(3)]
    T = [r for r in T if all(bmin[i] <= r["xyz"][i] <= bmax[i] for i in range(3))]

    out = Counter()
    sb_pairs, hb_pairs = [], []
    for b in B:
        seg = "anchor" if b["rn"] <= ANCHOR_HI else "designed"
        b_cat = b["at"] in CATION.get(b["res"], ())
        b_ani = b["at"] in ANION.get(b["res"], ())
        b_apo = (ONE.get(b["res"]) in APOLAR_RES and b["el"] == "C"
                 and b["at"] not in BACKBONE)
        b_pol = b["el"] in ("N", "O")
        for t in T:
            dd = d2(b["xyz"], t["xyz"])
            if dd > 20.25: continue                      # 4.5 A, the widest cutoff
            if dd <= 16.0:                               # 4.0 A salt bridge
                if (b_cat and t["at"] in ANION.get(t["res"], ())) or \
                   (b_ani and t["at"] in CATION.get(t["res"], ())):
                    out["salt_bridges"] += 1; out[f"sb_{seg}"] += 1
                    sb_pairs.append(f"{ONE.get(b['res'],'X')}{b['rn']}-"
                                    f"{ONE.get(t['res'],'X')}{t['rn']}")
            if dd <= 12.25 and b_pol and t["el"] in ("N", "O"):   # 3.5 A H-bond
                out["hbonds"] += 1; out[f"hb_{seg}"] += 1
            if b_apo and ONE.get(t["res"]) in APOLAR_RES and t["el"] == "C" \
                    and t["at"] not in BACKBONE:
                out["hydrophobic_contacts"] += 1; out[f"hc_{seg}"] += 1
    bseq = {r["rn"]: ONE.get(r["res"], "X") for r in B}
    des = "".join(bseq[k] for k in sorted(bseq) if k > ANCHOR_HI)
    anc = "".join(bseq[k] for k in sorted(bseq) if k <= ANCHOR_HI)
    c = Counter(des)
    out["designed_charge"] = c["K"] + c["R"] - c["D"] - c["E"]
    out["binder_len"] = len(bseq)
    out["anchor_seq"] = anc
    out["unique_salt_bridges"] = len(set(sb_pairs))
    out["salt_bridge_list"] = ";".join(sorted(set(sb_pairs))[:14])
    return out


ap = argparse.ArgumentParser()
ap.add_argument("--pred-dir", required=True)
ap.add_argument("--out-csv", required=True)
ap.add_argument("--flat", action="store_true",
                help="pred-dir holds *.cif directly instead of one "
                     "subdirectory per prediction")
a = ap.parse_args()

if a.flat:
    targets = [(os.path.basename(c).replace("_complex.cif", "").replace(".cif", ""), c)
               for c in sorted(glob.glob(os.path.join(a.pred_dir, "*.cif")))]
else:
    targets = []
    for d in sorted(glob.glob(os.path.join(a.pred_dir, "*"))):
        if not os.path.isdir(d): continue
        cif = glob.glob(os.path.join(d, "*_model_0.cif"))
        if cif: targets.append((os.path.basename(d), cif[0]))
rows = []
for nm, cifpath in targets:
    r = analyse(cifpath)
    if r is None: continue
    rec = {"name": nm}
    for k in ("salt_bridges", "unique_salt_bridges", "sb_designed", "sb_anchor",
              "hbonds", "hb_designed", "hb_anchor", "hydrophobic_contacts",
              "hc_designed", "hc_anchor", "designed_charge", "binder_len"):
        rec[k] = r.get(k, 0)
    rec["anchor_seq"] = r.get("anchor_seq", "")
    rec["salt_bridge_list"] = r.get("salt_bridge_list", "")
    rows.append(rec)

if not rows:
    print("no structures analysed"); sys.exit(0)
cols = list(rows[0])
os.makedirs(os.path.dirname(a.out_csv) or ".", exist_ok=True)
with open(a.out_csv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)

def med(k): return st.median([r[k] for r in rows])
print(f"analysed {len(rows)} structures")
print(f"  salt bridges (unique pairs) : median {med('unique_salt_bridges'):.0f}  "
      f"max {max(r['unique_salt_bridges'] for r in rows)}")
print(f"    of those in designed region: median {med('sb_designed'):.0f}")
print(f"  H-bond atom pairs <3.5A     : median {med('hbonds'):.0f}  (no angle term)")
print(f"  hydrophobic C-C <4.5A       : median {med('hydrophobic_contacts'):.0f}")
print(f"  q(11-140)                   : median {med('designed_charge'):+.0f}  "
      f"min {min(r['designed_charge'] for r in rows):+.0f}  "
      f"max {max(r['designed_charge'] for r in rows):+.0f}")
anc = Counter(r["anchor_seq"] for r in rows)
print(f"  distinct redesigned anchors : {len(anc)}")
for s, n in anc.most_common(5):
    print(f"      {s}  x{n}")
print(f"  pos1==H / pos3 in D,E       : "
      f"{sum(1 for r in rows if r['anchor_seq'][:1]=='H')}/{len(rows)} , "
      f"{sum(1 for r in rows if r['anchor_seq'][2:3] in ('D','E'))}/{len(rows)}")
print(f"wrote {a.out_csv}")
