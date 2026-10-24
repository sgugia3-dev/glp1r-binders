#!/usr/bin/env python
"""Interface composition and charge complementarity on the packaged best sets.

Sequence-level charge describes the whole binder; what governs binding is the
patch that actually touches the receptor. For each of the 250 packaged complexes
this takes the binder residues with any heavy atom within 5 A of any receptor
heavy atom (and vice versa) and reports the composition and net charge of both
sides, so binder-vs-epitope complementarity can be read directly.

The anchor motif (binder residues 1-10) is reported separately from the designed
interface, since it is identical in every design and always buried in the TMD.
"""
import csv, glob, os, statistics as st, sys, warnings
from collections import Counter
warnings.filterwarnings("ignore")
from Bio.PDB import MMCIFParser, NeighborSearch

G, OUT = sys.argv[1], sys.argv[2]
CUT = 5.0
AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
CLASSES = {"hydrophobic": set("AVLIMFWC"), "polar": set("STNQYH"),
           "positive": set("KR"), "negative": set("DE"), "special": set("GP")}
SETS = [("p2_best", "P2"), ("p3_best", "P3"), ("p10_best", "P10"),
        ("p12_best", "P12"), ("p13_best", "P13")]
parser = MMCIFParser(QUIET=True)


def compo(seq):
    n = len(seq) or 1
    c = Counter(seq)
    d = {f"f_{k}": sum(c[a] for a in s) / n for k, s in CLASSES.items()}
    d["net_charge"] = sum(c[a] for a in "KR") - sum(c[a] for a in "DE")
    d["n_res"] = len(seq)
    return d


rows = []
for d, lbl in SETS:
    cifs = sorted(glob.glob(os.path.join(G, d, "structures", "*.cif")))
    if not cifs:
        print(f"{lbl}: no structures"); continue
    for f in cifs:
        s = parser.get_structure("x", f)[0]
        chains = sorted(s, key=lambda c: -len([r for r in c if r.id[0] == " "]))
        tgt, bnd = chains[0], chains[1]        # 409 then 140
        tat = [a for a in tgt.get_atoms() if a.element != "H"]
        bat = [a for a in bnd.get_atoms() if a.element != "H"]
        ns_t, ns_b = NeighborSearch(tat), NeighborSearch(bat)
        bi, ti = set(), set()
        for a in bat:
            if ns_t.search(a.coord, CUT):
                bi.add(a.get_parent().id[1])
        for a in tat:
            if ns_b.search(a.coord, CUT):
                ti.add(a.get_parent().id[1])
        bres = {r.id[1]: AA3.get(r.get_resname().upper(), "X") for r in bnd if r.id[0] == " "}
        tres = {r.id[1]: AA3.get(r.get_resname().upper(), "X") for r in tgt if r.id[0] == " "}
        off = min(tres)                        # target renumbered after the binder
        b_all = "".join(bres[i] for i in sorted(bi) if i in bres)
        b_des = "".join(bres[i] for i in sorted(bi) if i in bres and i > 10)
        b_anc = "".join(bres[i] for i in sorted(bi) if i in bres and i <= 10)
        t_all = "".join(tres[i] for i in sorted(ti) if i in tres)
        r = {"set": lbl, "design": os.path.basename(f).replace("_complex.cif", "")}
        for tag, seq in (("bind_iface", b_all), ("bind_designed_iface", b_des),
                         ("bind_anchor_iface", b_anc), ("recep_iface", t_all)):
            for k, v in compo(seq).items():
                r[f"{tag}_{k}"] = round(v, 4) if isinstance(v, float) else v
        r["recep_iface_resnums"] = ";".join(str(i - off + 1) for i in sorted(ti))
        rows.append(r)
    print(f"{lbl}: {len(cifs)} complexes analysed")

cols = ["set", "design"] + [c for c in rows[0] if c not in ("set", "design")]
with open(os.path.join(OUT, "interface_properties.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
print(f"\nwrote interface_properties.csv ({len(rows)} rows)")

print("\n" + "=" * 100)
print(f"INTERFACE COMPOSITION, median over 50 designs per set (heavy-atom cutoff {CUT} A)")
print("=" * 100)
print(f"{'set':5s} {'side':22s} {'nres':>5s} {'charge':>7s} {'hydrophob':>10s} "
      f"{'polar':>7s} {'pos':>6s} {'neg':>6s}")
for d, lbl in SETS:
    sel = [r for r in rows if r["set"] == lbl]
    if not sel: continue
    for tag, nm in (("bind_designed_iface", "binder designed face"),
                    ("bind_anchor_iface", "binder anchor (fixed)"),
                    ("recep_iface", "receptor epitope")):
        f = lambda k: st.median([r[f"{tag}_{k}"] for r in sel])
        print(f"{lbl:5s} {nm:22s} {f('n_res'):5.0f} {f('net_charge'):7.1f} "
              f"{f('f_hydrophobic'):10.3f} {f('f_polar'):7.3f} "
              f"{f('f_positive'):6.3f} {f('f_negative'):6.3f}")
    print()

print("=" * 100)
print("CHARGE COMPLEMENTARITY: binder designed face vs receptor epitope")
print("=" * 100)
print(f"{'set':5s} {'binder q':>9s} {'receptor q':>11s} {'sum':>7s} "
      f"{'opposite sign?':>15s}  interpretation")
for d, lbl in SETS:
    sel = [r for r in rows if r["set"] == lbl]
    if not sel: continue
    bq = st.median([r["bind_designed_iface_net_charge"] for r in sel])
    tq = st.median([r["recep_iface_net_charge"] for r in sel])
    opp = (bq > 0) != (tq > 0) and bq != 0 and tq != 0
    interp = ("complementary" if opp else
              "like-charged - electrostatically unfavourable" if bq * tq > 0 else "one side neutral")
    print(f"{lbl:5s} {bq:9.1f} {tq:11.1f} {bq+tq:7.1f} {str(opp):>15s}  {interp}")
