#!/usr/bin/env python
"""Do the binders self-associate when given a copy of themselves?

Self-association is a liability: a binder that prefers its own surface is less
available for the receptor and often indicates a sticky designed face. The
question is not whether Boltz draws two chains near each other -- with two chains
and nothing else in the box it generally will -- but how confident that interface
is RELATIVE to the same binder's interface with the receptor, which was measured
at ipTM 0.92 median in the p15 complexes.

Reports, per design: the A-B pair ipTM, interface pLDDT, the number of
heavy-atom contacts between the two chains, and how that compares with the same
design's receptor interface.
"""
import csv, glob, json, os, re, statistics as st, sys
from collections import Counter

DIM, MAN, CPLX_A, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
CUT2 = 5.0 ** 2


def heavy(path):
    names, out = [], {}
    for l in open(path):
        if l.startswith("_atom_site."):
            names.append(l.strip().split(".")[1]); continue
        if not (l.startswith("ATOM") or l.startswith("HETATM")): continue
        t = l.split()
        if len(t) != len(names): continue
        d = dict(zip(names, t))
        if (d.get("type_symbol") or "").strip() == "H": continue
        ch = d.get("auth_asym_id") or d.get("label_asym_id")
        try:
            rn = int(d.get("auth_seq_id") or d.get("label_seq_id"))
            xyz = (float(d["Cartn_x"]), float(d["Cartn_y"]), float(d["Cartn_z"]))
        except (TypeError, ValueError, KeyError):
            continue
        out.setdefault(ch, []).append((rn, xyz))
    return out


# receptor-interface reference from the p15 complexes
ref = {}
for f in sorted(glob.glob(os.path.join(CPLX_A, "p15_1099477_sh*_deviation.csv"))):
    pass
cplx_iptm = {}
for d in sorted(glob.glob(os.path.join(CPLX_A, "..", "..", "x"))):
    pass

man = json.load(open(MAN))["monomers"]
rows = []
for m in man:
    dd = os.path.join(DIM, m["mono"].replace("mono_", "dimer_"))
    if not os.path.isdir(dd):
        dd = os.path.join(DIM, m["mono"])
    cif = glob.glob(os.path.join(dd, "*_model_0.cif"))
    conf = glob.glob(os.path.join(dd, "confidence_*.json"))
    if not cif or not conf: continue
    c = json.load(open(conf[0]))
    at = heavy(cif[0])
    chs = sorted(at, key=lambda k: -len(at[k]))
    ncont = nres = 0
    if len(chs) >= 2:
        A, B = at[chs[0]], at[chs[1]]
        bmin = [min(p[1][i] for p in B) - 6 for i in range(3)]
        bmax = [max(p[1][i] for p in B) + 6 for i in range(3)]
        Af = [p for p in A if all(bmin[i] <= p[1][i] <= bmax[i] for i in range(3))]
        rset = set()
        for rn, p in Af:
            for rn2, q in B:
                if (p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2 <= CUT2:
                    ncont += 1; rset.add(rn)
        nres = len(rset)
    pci = c.get("pair_chains_iptm") or {}
    ab = None
    try:
        ks = sorted(pci.keys())
        ab = float(pci[ks[0]][ks[1]])
    except Exception:
        pass
    rows.append({"mono": m["mono"], "n_designs_sharing": m["n_complex_designs"],
                 "members": ";".join(x["name"] for x in m["members"][:3]),
                 "source": m["members"][0]["name"].split("_")[0],
                 "iptm_self": round(float(c.get("iptm", -1)), 4),
                 "pair_iptm_AB": round(ab, 4) if ab is not None else "",
                 "plddt": round(float(c.get("complex_plddt", -1)), 4),
                 "iplddt": round(float(c.get("complex_iplddt", -1)), 4),
                 "ptm": round(float(c.get("ptm", -1)), 4),
                 "contacts_AB": ncont, "iface_residues_A": nres})

if not rows:
    print("nothing analysed"); sys.exit(1)
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

num = lambda k: sorted(r[k] for r in rows if isinstance(r[k], float) and r[k] >= 0)
print(f"analysed {len(rows)} two-chain predictions\n")
print(f"{'metric':22s} {'p10':>8s} {'median':>8s} {'p90':>8s}")
for k in ("iptm_self", "pair_iptm_AB", "iplddt", "plddt", "ptm"):
    v = num(k)
    if v: print(f"{k:22s} {v[len(v)//10]:8.3f} {st.median(v):8.3f} {v[9*len(v)//10]:8.3f}")
cv = sorted(r["contacts_AB"] for r in rows)
rv = sorted(r["iface_residues_A"] for r in rows)
print(f"{'contacts A-B':22s} {cv[len(cv)//10]:8d} {int(st.median(cv)):8d} {cv[9*len(cv)//10]:8d}")
print(f"{'interface residues':22s} {rv[len(rv)//10]:8d} {int(st.median(rv)):8d} {rv[9*len(rv)//10]:8d}")

print("\nCOMPARISON with the receptor interface of the same designs")
print("  p15 binder-receptor ipTM : median 0.923  (from the complex predictions)")
v = num("iptm_self")
print(f"  binder-binder ipTM       : median {st.median(v):.3f}")
print(f"  difference               : {0.923 - st.median(v):+.3f} in favour of the receptor"
      if st.median(v) < 0.923 else
      f"  difference               : {st.median(v) - 0.923:+.3f} in favour of SELF - a liability")

print("\nhow many look like confident self-dimers?")
for thr in (0.80, 0.85, 0.90):
    n = sum(1 for x in v if x >= thr)
    print(f"  ipTM(self) >= {thr:.2f} : {n:4d}/{len(v)} ({100*n/len(v):3.0f}%)")
low = [r for r in rows if isinstance(r["iptm_self"], float) and r["iptm_self"] < 0.60]
print(f"  ipTM(self) <  0.60 : {len(low):4d}/{len(rows)} "
      f"({100*len(low)/max(len(rows),1):3.0f}%)  <- least self-associating")
print("\nby source (median self ipTM, median A-B contacts):")
for s in ("p3", "p10", "p12", "p13"):
    sub = [r for r in rows if r["source"] == s]
    if not sub: continue
    sv = [r["iptm_self"] for r in sub if isinstance(r["iptm_self"], float)]
    print(f"  {s:4s} n={len(sub):3d}  ipTM {st.median(sv):.3f}  "
          f"contacts {int(st.median([r['contacts_AB'] for r in sub]))}")
print(f"\nwrote {OUT}")
