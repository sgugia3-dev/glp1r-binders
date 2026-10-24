#!/usr/bin/env python
"""Which receptor regions does each binder segment contact, and what is their charge?

The instruction is to make residues above 10 complementary to the ECD. That is
worth verifying rather than assuming: if binder 11-140 really contacts the ECD
(acidic, net -5.6, pI 4.88) while the 1-10 anchor contacts the TMD (basic, net
+6.6, pI 8.96), then a positive designed region and an acidic anchor are both
correct, and the charge targets follow from the measurement.

Reports, per packaged set, the receptor residues contacted by binder 1-10 and by
binder 11-140 separately, split into ECD (1-117) and TMD (118-409).
"""
import glob, os, statistics as st, sys, warnings
from collections import Counter
warnings.filterwarnings("ignore")
from Bio.PDB import MMCIFParser, NeighborSearch

G = sys.argv[1]
CUT = 5.0
ECD_MAX = 117
AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
SETS = [("p3_best", "P3"), ("p10_best", "P10"), ("p12_best", "P12"), ("p13_best", "P13")]
parser = MMCIFParser(QUIET=True)


def q(seq):
    c = Counter(seq)
    return c["K"] + c["R"] - c["D"] - c["E"]


print(f"{'set':5s} {'binder seg':12s} {'region':10s} {'n_res':>6s} {'charge':>7s} "
      f"{'f_neg':>6s} {'f_pos':>6s} {'f_hydr':>7s}")
agg = {}
for d, lbl in SETS:
    cifs = sorted(glob.glob(os.path.join(G, d, "structures", "*.cif")))
    if not cifs: continue
    acc = {(seg, reg): [] for seg in ("1-10", "11-140") for reg in ("ECD", "TMD")}
    for f in cifs:
        m = parser.get_structure("x", f)[0]
        ch = sorted(m, key=lambda c: -len([r for r in c if r.id[0] == " "]))
        tgt, bnd = ch[0], ch[1]
        tat = [a for a in tgt.get_atoms() if a.element != "H"]
        ns_t = NeighborSearch(tat)
        tres = {r.id[1]: AA3.get(r.get_resname().upper(), "X") for r in tgt if r.id[0] == " "}
        off = min(tres) - 1                       # target renumbered after the binder
        for seg, lo, hi in (("1-10", 1, 10), ("11-140", 11, 140)):
            hitset = set()
            for a in bnd.get_atoms():
                if a.element == "H": continue
                rn = a.get_parent().id[1]
                if not (lo <= rn <= hi): continue
                for b in ns_t.search(a.coord, CUT):
                    hitset.add(b.get_parent().id[1])
            for reg, rlo, rhi in (("ECD", 1, ECD_MAX), ("TMD", ECD_MAX + 1, 409)):
                s = "".join(tres[i] for i in sorted(hitset)
                            if i in tres and rlo <= (i - off) <= rhi)
                acc[(seg, reg)].append(s)
    for seg in ("1-10", "11-140"):
        for reg in ("ECD", "TMD"):
            ss = acc[(seg, reg)]
            ns = [len(x) for x in ss]
            if not any(ns):
                print(f"{lbl:5s} {seg:12s} {reg:10s} {0:6.0f} {'-':>7s} {'-':>6s} {'-':>6s} {'-':>7s}")
                continue
            qs = [q(x) for x in ss]
            pooled = "".join(ss)
            c = Counter(pooled); n = len(pooled) or 1
            print(f"{lbl:5s} {seg:12s} {reg:10s} {st.median(ns):6.0f} {st.median(qs):7.1f} "
                  f"{(c['D']+c['E'])/n:6.3f} {(c['K']+c['R'])/n:6.3f} "
                  f"{sum(c[a] for a in 'AVLIMFWC')/n:7.3f}")
            agg[(lbl, seg, reg)] = (st.median(ns), st.median(qs))
    print()

print("=" * 78)
print("CONCLUSION for the charge target")
print("=" * 78)
for seg in ("11-140", "1-10"):
    e = [v[1] for k, v in agg.items() if k[1] == seg and k[2] == "ECD"]
    t = [v[1] for k, v in agg.items() if k[1] == seg and k[2] == "TMD"]
    ne = [v[0] for k, v in agg.items() if k[1] == seg and k[2] == "ECD"]
    nt = [v[0] for k, v in agg.items() if k[1] == seg and k[2] == "TMD"]
    print(f"  binder {seg:7s} contacts  ECD: {st.median(ne) if ne else 0:5.0f} res, "
          f"charge {st.median(e) if e else 0:+5.1f}   |   TMD: {st.median(nt) if nt else 0:5.0f} res, "
          f"charge {st.median(t) if t else 0:+5.1f}")
print("\n  complementarity requires the binder segment to carry the OPPOSITE sign")
print("  to the surface it touches.")
