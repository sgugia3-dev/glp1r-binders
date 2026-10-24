#!/usr/bin/env python
"""Choose the hotspot atom for Leu95, calibrated against atoms of known performance.

This campaign's measured result is that burial, not distance, limits engagement:
  A42:OG1  96% engaged      A86:CD2  94%      A27:CD1  80% (the weak one)
So the useful question is not "is Leu95 near the others" but "which of its atoms
is as exposed as OG1/CD2 rather than as buried as CD1".

Burial proxy: count of OTHER receptor heavy atoms within 8 A (excluding the
residue's own atoms). Lower = more exposed. Calibrated on the three known atoms.
"""
import itertools, math, sys
from collections import defaultdict

pdb = sys.argv[1]
KNOWN = {(27, "CD1"): "80% (weak)", (42, "OG1"): "96%", (86, "CD2"): "94%"}
TARGETS = [95]

res = defaultdict(dict)
allatoms = []
pep = []
for l in open(pdb):
    if not l.startswith("ATOM"): continue
    ch, rn, an = l[21], int(l[22:26]), l[12:16].strip()
    el = l[76:78].strip()
    if el == "H": continue
    xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    if ch == "A":
        res[rn][an] = xyz
        allatoms.append((rn, an, xyz))
    elif ch == "B" and an == "CA":
        pep.append(xyz)

name = {}
for l in open(pdb):
    if l.startswith("ATOM") and l[21] == "A" and l[12:16].strip() == "CA":
        name[int(l[22:26])] = l[17:20].strip()

d = lambda p, q: math.dist(p, q)


def burial(rn, an):
    p = res[rn][an]
    return sum(1 for r2, a2, x2 in allatoms if r2 != rn and d(p, x2) < 8.0)


print(f"residue 95 identity: {name.get(95)}")
if name.get(95) != "LEU":
    print(f"  WARNING: residue 95 is {name.get(95)}, not LEU - check the numbering")
print()
print("CALIBRATION: burial of the atoms whose engagement rate we measured")
print(f"  {'atom':12s} {'burial':>7s}  measured engagement")
cal = []
for (rn, an), perf in KNOWN.items():
    b = burial(rn, an)
    cal.append((b, perf))
    print(f"  A{rn}:{an:8s} {b:7d}  {perf}")
good = [b for b, p in cal if "weak" not in p]
bad = [b for b, p in cal if "weak" in p]
print(f"  -> exposed/high-performing: {good}, buried/weak: {bad}")

print(f"\nCANDIDATE ATOMS on residue 95 ({name.get(95)}), least buried first")
cands = sorted(res[95], key=lambda a: burial(95, a))
print(f"  {'atom':8s} {'burial':>7s} {'d(A27)':>8s} {'d(A42)':>8s} {'d(A86)':>8s} "
      f"{'d(anchor)':>10s}  verdict")
ref = {27: res[27]["CD1"], 42: res[42]["OG1"], 86: res[86]["CD2"]}
for an in cands:
    if an in ("N", "C", "O"): continue
    b = burial(95, an)
    ds = {k: d(res[95][an], v) for k, v in ref.items()}
    da = min(d(res[95][an], p) for p in pep)
    v = "as exposed as the good ones" if b <= max(good) else (
        "more buried than A27:CD1" if b > max(bad) else "intermediate")
    print(f"  {an:8s} {b:7d} {ds[27]:8.1f} {ds[42]:8.1f} {ds[86]:8.1f} {da:10.1f}  {v}")

print("\nSET GEOMETRY with residue 95 added (using its least-buried side-chain atom)")
sc = [a for a in cands if a not in ("N", "C", "O", "CA")]
best = sc[0] if sc else "CA"
sel = {27: res[27]["CD1"], 42: res[42]["OG1"], 86: res[86]["CD2"], 95: res[95][best]}
print(f"  chosen atom for 95: {best}")
for x, y in itertools.combinations(sorted(sel), 2):
    print(f"  A{x} - A{y}: {d(sel[x], sel[y]):6.1f} A")
span = max(d(sel[x], sel[y]) for x, y in itertools.combinations(sorted(sel), 2))
print(f"  max span: {span:.1f} A")
print(f"  reference spans: P12 A27/42/86 = 16.5 A, P2 A24/25/27/42 = 11.6 A, P3 = 8.5 A")
print(f"  anchor to A95:{best}: {min(d(res[95][best], p) for p in pep):.1f} A")
