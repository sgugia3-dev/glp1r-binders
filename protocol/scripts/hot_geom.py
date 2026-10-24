#!/usr/bin/env python
"""Can one binder reach hotspots A27, A42 and A86 at once?

Compares the span of the requested combined set against the two patches that
already worked (P2 24/25/27/42 and P3 86-89), and reports each hotspot's
distance from the anchor, since the anchor position is fixed.
"""
import itertools, math, sys

pdb = sys.argv[1]
WANT = {27: "CD1", 42: "OG1", 86: "CD2"}   # the atoms P2 and P3 actually steered on

atoms = {}
pep = []
for l in open(pdb):
    if not l.startswith("ATOM"): continue
    rn, an, ch = int(l[22:26]), l[12:16].strip(), l[21]
    xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    if ch == "A": atoms.setdefault(rn, {})[an] = xyz
    elif ch == "B" and an == "CA": pep.append((rn, xyz))

d = lambda p, q: math.dist(p, q)
print("requested hotspot atoms")
sel = {}
for r, a in sorted(WANT.items()):
    res = atoms.get(r, {})
    status = "found" if a in res else "MISSING -> falling back to CA"
    sel[r] = res.get(a) or res.get("CA")
    print(f"  A{r} {res.get('_rn','')} want {a:4s} -> {status}")

print("\npairwise distances, selected atoms")
for x, y in itertools.combinations(sorted(sel), 2):
    print(f"  A{x} - A{y}: {d(sel[x], sel[y]):6.1f} A")

span = max(d(sel[x], sel[y]) for x, y in itertools.combinations(sorted(sel), 2))
print(f"\nmax span across the requested set: {span:.1f} A")
print("reference spans for the patches that worked (CA):")
for nm, grp in (("P2 24/25/27/42", [24, 25, 27, 42]), ("P3 86-89", [86, 87, 88, 89])):
    pts = [atoms[r]["CA"] for r in grp if r in atoms]
    s = max(d(p, q) for p, q in itertools.combinations(pts, 2))
    print(f"  {nm:16s} {s:5.1f} A")

print("\ndistance from the anchor (chain B 1-10) to each hotspot, nearest anchor CA:")
for r in sorted(sel):
    print(f"  A{r}: {min(d(sel[r], p) for _, p in pep):6.1f} A")

print("\nresidue identities")
AA = {}
for l in open(pdb):
    if l.startswith("ATOM") and l[21] == "A" and l[12:16].strip() == "CA":
        AA[int(l[22:26])] = l[17:20].strip()
for r in sorted(set(list(WANT) + [24, 25, 87, 88, 89])):
    print(f"  A{r}: {AA.get(r,'?')}" + ("   <- requested" if r in WANT else ""))
