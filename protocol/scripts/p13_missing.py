#!/usr/bin/env python
"""Which hotspot do the 3-of-4 designs drop?

P13 added A95 and A27 fell from 80% to 43%, while A42/A86/A95 all sit above 90%.
If the 619 designs at exactly 3/4 are almost all missing A27, then the set is
really a strong 3-point epitope plus one unreachable point, and the fix is to
re-specify or drop A27 rather than to steer harder.
"""
import csv, glob, os, re, statistics as st, sys
from collections import Counter

A, RUN = sys.argv[1], sys.argv[2]
HOTS = ["A27", "A42", "A86", "A95"]


def load(pat):
    out = {}
    for f in sorted(glob.glob(os.path.join(A, f"{RUN}_sh*_{pat}.csv"))):
        sh = re.search(r"_sh(\d+)_", f).group(1)
        for r in csv.DictReader(open(f)): out[(sh, r["name"])] = r
    return out


dev, epi = load("deviation"), load("epitope")
keys = sorted(set(dev) & set(epi))
rows = []
for k in keys:
    e, d = epi[k], dev[k]
    hit = {h: (float(e.get(f"d_des_{h}", 99) or 99) < 5.0) for h in HOTS}
    rows.append({"hits": sum(hit.values()), "hit": hit,
                 "motif": float(d["motif_rmsd_vs_af3"])})
n = len(rows)
print(f"{n} refolds\n")

for want in (3, 4, 2):
    sub = [r for r in rows if r["hits"] == want]
    if not sub: continue
    print(f"designs with exactly {want}/4 hotspots: {len(sub)} ({100*len(sub)/n:.0f}%)")
    miss = Counter()
    for r in sub:
        for h in HOTS:
            if not r["hit"][h]: miss[h] += 1
    for h, c in miss.most_common():
        print(f"    missing {h}: {c:4d} ({100*c/len(sub):3.0f}% of this group)")
    combo = Counter("+".join(h for h in HOTS if r["hit"][h]) for r in sub)
    for c, k2 in combo.most_common(3):
        print(f"    set {c:28s} {k2:4d}")
    print()

print("if A27 were dropped from the set, how many designs hold A42+A86+A95?")
triple = [r for r in rows if r["hit"]["A42"] and r["hit"]["A86"] and r["hit"]["A95"]]
print(f"  {len(triple)}/{n} ({100*len(triple)/n:.0f}%)")
print(f"  compare P12's 3-point A27/42/86 full engagement: 75%")
print(f"  and the current 4-point full engagement:         {100*sum(1 for r in rows if r['hits']==4)/n:.0f}%")

print("\nanchor placement, by whether A27 was reached")
for lbl, sel in (("A27 reached", [r for r in rows if r["hit"]["A27"]]),
                 ("A27 missed", [r for r in rows if not r["hit"]["A27"]])):
    m = [r["motif"] for r in sel]
    print(f"  {lbl:14s} n={len(sel):4d}  median motif RMSD {st.median(m):.2f} A  "
          f"under 2A: {100*sum(1 for x in m if x<2)/len(m):3.0f}%")
