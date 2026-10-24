#!/usr/bin/env python
"""How extended are the P1 designs? rfd3 already records radius_of_gyration per design.

Reference: a compact globular protein has Rg ~ 2.2 * N^0.38 A. The ratio
observed/expected is the compactness score; ~1.0 is globular, >1.3 is extended.
"""
import glob, json, os, statistics as st, sys

rows = []
for j in glob.glob(os.path.join(sys.argv[1], "sh*", "diff", "*.json")):
    try: d = json.load(open(j))
    except Exception: continue
    m = d.get("metrics", {}); sp = d.get("specification", {})
    sc = sp.get("extra", {}).get("sampled_contig", "")
    im = d.get("diffused_index_map", {})
    nb = len([v for k, v in im.items() if k.startswith("B") and v.startswith("A")])
    flanks = sum(int(x[:-1]) for x in sc.split(",/0,")[0].split(",") if x.endswith("P"))
    L = flanks + nb
    rg = m.get("radius_of_gyration")
    if not rg or not L: continue
    exp = 2.2 * (L ** 0.38)
    rows.append({"L": L, "rg": rg, "exp": exp, "ratio": rg / exp,
                 "helix": m.get("helix_fraction", 0), "ss": m.get("num_ss_elements", 0),
                 "dev": m.get("max_ca_deviation", 99), "brk": m.get("n_chainbreaks", 9)})

print(f"designs with metrics: {len(rows)}")
if not rows: sys.exit(1)
r = [x["ratio"] for x in rows]
print(f"\nRg / expected-globular:  min {min(r):.2f}  median {st.median(r):.2f}  max {max(r):.2f}")
print(f"  compact   (<=1.15): {sum(1 for x in r if x<=1.15):4d} ({100*sum(1 for x in r if x<=1.15)/len(r):.0f}%)")
print(f"  moderate  (<=1.30): {sum(1 for x in r if x<=1.30):4d} ({100*sum(1 for x in r if x<=1.30)/len(r):.0f}%)")
print(f"  EXTENDED   (>1.30): {sum(1 for x in r if x>1.30):4d} ({100*sum(1 for x in r if x>1.30)/len(r):.0f}%)")
print(f"  very ext   (>1.60): {sum(1 for x in r if x>1.60):4d} ({100*sum(1 for x in r if x>1.60)/len(r):.0f}%)")
print(f"\nabsolute Rg: median {st.median([x['rg'] for x in rows]):.1f} A, "
      f"length median {st.median([x['L'] for x in rows]):.0f}")
print("\nby binder length:")
print("%8s %5s %8s %8s %7s" % ("len band", "n", "Rg med", "exp", "ratio"))
for lo, hi in ((50,65),(66,75),(76,85),(86,95),(96,110)):
    s = [x for x in rows if lo <= x["L"] <= hi]
    if len(s) < 5: continue
    print("%8s %5d %8.1f %8.1f %7.2f" % (f"{lo}-{hi}", len(s),
          st.median([x["rg"] for x in s]), st.median([x["exp"] for x in s]),
          st.median([x["ratio"] for x in s])))
# does compactness co-vary with anything we already screen on?
comp = [x for x in rows if x["ratio"] <= 1.15]; ext = [x for x in rows if x["ratio"] > 1.30]
if comp and ext:
    print(f"\ncompact  (n={len(comp)}): helix {st.median([x['helix'] for x in comp]):.2f}  "
          f"SS elems {st.median([x['ss'] for x in comp]):.0f}  "
          f"motif-held {sum(1 for x in comp if x['dev']<=2 and x['brk']<=1)}/{len(comp)}")
    print(f"extended (n={len(ext)}): helix {st.median([x['helix'] for x in ext]):.2f}  "
          f"SS elems {st.median([x['ss'] for x in ext]):.0f}  "
          f"motif-held {sum(1 for x in ext if x['dev']<=2 and x['brk']<=1)}/{len(ext)}")
