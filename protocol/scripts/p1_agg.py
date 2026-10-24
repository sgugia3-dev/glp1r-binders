#!/usr/bin/env python
"""Aggregate the 20 shards: motif deviation + hotspot engagement, and rank candidates."""
import csv, glob, statistics as st, sys, json

# rfd3 names a design {json_basename}_{key}_{batch}_{model}, so with one shared input JSON
# EVERY shard emits the same 90 names and only the seed differs. Keying on the name alone
# collapsed 1800 measurements into 90. Key on (shard, name) instead; the shard comes from
# the filename. Fix for next time: pass a per-shard `global_prefix=` to rfd3.
import re
dev, epi = {}, {}
def shard_of(path):
    m = re.search(r"_sh(\d+)_", path)
    return m.group(1) if m else "?"
for f in sorted(glob.glob(sys.argv[1] + "/p1_*_deviation.csv")):
    sh = shard_of(f)
    for r in csv.DictReader(open(f)): dev[(sh, r["name"])] = r
for f in sorted(glob.glob(sys.argv[1] + "/p1_*_epitope.csv")):
    sh = shard_of(f)
    for r in csv.DictReader(open(f)): epi[(sh, r["name"])] = r
print(f"deviation rows {len(dev)}, epitope rows {len(epi)}, joined {len(set(dev)&set(epi))}")

HOT = [24, 25, 27, 42]
rows = []
for key in sorted(set(dev) & set(epi)):
    sh, n = key
    d, e = dev[key], epi[key]
    try:
        rows.append({"name": n, "shard": sh, "len": int(d["binder_len"]),
                     "motif": float(d["motif_rmsd_vs_af3"]),
                     "rec": float(d["receptor_rmsd"]),
                     "hits": int(e["n_hotspots_hit"]),
                     "mind": float(e["min_hotspot_dist"]) if e["min_hotspot_dist"] else 99.0,
                     **{f"A{h}": float(e.get(f"d_des_A{h}", 99)) for h in HOT}})
    except (ValueError, KeyError):
        continue
if not rows: print("nothing joined"); sys.exit(1)

m = [r["motif"] for r in rows]; rc = [r["rec"] for r in rows]; h = [r["hits"] for r in rows]
print(f"\nmotif RMSD    : min {min(m):.2f} med {st.median(m):.2f} max {max(m):.2f}   under 2 A: {sum(1 for x in m if x<2)}/{len(m)}")
print(f"receptor RMSD : min {min(rc):.2f} med {st.median(rc):.2f} max {max(rc):.2f}")
print(f"hotspots hit  : median {st.median(h)}  max {max(h)}  (of {len(HOT)})")
for k in range(len(HOT), -1, -1):
    c = sum(1 for x in h if x == k)
    if c: print(f"   {k}/{len(HOT)} hotspots: {c:4d} designs")
print(f"\nper-hotspot contact rate (<=5 A from the designed region):")
for hh in HOT:
    v = [r[f"A{hh}"] for r in rows]
    print(f"   A{hh}: {sum(1 for x in v if x<=5):4d}/{len(v)}  median {st.median(v):.1f} A")

# candidates: all 4 hotspots AND motif held AND receptor frame trustworthy
cand = [r for r in rows if r["hits"] == len(HOT) and r["motif"] < 1.5 and r["rec"] < 5.0]
cand.sort(key=lambda r: (r["mind"], r["motif"]))
print(f"\n=== candidates: {len(HOT)}/{len(HOT)} hotspots, motif <1.5 A, receptor <5 A ===")
print(f"found {len(cand)} of {len(rows)}\n")
print("%-40s %3s %4s %7s %7s %8s" % ("design", "sh", "len", "motif", "recRMSD", "minDist"))
for r in cand[:20]:
    print("%-40s %3s %4d %7.2f %7.2f %8.2f" % (r["name"].replace("glp1r_ecd4_glp1r_", "")[:40],
                                               r["shard"], r["len"], r["motif"], r["rec"], r["mind"]))
json.dump({"n_scored": len(rows), "n_candidates": len(cand),
           "motif_median": st.median(m), "hotspots": HOT,
           "candidates": cand[:60]}, open(sys.argv[1] + "/p1_summary.json", "w"), indent=1)
print(f"\nwrote {sys.argv[1]}/p1_summary.json")
