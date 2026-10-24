#!/usr/bin/env python
"""Gather the ranked candidates into one flat, archived directory.

Names collide across shards (rfd3 derives them from the input JSON), so every file is
renamed `rank##_sh<shard>_<name>` to make it unique and self-describing. Also emits a
single table with the metrics AND the binder sequence, which is what is actually needed
downstream.
"""
import csv, glob, json, os, shutil, sys

analysis, run_root, out = sys.argv[1], sys.argv[2], sys.argv[3]
summ = json.load(open(os.path.join(analysis, "p1_summary.json")))
cands = summ["candidates"]
print(f"candidates in summary: {len(cands)}")

# binder sequences live in each shard's manifest
seqs = {}
for man_p in glob.glob(os.path.join(run_root, "sh*", "manifest.json")):
    sh = os.path.basename(os.path.dirname(man_p)).replace("sh", "")
    for y in json.load(open(man_p))["yamls"]:
        seqs[(sh, y["name"])] = y

os.makedirs(out, exist_ok=True)
rows, missing = [], 0
for i, c in enumerate(cands, 1):
    sh, name = str(c["shard"]), c["name"]
    hits = glob.glob(os.path.join(run_root, f"sh{sh}", "**", "predictions", name,
                                  f"{name}_model_0.cif"), recursive=True)
    if not hits:
        missing += 1; continue
    stem = f"rank{i:02d}_sh{sh}_{name}"
    shutil.copy2(hits[0], os.path.join(out, stem + ".cif"))
    conf = glob.glob(os.path.join(os.path.dirname(hits[0]), "confidence_*.json"))
    cj = {}
    if conf:
        shutil.copy2(conf[0], os.path.join(out, stem + "_confidence.json"))
        cj = json.load(open(conf[0]))
    y = seqs.get((sh, name), {})
    rows.append({"rank": i, "file": stem + ".cif", "shard": sh, "design": name,
                 "binder_len": c["len"], "motif_rmsd": c["motif"],
                 "receptor_rmsd": c["rec"], "hotspots_hit": c["hits"],
                 "min_hotspot_dist": c["mind"],
                 "iptm": round(cj.get("iptm", -1), 3), "ptm": round(cj.get("ptm", -1), 3),
                 "complex_plddt": round(cj.get("complex_plddt", -1), 3),
                 "binder_seq": y.get("binder_seq", "")})

if not rows: print("FATAL: collected nothing"); sys.exit(3)
cols = list(rows[0].keys())
tbl = os.path.join(out, "TOP_CANDIDATES.csv")
with open(tbl, "w") as fh:
    fh.write(",".join(cols) + "\n")
    for r in rows: fh.write(",".join(str(r[c]) for c in cols) + "\n")
print(f"collected {len(rows)} structures ({missing} missing) -> {out}")
print(f"table: {tbl}\n")
print("%-5s %-26s %4s %7s %8s %7s" % ("rank", "design", "len", "motif", "minDist", "ipTM"))
for r in rows[:10]:
    print("%-5d %-26s %4d %7.2f %8.2f %7.3f" % (r["rank"], r["design"].replace("glp1r_ecd4_glp1r_", "")[:26],
          r["binder_len"], r["motif_rmsd"], r["min_hotspot_dist"], r["iptm"]))
import statistics as st
print(f"\nipTM over the set: min {min(r['iptm'] for r in rows):.3f} "
      f"med {st.median([r['iptm'] for r in rows]):.3f} max {max(r['iptm'] for r in rows):.3f}")
