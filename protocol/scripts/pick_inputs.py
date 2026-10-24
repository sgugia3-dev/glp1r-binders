#!/usr/bin/env python
"""Select the best low-charge binders from p14 as inputs for the anchor repair.

Criteria: in the q(11-140) +4..+9 band, gated, and contacting every hotspot of
their source lineage. Deduplicated by backbone, best score first, so each input
is a distinct fold rather than sibling MPNN sequences of one.
"""
import csv, glob, os, re, sys
A, ROOT, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
rows = [r for r in csv.DictReader(open(os.path.join(A, "p14_band_q4to9.csv")))
        if r["gated"] == "True" and r["full"] == "True" and r["score"]]
rows.sort(key=lambda r: -float(r["score"]))
seen, out = set(), []
for r in rows:
    if r["bb_tag"] in seen: continue
    d = os.path.join(ROOT, "sh" + r["shard"], "boltz_results_yaml",
                     "predictions", r["name"])
    cif = glob.glob(os.path.join(d, "*_model_0.cif"))
    if not cif: continue
    seen.add(r["bb_tag"]); out.append((r, cif[0]))
with open(OUT, "w") as fh:
    for r, c in out: fh.write(c + "\n")
from collections import Counter
print(f"{len(rows)} gated+full band rows -> {len(out)} distinct backbones")
print(f"  source mix: {dict(Counter(r['src'] for r, _ in out))}")
print(f"  score range {float(out[-1][0]['score']):.3f} .. {float(out[0][0]['score']):.3f}")
print(f"  q(11-140) range {min(int(r['q_des']) for r,_ in out)} .. "
      f"{max(int(r['q_des']) for r,_ in out)}")
print(f"wrote {OUT}")
