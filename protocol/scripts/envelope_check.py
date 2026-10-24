#!/usr/bin/env python
"""Score the GLP-1R p15 designs on the four composition axes from the EGFR submission.

The envelope comes from the EGFR repository (analysis/composition_envelope.md) and is
the observed range of ~10 experimentally confirmed sub-micromolar all-alpha EGFR
binders. That document is explicit that it is a bounding box around a handful of
successes, not a fitted model -- and the companion analysis states outright that
composition filters must be derived from the target's own cohort and never imported.
So this is a yardstick, not a gate: it says where GLP-1R designs sit relative to a
region where a different target's binders are known to have worked.
"""
import csv, statistics as st, sys
from collections import Counter
from Bio.SeqUtils.ProtParam import ProteinAnalysis

BOUNDS = {"alanine": (None, 0.104, 0.085),
          "hydrophobic_AVILMFWY": (0.329, 0.42, 0.376),
          "net_charge_pH7": (-19.3, -5.0, None),
          "pI": (4.60, 5.00, None)}
HYD = set("AVILMFWY")

rows = list(csv.DictReader(open(sys.argv[1])))
region = sys.argv[2] if len(sys.argv) > 2 else "full"
vals = {k: [] for k in BOUNDS}
for r in rows:
    s = r["seq"] if region == "full" else r["seq"][10:]
    p = ProteinAnalysis(s); c = Counter(s); n = len(s)
    vals["alanine"].append(c["A"] / n)
    vals["hydrophobic_AVILMFWY"].append(sum(c[a] for a in HYD) / n)
    vals["net_charge_pH7"].append(p.charge_at_pH(7.0))
    vals["pI"].append(p.isoelectric_point())

print(f"{len(rows)} designs, {region} chain\n")
print(f"{'axis':24s} {'EGFR envelope':>22s} {'GLP-1R median':>14s} "
      f"{'range':>18s} {'inside':>8s}")
for k, (lo, hi, cen) in BOUNDS.items():
    v = sorted(vals[k])
    inside = sum(1 for x in v if (lo is None or x >= lo) and (hi is None or x <= hi))
    env = (f"<= {hi}" if lo is None else
           f"{lo} .. {hi}")
    print(f"{k:24s} {env:>22s} {st.median(v):14.3f} "
          f"{f'{min(v):.2f} .. {max(v):.2f}':>18s} {inside:4d}/{len(v)}")
    if cen is not None:
        print(f"{'':24s} {'target ' + str(cen):>22s} "
              f"{'(x%.1f the target)' % (st.median(v)/cen):>14s}")
allin = 0
for i in range(len(rows)):
    ok = all((BOUNDS[k][0] is None or vals[k][i] >= BOUNDS[k][0]) and
             (BOUNDS[k][1] is None or vals[k][i] <= BOUNDS[k][1]) for k in BOUNDS)
    allin += ok
print(f"\ninside all four axes: {allin}/{len(rows)}")
print("\nhighest alanine fraction reported in ANY observed binder (external figure"
      "\nquoted in the EGFR repo): 0.158")
over = sum(1 for x in vals["alanine"] if x > 0.158)
print(f"GLP-1R designs above that: {over}/{len(rows)}")
