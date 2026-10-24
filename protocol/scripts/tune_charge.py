"""Tune the charge to land near +5/+6 instead of +12.

Measured on shard 0 of 1091666 with --bulk-charge-bias 1.0:
    q(11-140) median +12, min +3, max +26
Two problems: the centre is twice the +6 target, and the spread is far too wide.

Fixes:
  1. bulk-charge-bias 1.0 -> 0.65. The baseline across earlier campaigns was -5.7,
     so bias 1.0 produced a swing of +17.7. Reaching +6 needs about +11.7, i.e.
     roughly two thirds of that effect.
  2. A FLAT-BOTTOMED charge window in the selection score. Previously the penalty
     was linear from the target, so a candidate at +20 with two extra salt bridges
     still won. Now anything inside target +- window is penalty-free and everything
     outside is penalised hard, so selection cannot trade charge away for contacts.
"""
import json, glob, os
from collections import Counter

p = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts/h6_mpnn_charge.py"
s = open(p).read()

if "charge_window" not in s:
    old = '''    ap.add_argument("--charge-weight", type=float, default=0.35,
                    help="penalty per unit of deviation from --target-charge")'''
    new = '''    ap.add_argument("--charge-weight", type=float, default=1.0,
                    help="penalty per unit of deviation beyond --charge-window")
    ap.add_argument("--charge-window", type=float, default=2.0,
                    help="deviation from --target-charge that costs nothing; beyond "
                         "this the penalty is steep, so selection cannot buy salt "
                         "bridges by overshooting the charge")'''
    assert s.count(old) == 1, s.count(old)
    s = s.replace(old, new)

    old = '''                          "score": sb + 0.4 * iface_q + 0.2 * hyd
                                   - a.charge_weight * abs(des_q - a.target_charge)})'''
    new = '''                          "score": sb + 0.4 * iface_q + 0.2 * hyd
                                   - a.charge_weight
                                     * max(0.0, abs(des_q - a.target_charge)
                                                - a.charge_window)})'''
    assert s.count(old) == 1, s.count(old)
    s = s.replace(old, new)
    open(p, "w").write(s)
    import ast; ast.parse(open(p).read())
    print("patched h6: flat-bottomed charge window in the selection score")
else:
    print("h6 already has the charge window")

# --- verify the scoring maths offline, no job needed ---------------------
def score(sb, iface_q, hyd, des_q, target=6.0, w=1.0, win=2.0):
    return sb + 0.4 * iface_q + 0.2 * hyd - w * max(0.0, abs(des_q - target) - win)

print("\nselection check: a high-charge candidate must no longer beat an on-target one")
cases = [("on target  q=+6, sb=6", score(6, 3, 20, 6)),
         ("on target  q=+8, sb=6", score(6, 3, 20, 8)),
         ("overshoot  q=+12, sb=8", score(8, 5, 22, 12)),
         ("overshoot  q=+20, sb=10", score(10, 6, 24, 20)),
         ("undershoot q=-2, sb=7", score(7, 2, 21, -2))]
for lbl, v in sorted(cases, key=lambda x: -x[1]):
    print(f"  {lbl:26s} score {v:7.2f}")
best = max(cases, key=lambda x: x[1])
print(f"  winner: {best[0]}")
assert "on target" in best[0] or "q=+12" not in best[0], "window is not biting"

# --- report full-chain vs designed-region charge from the smoke ----------
print("\nwhat the two charge definitions mean, from smoke 1091599:")
for f in sorted(glob.glob("/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/"
                          "p14_1091599/sh0/mpnn/h6_mpnn.json")):
    for r in json.load(open(f))["designs"]:
        sq = r["seq"]
        q = lambda t: Counter(t)["K"] + Counter(t)["R"] - Counter(t)["D"] - Counter(t)["E"]
        print(f"  q(full 140)={q(sq):+3d}   q(11-140)={q(sq[10:]):+3d}   "
              f"q(anchor)={q(sq[:10]):+3d}   anchor={sq[:10]}")
print("  the anchor contributes about -2, so a designed-region target of +6")
print("  gives a full-chain net charge of roughly +4.")
