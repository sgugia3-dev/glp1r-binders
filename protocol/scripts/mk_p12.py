"""Production batch for the combined A27/A42/A86 set, derived from p10_combo.sbatch.

Sized to P2's scale (20 shards x 20 designs x 3 folds = 1200 folds) so the hit
rate is directly comparable to P2's 1200 and P3's 480 rather than to the 90-fold
smoke test. Only the run tag, wall time and deadline change.
"""
import os
S = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts"
s = open(os.path.join(S, "p10_combo.sbatch")).read()
subs = [
    ("#SBATCH --job-name=p10-combo", "#SBATCH --job-name=p12-prod"),
    ("/logs/p10_combo_%A_%a.out", "/logs/p12_prod_%A_%a.out"),
    ("#SBATCH --time=05:00:00", "#SBATCH --time=04:00:00"),
    ("DEADLINE_S=${DEADLINE_S:-15300}", "DEADLINE_S=${DEADLINE_S:-12600}"),
    ('runs/glp1r/p10_${RUN}/sh${SH}', 'runs/glp1r/p12_${RUN}/sh${SH}'),
    ("seed=$((9000 + SH))", "seed=$((12000 + SH))"),
    ("--seed $((71 + SH))", "--seed $((131 + SH))"),
    ("SEGRG-combo sh${SH}", "SEGRG-p12 sh${SH}"),
]
for old, new in subs:
    n = s.count(old)
    assert n >= 1, f"not found: {old!r}"
    s = s.replace(old, new)
    print(f"  {n}x {old[:52]!r} -> {new[:38]!r}")
n = s.count('p10_${RUN}_sh${SH}')
s = s.replace('p10_${RUN}_sh${SH}', 'p12_${RUN}_sh${SH}')
print(f"  {n}x analysis prefix p10_ -> p12_")
n2 = s.count('p10_${RUN}_outdirs.txt')
s = s.replace('p10_${RUN}_outdirs.txt', 'p12_${RUN}_outdirs.txt')
print(f"  {n2}x outdirs")
hdr = """#
# PRODUCTION batch for the combined hotspot set, confirming the smoke test.
# Smoke test (1089071, 90 folds) gave 82% full 3/3 engagement, 0 complete misses,
# and A86 engaged in 93% of designs despite sitting 16.5 A from A27. This runs the
# same config at P2's scale so the numbers are comparable, not smoke-test sized.
#
# Sizing: 20 shards is BOTH limits at once - the gpus QOS allows MaxJobsPU=20, and
# only 20 of the 37 free GPUs are reachable because nodes do not oversubscribe and
# wn07/11/16/18/23 have free GPUs but zero free CPUs.
# 20 shards x 20 designs x 3 folds/design = 400 designs, 1200 folds.
"""
s = s.replace("#\n# SMOKE TEST:", hdr + "#\n# (smoke-test note retained below)\n# SMOKE TEST:", 1)
dst = os.path.join(S, "p12_prod.sbatch")
open(dst, "w").write(s)
os.chmod(dst, 0o755)
for stale in ('p10_${RUN}', "p10-combo", "seed=$((9000"):
    assert stale not in s, f"stale: {stale}"
print("\nno stale p10 run references; wrote", dst)
