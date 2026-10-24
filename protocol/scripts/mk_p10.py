#!/usr/bin/env python
"""Derive the combined-hotspot smoke test from p3_all.sbatch by exact substitution.

Changing only what must change keeps the run comparable to P2 and P3: the input
JSON, the hotspot list passed to the two measurements, the run tag and the labels.
Every replacement is asserted, so a drifted template fails loudly instead of
silently producing a differently-configured run.
"""
import os, re

S = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts"
src = open(os.path.join(S, "p3_all.sbatch")).read()
subs = [
    ("#SBATCH --job-name=p3-glp-r87", "#SBATCH --job-name=p10-combo"),
    ("/logs/p3_all_%A_%a.out", "/logs/p10_combo_%A_%a.out"),
    ('OUT="$C26_SCRATCH/runs/glp1r/p3_${RUN}/sh${SH}"',
     'OUT="$C26_SCRATCH/runs/glp1r/p10_${RUN}/sh${SH}"'),
    ('inputs="$G/inputs/rfd3/w2_r86_89.json"',
     'inputs="$G/inputs/rfd3/w3_r27_42_86.json"'),
    ('seed=$((7000 + SH))', 'seed=$((9000 + SH))'),
    ('--seed $((37 + SH))', '--seed $((71 + SH))'),
    ('"$G/scripts/r_seg_measure.py" "$DIFF" 86,87,88,89 30',
     '"$G/scripts/r_seg_measure.py" "$DIFF" 27,42,86 30'),
    ("--hotspots 86,87,88,89", "--hotspots 27,42,86"),
    ('SEGRG-r87 sh${SH}', 'SEGRG-combo sh${SH}'),
]
out = src
for old, new in subs:
    n = out.count(old)
    assert n >= 1, f"template drifted, not found: {old!r}"
    out = out.replace(old, new)
    print(f"  {n}x  {old[:58]!r:62s} -> {new[:40]!r}")

# p3_${RUN}_sh${SH} appears in several analysis paths; rename them all
n = len(re.findall(r'p3_\$\{RUN\}_sh\$\{SH\}', out))
out = out.replace('p3_${RUN}_sh${SH}', 'p10_${RUN}_sh${SH}')
print(f"  {n}x  analysis file prefix p3_ -> p10_")
n2 = out.count('p3_${RUN}_outdirs.txt')
out = out.replace('p3_${RUN}_outdirs.txt', 'p10_${RUN}_outdirs.txt')
print(f"  {n2}x  outdirs file")

hdr = """#
# SMOKE TEST: the P2/P3 protocol exactly, with hotspots A27 + A42 + A86 - one
# patch drawn from each face. A27:CD1 and A42:OG1 come from P2, A86:CD2 from P3,
# so the atoms are the ones those campaigns actually steered on.
#
# Geometry measured before launching (scripts/hot_geom.py):
#   A27-A42  6.1 A, A42-A86 11.9 A, A27-A86 16.5 A -> max span 16.5 A
#   references: P2 24/25/27/42 span 11.6 A, P3 86-89 span 8.5 A
#   anchor to A27 37.4 A, to A42 37.2 A, to A86 47.1 A
# So A27+A42 are effectively one site and A86 a second site ~12-16 A away. The
# open question is whether one binder bridges both, or trades one for the other.
#
# Everything else - contig B1-10,130-130,/0,A1-409, infer_ori_strategy=hotspots,
# is_non_loopy, step_scale, gamma_0, MPNN settings, Boltz settings - is identical
# to p3_all.sbatch, so the hit rate is directly comparable to P2's 64% at 4/4
# and P3's 23% at 4/4 (note: 3 hotspots here, not 4, so compare the FRACTION).
"""
out = out.replace("#\n# All three stages in ONE job", hdr + "#\n# All three stages in ONE job", 1)
dst = os.path.join(S, "p10_combo.sbatch")
open(dst, "w").write(out)
os.chmod(dst, 0o755)
# nothing from the old run may survive
for stale in ("w2_r86_89", "86,87,88,89", "p3_${RUN}", "p3-glp-r87"):
    assert stale not in out, f"stale reference left behind: {stale}"
print("\nno stale P3 references remain")
print(f"wrote {dst}")
