#!/usr/bin/env python
"""Build the four-hotspot input (A27/A42/A86/A95) and the P13 campaign script.

A95:CD2 was chosen by burial, not convenience: at 30 neighbouring heavy atoms
within 8 A it is the most exposed atom in the set, against 41 for A86:CD2 (94%
engaged), 58 for A42:OG1 (96%) and 78 for A27:CD1 (80%, the weak one). It also
leaves the span unchanged at 16.5 A, sitting between the two existing sites.

A27:CD1 is deliberately left as-is so that A95 is the only change from P12.
"""
import json, os

G = "/lustre/home/gsgueglia/challenge2026/glp1r2026"
S = os.path.join(G, "scripts")
HOT = {"A27": "CD1", "A42": "OG1", "A86": "CD2", "A95": "CD2"}

src = os.path.join(G, "inputs/rfd3/w3_r27_42_86.json")
dst = os.path.join(G, "inputs/rfd3/w4_r27_42_86_95.json")
d = json.load(open(src))
out = {}
for k, v in d.items():
    nk = k.replace("w3_r27_42_86", "w4_r27_42_86_95")
    nv = dict(v)
    old = nv.get("select_hotspots")
    nv["select_hotspots"] = dict(HOT)
    out[nk] = nv
    print(f"  {k} -> {nk}")
    print(f"      {old}  ->  {nv['select_hotspots']}")
json.dump(out, open(dst, "w"), indent=2)
chk = json.load(open(dst))
assert len(chk) == len(d) == 5
for (ko, vo), (kn, vn) in zip(d.items(), chk.items()):
    for f in vo:
        if f == "select_hotspots": continue
        assert vo[f] == vn[f], (f, vo[f], vn[f])
    assert vn["select_hotspots"] == HOT
    assert vn["contig"] == "B1-10,130-130,/0,A1-409"
print(f"\nverified: only select_hotspots differs from P12; contig unchanged")
print(f"wrote {dst}\n")

s = open(os.path.join(S, "p12_prod.sbatch")).read()
subs = [
    ("#SBATCH --job-name=p12-prod", "#SBATCH --job-name=p13-prod"),
    ("/logs/p12_prod_%A_%a.out", "/logs/p13_prod_%A_%a.out"),
    ('runs/glp1r/p12_${RUN}/sh${SH}', 'runs/glp1r/p13_${RUN}/sh${SH}'),
    ('inputs/rfd3/w3_r27_42_86.json', 'inputs/rfd3/w4_r27_42_86_95.json'),
    ("seed=$((12000 + SH))", "seed=$((15000 + SH))"),
    ("--seed $((131 + SH))", "--seed $((191 + SH))"),
    ('"$DIFF" 27,42,86 30', '"$DIFF" 27,42,86,95 30'),
    ("--hotspots 27,42,86", "--hotspots 27,42,86,95"),
    ("SEGRG-p12 sh${SH}", "SEGRG-p13 sh${SH}"),
]
for old, new in subs:
    n = s.count(old)
    assert n >= 1, f"not found: {old!r}"
    s = s.replace(old, new)
    print(f"  {n}x {old[:48]!r} -> {new[:40]!r}")
n = s.count('p12_${RUN}_sh${SH}'); s = s.replace('p12_${RUN}_sh${SH}', 'p13_${RUN}_sh${SH}')
print(f"  {n}x analysis prefix p12_ -> p13_")
n2 = s.count('p12_${RUN}_outdirs.txt'); s = s.replace('p12_${RUN}_outdirs.txt', 'p13_${RUN}_outdirs.txt')
print(f"  {n2}x outdirs")
hdr = """#
# P13: P12's set plus Leu95. FOUR hotspots now, so ">=3 hits" means 3 of 4 again
# and the gate is directly comparable to P2 (64% at 4/4) and P3 (23% at 4/4),
# not to P12's 3-of-3.
#
# A95:CD2 picked by burial (scripts/pick_atom.py): 30 heavy-atom neighbours within
# 8 A, the most exposed atom in the set - cf. A86:CD2 41 (94% engaged), A42:OG1 58
# (96%), A27:CD1 78 (80%, the weak one). Span stays 16.5 A; A95 sits between the
# two existing sites (10.2 A from A42, 14.6 A from A86) and is the closest of the
# four to the anchor at 33.5 A.
#
# A27:CD1 is left unchanged so A95 is the only difference from P12.
"""
s = s.replace("#\n# PRODUCTION batch", hdr + "#\n# PRODUCTION batch", 1)
dst2 = os.path.join(S, "p13_prod.sbatch")
open(dst2, "w").write(s)
os.chmod(dst2, 0o755)
for stale in ('p12_${RUN}', "p12-prod", "w3_r27_42_86.json", "--hotspots 27,42,86\n"):
    assert stale not in s, f"stale: {stale!r}"
print("\nno stale p12 references; wrote", dst2)
