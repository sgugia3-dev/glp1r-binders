p = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts/best_dir.py"
s = open(p).read()
old = 'n4 = sum(1 for r in kept if int(r["hotspots"]) == 4)\nnhot = len(a.hotspots.split(","))'
new = ('nhot = len(a.hotspots.split(","))\n'
       '# count full engagement against THIS campaign\'s hotspot count, not a literal 4:\n'
       '# P10 has three hotspots, so "== 4" would report zero fully-engaged designs.\n'
       'n4 = sum(1 for r in kept if int(r["hotspots"]) == nhot)')
if new.split("\n")[-1] in s:
    print("already fixed"); raise SystemExit(0)
assert s.count(old) == 1, ("pattern", s.count(old))
open(p, "w").write(s.replace(old, new))
import ast; ast.parse(open(p).read())
print("fixed: full-engagement count now uses the campaign's own hotspot count")
