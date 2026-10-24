"""Let the hotspot count be stated explicitly.

best_dir infers the number of hotspots by counting commas in --hotspots, which
works when a campaign has one set but not for the p14 round, whose backbones come
from sources with 3 or 4 hotspots. With a prose description it inferred 1 and then
reported "0/50 of these contact all 1 hotspots", which is simply false.
"""
p = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts/best_dir.py"
s = open(p).read()
if "--hotspot-count" in s:
    print("already patched"); raise SystemExit(0)
old = 'ap.add_argument("--contig", default="")'
new = ('ap.add_argument("--contig", default="")\n'
       'ap.add_argument("--hotspot-count", type=int, default=None,\n'
       '                help="number of hotspots to call full engagement. Needed when "\n'
       '                     "--hotspots is prose or the round mixes sources with "\n'
       '                     "different set sizes, where counting commas is wrong.")')
assert s.count(old) == 1
s = s.replace(old, new)
old = 'nhot = len(a.hotspots.split(","))'
new = 'nhot = a.hotspot_count or len(a.hotspots.split(","))'
assert s.count(old) == 1
s = s.replace(old, new)
open(p, "w").write(s)
import ast; ast.parse(open(p).read())
print("patched best_dir: --hotspot-count override")
