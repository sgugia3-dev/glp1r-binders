p = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts/campaign_agg.py"
s = open(p).read()
old = 'print(f"{RUN}  combined hotspots A27 + A42 + A86   ({n} refolds)")'
new = ('print(f"{RUN}   hotspots "\n'
       '      + " + ".join(c.replace("d_des_", "") for c in sorted(hot_cols))\n'
       '      + f"   ({n} refolds)")')
if "sorted(hot_cols))" in s:
    print("already generic"); raise SystemExit(0)
assert s.count(old) == 1, s.count(old)
open(p, "w").write(s.replace(old, new))
import ast; ast.parse(open(p).read())
print("header now lists the campaign's actual hotspots, not a hardcoded trio")
