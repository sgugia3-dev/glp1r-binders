p = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts/h6_mpnn_charge.py"
s = open(p).read()
old = '    allxyz = [(x, y, z) for ats in at.values() for _, x, y, z in ats]'
new = ('    # at is {chain: {resnum: [(name, x, y, z)]}} - two levels of dict before the\n'
       '    # atom list, so both have to be traversed.\n'
       '    allxyz = [(x, y, z) for chd in at.values() for ats in chd.values()\n'
       '              for _, x, y, z in ats]')
assert s.count(old) == 1, s.count(old)
s = s.replace(old, new)
open(p, "w").write(s)
import ast; ast.parse(open(p).read())
print("fixed burial_map traversal")

# self-test the helper against a real structure before resubmitting
import glob, sys
sys.argv = ["x"]
ns = {}
exec(compile(s.split("def main(")[0], "h6", "exec"), ns)
f = sorted(glob.glob("/lustre/home/gsgueglia/challenge2026/glp1r2026/p12_best/structures/*.cif"))[0]
at, seq = ns["read_cif"](f)
bnd = sorted(at, key=lambda c: -len(at[c]))[1]
bm = ns["burial_map"](at, bnd)
v = sorted(bm.values())
print(f"self-test on {f.split('/')[-1]}: {len(bm)} positions, "
      f"burial min {min(v)} p25 {v[len(v)//4]} median {v[len(v)//2]} max {max(v)}")
exposed = sum(1 for x in v if x < 85)
print(f"  exposed (<85): {exposed}/{len(v)}")
assert len(bm) == 140 and 0 < min(v) < max(v) < 400, "burial values look wrong"
print("  helper verified")
