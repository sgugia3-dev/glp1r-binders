import os
S = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts"
# add P13 to the head-to-head table
p = os.path.join(S, "campaign_agg.py")
s = open(p).read()
old = '                     ("p12_1089085", "P12 A27/42/86 prod", 3)):'
new = ('                     ("p12_1089085", "P12 A27/42/86 prod", 3),\n'
       '                     ("p13_1089119", "P13 A27/42/86/95", 4)):')
if "p13_1089119" not in s:
    assert s.count(old) == 1, s.count(old)
    open(p, "w").write(s.replace(old, new))
    print("added P13 to the head-to-head comparison")
else:
    print("P13 already in the comparison")

# finish chain
f = open(os.path.join(S, "p12_finish.sbatch")).read()
for a, b in (("p12_1089085", "p13_1089119"), ("p12_RANKED", "p13_RANKED"),
             ("glp1r2026/p12_best", "glp1r2026/p13_best"), ("--prefix p12", "--prefix p13"),
             ("--campaign p12_1089085", "--campaign p13_1089119"),
             ("-J p12fin", "-J p13fin")):
    f = f.replace(a, b)
f = f.replace('--hotspots "A27:CD1, A42:OG1, A86:CD2"',
              '--hotspots "A27:CD1, A42:OG1, A86:CD2, A95:CD2"')
assert "A95:CD2" in f and "p13_1089119" in f and "p12" not in f.replace("p12_finish", "")
dst = os.path.join(S, "p13_finish.sbatch")
open(dst, "w").write(f); os.chmod(dst, 0o755)
print("wrote", dst)
