"""Let h7 also read a flat directory of complexes.

The packaged best sets store structures as <dir>/structures/*.cif, while Boltz
output nests one directory per prediction. Measuring the pre-redesign baseline
with the same code as the post-redesign round is the only way the before/after
interaction counts are comparable, so h7 needs to read both layouts.
"""
p = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts/h7_interactions.py"
s = open(p).read()
if "--flat" in s:
    print("already patched"); raise SystemExit(0)
old = 'ap.add_argument("--out-csv", required=True)'
new = ('ap.add_argument("--out-csv", required=True)\n'
       'ap.add_argument("--flat", action="store_true",\n'
       '                help="pred-dir holds *.cif directly instead of one "\n'
       '                     "subdirectory per prediction")')
assert s.count(old) == 1
s = s.replace(old, new)

old = '''dirs = [d for d in sorted(glob.glob(os.path.join(a.pred_dir, "*"))) if os.path.isdir(d)]
rows = []
for d in dirs:
    cif = glob.glob(os.path.join(d, "*_model_0.cif"))
    if not cif: continue
    r = analyse(cif[0])
    if r is None: continue
    rec = {"name": os.path.basename(d)}'''
new = '''if a.flat:
    targets = [(os.path.basename(c).replace("_complex.cif", "").replace(".cif", ""), c)
               for c in sorted(glob.glob(os.path.join(a.pred_dir, "*.cif")))]
else:
    targets = []
    for d in sorted(glob.glob(os.path.join(a.pred_dir, "*"))):
        if not os.path.isdir(d): continue
        cif = glob.glob(os.path.join(d, "*_model_0.cif"))
        if cif: targets.append((os.path.basename(d), cif[0]))
rows = []
for nm, cifpath in targets:
    r = analyse(cifpath)
    if r is None: continue
    rec = {"name": nm}'''
assert s.count(old) == 1
s = s.replace(old, new)
open(p, "w").write(s)
import ast; ast.parse(open(p).read())
print("patched h7: added --flat")
