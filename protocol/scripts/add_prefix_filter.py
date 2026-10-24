p = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts/h5_epitope.py"
s = open(p).read()
if "--name-prefix" in s:
    print("already patched"); raise SystemExit(0)
import shutil, os
if not os.path.exists(p + ".bak3"): shutil.copy(p, p + ".bak3")

old_arg = '    ap.add_argument("--manifest", required=True)'
new_arg = (old_arg + '\n'
           '    # the p14 round mixes backbones from campaigns with DIFFERENT hotspot sets in\n'
           '    # one prediction directory, so each set must be scored over only its own\n'
           '    # designs. Without this the glob would score every design four times.\n'
           '    ap.add_argument("--name-prefix", default=None,\n'
           '                    help="only score designs whose directory name starts with this")')
assert s.count(old_arg) == 1, s.count(old_arg)
s = s.replace(old_arg, new_arg)

old_glob = '    for cif in sorted(glob.glob(os.path.join(a.pred_dir, "*", "*_model_0.cif"))):'
new_glob = ('    cifs = sorted(glob.glob(os.path.join(a.pred_dir, "*", "*_model_0.cif")))\n'
            '    if a.name_prefix:\n'
            '        before = len(cifs)\n'
            '        cifs = [c for c in cifs\n'
            '                if os.path.basename(os.path.dirname(c)).startswith(a.name_prefix)]\n'
            '        print(f"name-prefix {a.name_prefix!r}: {len(cifs)}/{before} designs")\n'
            '    for cif in cifs:')
assert s.count(old_glob) == 1, s.count(old_glob)
s = s.replace(old_glob, new_glob)
open(p, "w").write(s)
import ast; ast.parse(open(p).read())
print("added --name-prefix to h5_epitope.py (backup .bak3); file parses")
