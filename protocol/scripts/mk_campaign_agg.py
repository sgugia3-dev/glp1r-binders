"""Generalize p10_agg.py so the campaign run id is an argument."""
import os
S = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts"
s = open(os.path.join(S, "p10_agg.py")).read()
old_a = "A = sys.argv[1]"
new_a = 'A = sys.argv[1]\nRUN = sys.argv[2] if len(sys.argv) > 2 else "p10_1089071"'
old_l = 'dev, epi = load("p10_1089071", "deviation"), load("p10_1089071", "epitope")'
new_l = 'dev, epi = load(RUN, "deviation"), load(RUN, "epitope")'
old_h = '''for run, lbl, nh in (("p2_1088883", "P2  A24/25/27/42", 4),
                     ("p3_1088996", "P3  A86-89", 4),
                     ("p10_1089071", "P10 A27/42/86", 3)):'''
new_h = '''for run, lbl, nh in (("p2_1088883", "P2  A24/25/27/42", 4),
                     ("p3_1088996", "P3  A86-89", 4),
                     ("p10_1089071", "P10 A27/42/86 smoke", 3),
                     ("p12_1089085", "P12 A27/42/86 prod", 3)):'''
for o, n in ((old_a, new_a), (old_l, new_l), (old_h, new_h)):
    assert s.count(o) == 1, (o[:40], s.count(o))
    s = s.replace(o, n)
s = s.replace('print(f"P10 combined hotspots A27 + A42 + A86   ({n} refolds)")',
              'print(f"{RUN}  combined hotspots A27 + A42 + A86   ({n} refolds)")')
s = s.replace('print("\\nNOTE: P10 has 3 hotspots, P2 and P3 have 4, so \'full\' is an easier bar")\nprint("for P10. The per-hotspot and bridging numbers above are the fair comparison.")',
              'print("\\nNOTE: the A27/42/86 runs have 3 hotspots, P2 and P3 have 4, so \'full\'")\nprint("is an easier bar for them. The per-hotspot and bridging numbers are the")\nprint("fair comparison.")')
dst = os.path.join(S, "campaign_agg.py")
open(dst, "w").write(s)
import ast; ast.parse(s)
print("wrote", dst)
