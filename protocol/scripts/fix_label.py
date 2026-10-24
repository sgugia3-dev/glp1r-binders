"""Stop the loop claiming convergence when it merely ran out of rounds.

In 1092117 the p12 backbone printed "round 3: w=0.79 -> median q(full)=+19
(converged)" when +19 is nowhere near the +5.5 +-1.5 window. The condition fired
on `or rnd == charge_rounds - 1`, so the last round was always labelled converged.

The final designs were still correct, because SELECTION does the real work: with 32
samples there is enough spread that the flat-bottomed window finds in-window
candidates even when the round median is off (p12 had 5/32 in-window and all three
kept designs came out at q(full)=+7). The loop only improves the odds of finding
them. But a log that says "converged" when it did not is worse than no log.
"""
p = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts/h6_mpnn_charge.py"
s = open(p).read()
old = '''            if abs(med - a.target_charge) <= a.charge_window or rnd == a.charge_rounds - 1:
                if rnd:
                    print(f"    round {rnd+1}: w={w_bulk:.2f} -> median "
                          f"q({a.charge_basis})={med:+.0f}  (converged)")
                break'''
new = '''            hit = abs(med - a.target_charge) <= a.charge_window
            if hit or rnd == a.charge_rounds - 1:
                if rnd or hit:
                    n_in = sum(1 for q in qsel
                               if abs(q - a.target_charge) <= a.charge_window)
                    state = ("ON TARGET" if hit else
                             f"OFF TARGET after {a.charge_rounds} rounds")
                    print(f"    round {rnd+1}: w={w_bulk:.2f} -> median "
                          f"q({a.charge_basis})={med:+.0f}  [{state}; "
                          f"{n_in}/{len(qsel)} candidates in window]")
                break'''
assert s.count(old) == 1, s.count(old)
open(p, "w").write(s.replace(old, new))
import ast; ast.parse(open(p).read())
print("fixed: the log now distinguishes ON TARGET from OFF TARGET after N rounds")
