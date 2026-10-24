"""Close the loop on charge instead of guessing a bias strength.

Calibration 1091958 at --bulk-charge-bias 0.85 gave q(full) median +11 against a
+5.5 target, and per backbone it was hopeless to tune with one number:

    backbone   exposed bulk positions   best q(full)   in-window
    p10_01     40                       +7             7/32
    p13_01     44                       +7            15/32
    p12_01     47                      +11             0/32
    p3_01      69                      +18             0/32

Charge scales with how many positions get swapped, and that count is a property of
the backbone's own surface, so no single bias strength can serve all of them.

So the bias is now fitted per backbone. Sample once, measure the achieved median,
and solve for the strength that would have hit the target, assuming charge responds
linearly in the bias weight:

    q(w) ~= Q0 + k*w ,  Q0 = unbiased baseline, k fitted from the first round

Q0 is the measured pre-redesign baseline: designed region -5.7 and native anchor
-2, so -7.7 for the full chain. Up to --charge-rounds attempts, stopping as soon as
the median lands inside the window.
"""
p = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts/h6_mpnn_charge.py"
s = open(p).read()
if "charge_rounds" in s:
    print("already patched"); raise SystemExit(0)

old = '''    ap.add_argument("--burial-cut", type=float, default=85.0,'''
new = '''    ap.add_argument("--charge-rounds", type=int, default=3,
                    help="refit the bulk bias per backbone up to this many times")
    ap.add_argument("--baseline-charge", type=float, default=-7.7,
                    help="measured unbiased full-chain net charge, used to fit the "
                         "bias strength from one sampling round")
    ap.add_argument("--burial-cut", type=float, default=85.0,'''
assert s.count(old) == 1
s = s.replace(old, new)

# wrap the sampling in a fitting loop ------------------------------------
old = '''        inp = {"bias": None, "bias_per_residue": bias_pr,
               "omit": ["CYS"], "omit_per_residue": omit_pr,
               "structure_path": rec["cif"], "name": rec["tag"],
               "seed": a.seed, "batch_size": a.nseq, "number_of_batches": 1,
               "temperature": a.temperature, "designed_residues": designed}
        outs = engine.run(input_dicts=[inp], atom_arrays=None)

        cands = []'''
new = '''        bulk_pos = [k for k in range(ANCHOR_HI + 1, blen + 1)
                    if k not in iface_pos and f"{ch}{k}" not in bias_pr
                    and bmap.get(k, 999) < a.burial_cut]
        w_bulk = a.bulk_charge_bias
        cands = []
        for rnd in range(max(1, a.charge_rounds)):
            for k in bulk_pos:
                bias_pr[f"{ch}{k}"] = {"ASP": -w_bulk, "GLU": -w_bulk,
                                       "LYS": w_bulk, "ARG": w_bulk}
            inp = {"bias": None, "bias_per_residue": bias_pr,
                   "omit": ["CYS"], "omit_per_residue": omit_pr,
                   "structure_path": rec["cif"], "name": rec["tag"],
                   "seed": a.seed + 97 * rnd, "batch_size": a.nseq,
                   "number_of_batches": 1, "temperature": a.temperature,
                   "designed_residues": designed}
            outs = engine.run(input_dicts=[inp], atom_arrays=None)
            cands = []'''
assert s.count(old) == 1
s = s.replace(old, new)

# reindent the candidate-collection block by 4 spaces --------------------
start = s.index('            cands = []\n', s.index('for rnd in range'))
end = s.index('        passing = list(cands)')
block = s[start + len('            cands = []\n'):end]
reindented = "".join(("    " + l if l.strip() else l) for l in block.splitlines(True))
loop_tail = '''            qsel = sorted((c["full_q"] if a.charge_basis == "full" else c["des_q"])
                          for c in cands)
            if not qsel: break
            med = qsel[len(qsel) // 2]
            if abs(med - a.target_charge) <= a.charge_window or rnd == a.charge_rounds - 1:
                if rnd:
                    print(f"    round {rnd+1}: w={w_bulk:.2f} -> median "
                          f"q({a.charge_basis})={med:+.0f}  (converged)")
                break
            # fit k from this round and solve for the weight that hits the target
            denom = (med - a.baseline_charge)
            k = denom / w_bulk if w_bulk > 1e-6 else 0.0
            w_new = ((a.target_charge - a.baseline_charge) / k) if abs(k) > 1e-6 else w_bulk
            w_new = max(0.05, min(3.0, w_new))
            print(f"    round {rnd+1}: w={w_bulk:.2f} gave median {med:+.0f}, "
                  f"refitting to w={w_new:.2f} ({len(bulk_pos)} bulk positions)")
            w_bulk = w_new
'''
s = s[:start + len('            cands = []\n')] + reindented + loop_tail + s[end:]

# the old unconditional bulk-bias block is now dead; remove it ----------
old_dead = '''        bmap = rec["bmap"]
        n_bulk = 0
        for k in range(ANCHOR_HI + 1, blen + 1):
            if k in iface_pos or f"{ch}{k}" in bias_pr: continue
            if bmap.get(k, 999) >= a.burial_cut: continue
            w2 = a.bulk_charge_bias
            bias_pr[f"{ch}{k}"] = {"ASP": -w2, "GLU": -w2, "LYS": w2, "ARG": w2}
            n_bulk += 1
        print(f"  {rec['tag']}: {len(iface_pos)} interface positions biased "
              f"({dict(kinds)}), {n_bulk} exposed bulk positions charge-swapped")'''
new_live = '''        bmap = rec["bmap"]
        print(f"  {rec['tag']}: {len(iface_pos)} interface positions biased "
              f"({dict(kinds)})")'''
assert s.count(old_dead) == 1
s = s.replace(old_dead, new_live)

open(p, "w").write(s)
import ast; ast.parse(open(p).read())
print("patched h6: per-backbone closed-loop charge fitting")

# sanity-check the fitting maths
Q0, target = -7.7, 5.5
for n, med, w in ((40, 9, 0.85), (47, 11, 0.85), (69, 18, 0.85), (44, 7, 0.85)):
    k = (med - Q0) / w
    wn = (target - Q0) / k
    print(f"  {n:3d} bulk pos: median {med:+3d} at w={w:.2f} -> refit w={wn:.2f}")
