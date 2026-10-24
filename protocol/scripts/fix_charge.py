"""Strengthen the charge modulation in h6, and make selection a window not a filter.

Smoke test 1091026 held both anchor constraints 3/3 and gave 7 salt bridges, but
only 3 of 48 sampled sequences reached net charge >= 0 and the median landed at
just +1. The cause: only ~36 of 130 designed positions touch the receptor and get
biased, while the other ~94 keep MPNN's default composition, which the campaign-wide
analysis measured at 16-19% glutamate. The bulk drags the net charge negative.

Two changes:

1. EXPOSED NON-INTERFACE POSITIONS now get a charge swap (D/E down, K/R up).
   Buried positions are left alone so the hydrophobic core survives. Burial is the
   count of heavy atoms within 8 A of CB (CA for Gly); measured on these backbones
   the distribution is p25=75, median=94, max=117, so the default threshold of 85
   treats roughly the exposed 40% as swappable.

2. CHARGE IS A TARGET WINDOW IN THE SCORE, not a hard accept/reject. The hard
   filter discarded 45 of 48 sequences and then fell back to keeping whatever had
   the least-negative charge, which wasted the sampling. Scoring the distance to a
   target charge keeps every candidate usable and still prefers the right ones.
"""
import os, re

p = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts/h6_mpnn_charge.py"
s = open(p).read()
if "burial" in s:
    print("already patched"); raise SystemExit(0)
if not os.path.exists(p + ".bak"):
    import shutil; shutil.copy(p, p + ".bak")

# --- burial helper ---------------------------------------------------------
anchor = "def interface_map("
burial_fn = '''def burial_map(at, bnd):
    """{binder resnum: heavy atoms within 8 A of CB (CA for Gly)}.

    A charge swap is safe on a solvent-exposed position and damaging on a buried
    one, so the bulk-charge bias is applied only to the exposed set.
    """
    allxyz = [(x, y, z) for ats in at.values() for _, x, y, z in ats]
    ref = {}
    for rn, ats in at[bnd].items():
        cb = next((t for t in ats if t[0] == "CB"), None)
        ca = next((t for t in ats if t[0] == "CA"), None)
        t = cb or ca
        if t: ref[rn] = t[1:]
    out = {}
    for rn, p0 in ref.items():
        out[rn] = sum(1 for q in allxyz
                      if (p0[0]-q[0])**2 + (p0[1]-q[1])**2 + (p0[2]-q[2])**2 < 64.0)
    return out


'''
assert s.count(anchor) == 1
s = s.replace(anchor, burial_fn + anchor, 1)

# --- new args --------------------------------------------------------------
old = '    ap.add_argument("--min-designed-charge", type=float, default=0.0,\n' \
      '                    help="reject sequences whose 11-end net charge is below this")'
new = '''    ap.add_argument("--target-charge", type=float, default=6.0,
                    help="preferred net charge of residues 11-end; scored, not enforced")
    ap.add_argument("--charge-weight", type=float, default=0.35,
                    help="penalty per unit of deviation from --target-charge")
    ap.add_argument("--bulk-charge-bias", type=float, default=1.0,
                    help="D/E down and K/R up at exposed non-interface positions")
    ap.add_argument("--burial-cut", type=float, default=85.0,
                    help="positions with burial below this count as exposed")'''
assert s.count(old) == 1
s = s.replace(old, new)

# --- compute burial alongside the interface map ----------------------------
old = '        imap = interface_map(at, seq, bnd, tgt)\n'
new = '        imap = interface_map(at, seq, bnd, tgt)\n        bmap = burial_map(at, bnd)\n'
assert s.count(old) == 1
s = s.replace(old, new)
old = '                     "bnd": bnd, "tgt": tgt, "blen": blen, "imap": imap,'
new = '                     "bnd": bnd, "tgt": tgt, "blen": blen, "imap": imap,\n                     "bmap": bmap,'
assert s.count(old) == 1
s = s.replace(old, new)

# --- apply the bulk bias ---------------------------------------------------
old = '''        print(f"  {rec['tag']}: {len(iface_pos)} interface positions biased "
              f"({dict(kinds)})")'''
new = '''        # exposed positions that do NOT touch the receptor: swap charge safely.
        # This is what actually moves the net charge; the interface set is far too
        # small (~36 of 130) to overcome MPNN's glutamate-heavy bulk composition.
        bmap = rec["bmap"]
        n_bulk = 0
        for k in range(ANCHOR_HI + 1, blen + 1):
            if k in iface_pos or f"{ch}{k}" in bias_pr: continue
            if bmap.get(k, 999) >= a.burial_cut: continue
            w2 = a.bulk_charge_bias
            bias_pr[f"{ch}{k}"] = {"ASP": -w2, "GLU": -w2, "LYS": w2, "ARG": w2}
            n_bulk += 1
        print(f"  {rec['tag']}: {len(iface_pos)} interface positions biased "
              f"({dict(kinds)}), {n_bulk} exposed bulk positions charge-swapped")'''
assert s.count(old) == 1
s = s.replace(old, new)

# --- selection: window, not filter ----------------------------------------
old = '''            cands.append({"seq": binder, "des_q": des_q, "iface_q": iface_q,
                          "sb_pot": sb, "hyd_pot": hyd,
                          "score": sb + 0.4 * iface_q + 0.2 * hyd})'''
new = '''            cands.append({"seq": binder, "des_q": des_q, "iface_q": iface_q,
                          "sb_pot": sb, "hyd_pot": hyd,
                          "score": sb + 0.4 * iface_q + 0.2 * hyd
                                   - a.charge_weight * abs(des_q - a.target_charge)})'''
assert s.count(old) == 1
s = s.replace(old, new)

old = '''        passing = [c for c in cands if c["des_q"] >= a.min_designed_charge]
        if not passing:
            best_q = max(c["des_q"] for c in cands)
            print(f"    note {rec['tag']}: none reached net charge "
                  f">= {a.min_designed_charge:g} (best {best_q:+.0f}); keeping the top by charge")
            passing = sorted(cands, key=lambda c: -c["des_q"])[: a.keep]
            stats["charge_target_missed"] += 1
        passing.sort(key=lambda c: -c["score"])'''
new = '''        passing = list(cands)
        n_pos = sum(1 for c in cands if c["des_q"] > 0)
        if n_pos == 0:
            print(f"    WARNING {rec['tag']}: no sampled sequence reached a positive "
                  f"net charge (best {max(c['des_q'] for c in cands):+.0f})")
            stats["no_positive_charge"] += 1
        passing.sort(key=lambda c: -c["score"])'''
assert s.count(old) == 1
s = s.replace(old, new)

old = '''        print(f"    kept {min(len(passing), a.keep)}/{len(cands)}  "
              f"best: q(11-140)={passing[0]['des_q']:+.0f} "
              f"iface_q={passing[0]['iface_q']:+.0f} sb={passing[0]['sb_pot']}")'''
new = '''        print(f"    kept {min(len(passing), a.keep)}/{len(cands)} "
              f"({n_pos}/{len(cands)} sampled were net positive)  "
              f"best: q(11-140)={passing[0]['des_q']:+.0f} "
              f"iface_q={passing[0]['iface_q']:+.0f} sb={passing[0]['sb_pot']}")'''
assert s.count(old) == 1
s = s.replace(old, new)

open(p, "w").write(s)
import ast; ast.parse(open(p).read())
print("patched h6: burial-aware bulk charge swap + charge-window selection")
