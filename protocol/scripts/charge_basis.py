"""Make the controlled quantity the FULL-CHAIN net charge, as asked.

The anchor is now designed and comes out at -3 to -5 (His1 plus the forced acid at
position 3 plus the mild acidic bias on 2 and 4-10, all of which is correct for the
basic TMD surface it touches). So the designed region and the whole chain are two
different targets, about 4 units apart:

    smoke 1091599:  q(full 140) +5..+11     q(11-140) +10..+14     q(anchor) -3..-5

"Total charge" names the whole molecule, so selection now scores q(full 140) by
default and the target is +5.5 with a +-1.5 window, i.e. +4 to +7. The designed
region will sit near +9-10 as a consequence, which is what keeps the ECD-facing
surface complementary while the anchor stays acidic for the TMD.

--charge-basis designed restores the old behaviour if the designed region is the
quantity of interest instead.
"""
p = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts/h6_mpnn_charge.py"
s = open(p).read()

if "charge_basis" not in s:
    old = '''    ap.add_argument("--charge-window", type=float, default=2.0,'''
    new = '''    ap.add_argument("--charge-basis", choices=["full", "designed"], default="full",
                    help="which net charge --target-charge refers to: the whole "
                         "140-residue chain (default) or only residues 11-end. The "
                         "designed anchor contributes about -4, so the two differ.")
    ap.add_argument("--charge-window", type=float, default=2.0,'''
    assert s.count(old) == 1
    s = s.replace(old, new)

    # compute the full-chain charge next to the designed-region one
    old = '''            des = binder[ANCHOR_HI:]
            c = Counter(des)
            des_q = c["K"] + c["R"] - c["D"] - c["E"]'''
    new = '''            des = binder[ANCHOR_HI:]
            c = Counter(des)
            des_q = c["K"] + c["R"] - c["D"] - c["E"]
            cf = Counter(binder)
            full_q = cf["K"] + cf["R"] - cf["D"] - cf["E"]
            anc_q = full_q - des_q
            q_sel = full_q if a.charge_basis == "full" else des_q'''
    assert s.count(old) == 1
    s = s.replace(old, new)

    old = '''            cands.append({"seq": binder, "des_q": des_q, "iface_q": iface_q,
                          "sb_pot": sb, "hyd_pot": hyd,'''
    new = '''            cands.append({"seq": binder, "des_q": des_q, "full_q": full_q,
                          "anc_q": anc_q, "iface_q": iface_q,
                          "sb_pot": sb, "hyd_pot": hyd,'''
    assert s.count(old) == 1
    s = s.replace(old, new)

    old = '''                                   - a.charge_weight
                                     * max(0.0, abs(des_q - a.target_charge)
                                                - a.charge_window)})'''
    new = '''                                   - a.charge_weight
                                     * max(0.0, abs(q_sel - a.target_charge)
                                                - a.charge_window)})'''
    assert s.count(old) == 1
    s = s.replace(old, new)

    old = '''        n_pos = sum(1 for c in cands if c["des_q"] > 0)'''
    new = '''        n_pos = sum(1 for c in cands if c["des_q"] > 0)
        bk = "full_q" if a.charge_basis == "full" else "des_q"
        n_win = sum(1 for c in cands
                    if abs(c[bk] - a.target_charge) <= a.charge_window)'''
    assert s.count(old) == 1
    s = s.replace(old, new)

    old = '''        print(f"    kept {min(len(passing), a.keep)}/{len(cands)} "
              f"({n_pos}/{len(cands)} sampled were net positive)  "
              f"best: q(11-140)={passing[0]['des_q']:+.0f} "
              f"iface_q={passing[0]['iface_q']:+.0f} sb={passing[0]['sb_pot']}")'''
    new = '''        print(f"    kept {min(len(passing), a.keep)}/{len(cands)} "
              f"({n_win}/{len(cands)} inside the charge window)  best: "
              f"q(full)={passing[0]['full_q']:+.0f} q(11-140)={passing[0]['des_q']:+.0f} "
              f"q(anchor)={passing[0]['anc_q']:+.0f} sb={passing[0]['sb_pot']}")'''
    assert s.count(old) == 1
    s = s.replace(old, new)

    old = '''                         "rank_in_backbone": i, "designed_charge": c["des_q"],'''
    new = '''                         "rank_in_backbone": i, "designed_charge": c["des_q"],
                         "full_charge": c["full_q"], "anchor_charge": c["anc_q"],'''
    assert s.count(old) == 1
    s = s.replace(old, new)

    old = '''    qs = [r["designed_charge"] for r in rows]'''
    new = '''    qs = [r["designed_charge"] for r in rows]
    qf = sorted(r["full_charge"] for r in rows)'''
    assert s.count(old) == 1
    s = s.replace(old, new)

    old = '''    print(f"positive-charge designs: {sum(1 for q in qs if q > 0)}/{len(rows)}")'''
    new = '''    print(f"q(full 140)    : median {qf[len(qf)//2]:+.0f}  min {qf[0]:+.0f}  "
          f"max {qf[-1]:+.0f}   <-- the controlled quantity "
          f"(basis={a.charge_basis}, target {a.target_charge:+.1f} "
          f"+-{a.charge_window:g})")
    inwin = sum(1 for q in qf if abs(q - a.target_charge) <= a.charge_window)
    print(f"inside the window: {inwin}/{len(rows)}")
    print(f"positive-charge designs: {sum(1 for q in qs if q > 0)}/{len(rows)}")'''
    assert s.count(old) == 1
    s = s.replace(old, new)

    open(p, "w").write(s)
    import ast; ast.parse(open(p).read())
    print("patched h6: --charge-basis full|designed, full-chain charge reported")
else:
    print("already patched")
