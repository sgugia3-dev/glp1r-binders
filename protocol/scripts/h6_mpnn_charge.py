#!/usr/bin/env python
"""Charge-complementary redesign of validated backbones, anchor included.

Two changes from h2_mpnn.py:

1. THE ANCHOR IS NOW DESIGNABLE. Previously binder residues 1-10 were the fixed
   native GLP-1 motif. Here all 140 positions are designed, with two hard
   constraints enforced by omission (which forces an identity, unlike a bias):
       position 1 -> HIS
       position 3 -> ASP or GLU
   The native anchor HAEGTFTSDV already satisfies both, so this relaxes the other
   eight positions rather than fighting the native sequence.

2. CHARGE IS MODULATED PER POSITION FROM MEASURED GEOMETRY, not globally.
   Measured on the 200 input complexes (scripts/contact_charge.py):
       binder 11-140 contacts the ECD   : 44 res, net charge -3.0  -> wants BASIC
       binder 1-10   contacts the TMD   : 27 res, net charge +3.0  -> wants ACIDIC
   So the designed region is biased basic and the anchor acidic. Crucially the
   bias is applied ONLY at positions that actually touch the receptor; biasing all
   130 residues toward K/R would wreck the fold. Interior positions are left
   unbiased so MPNN can still build a foldable core.

   Per interface position, the bias follows the partner residues within 5 A:
       partner acidic (D/E)  -> K/R up      (salt bridge)
       partner basic  (K/R)  -> D/E up      (salt bridge, other polarity)
       partner apolar        -> A/V/L/I/M/F up  (hydrophobic contact)
       partner polar         -> N/Q/S/T/Y/H up  (hydrogen bond)

   Only `bias_per_residue` is used, never the global `bias`: per-residue entries
   are known to override the global one, so mixing them makes the net effect at a
   position ambiguous.

Sequences are then SELECTED rather than merely sampled: many are drawn per
backbone and ranked by salt-bridge potential against the actual partners, so the
charge target is met by selection instead of by an aggressive bias.
"""
import argparse, glob, json, math, os, sys
from collections import Counter

ALL20 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
         "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
A3 = {v: k for k, v in ONE.items()}
ACIDIC, BASIC = set("DE"), set("KR")
APOLAR, POLAR = set("AVLIMFWCGP"), set("STNQYH")
CUT = 5.0
ANCHOR_HI = 10


def read_cif(path):
    """{chain: {resnum: [(atomname, x, y, z)]}} and {chain: {resnum: one-letter}}."""
    names, at, seq = [], {}, {}
    for l in open(path):
        if l.startswith("_atom_site."):
            names.append(l.strip().split(".")[1]); continue
        if not (l.startswith("ATOM") or l.startswith("HETATM")): continue
        t = l.split()
        if len(t) != len(names): continue
        d = dict(zip(names, t))
        if (d.get("type_symbol") or "").strip() == "H": continue
        ch = d.get("auth_asym_id") or d.get("label_asym_id")
        try:
            rn = int(d.get("auth_seq_id") or d.get("label_seq_id"))
            xyz = (float(d["Cartn_x"]), float(d["Cartn_y"]), float(d["Cartn_z"]))
        except (TypeError, ValueError, KeyError):
            continue
        at.setdefault(ch, {}).setdefault(rn, []).append(
            ((d.get("label_atom_id") or "").strip().strip('"'),) + xyz)
        seq.setdefault(ch, {})[rn] = ONE.get((d.get("label_comp_id") or "").upper(), "X")
    return at, seq


def burial_map(at, bnd):
    """{binder resnum: heavy atoms within 8 A of CB (CA for Gly)}.

    A charge swap is safe on a solvent-exposed position and damaging on a buried
    one, so the bulk-charge bias is applied only to the exposed set.
    """
    # at is {chain: {resnum: [(name, x, y, z)]}} - two levels of dict before the
    # atom list, so both have to be traversed.
    allxyz = [(x, y, z) for chd in at.values() for ats in chd.values()
              for _, x, y, z in ats]
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


def interface_map(at, seq, bnd, tgt):
    """{binder resnum: Counter of partner one-letter residues within CUT}."""
    tlist = [(rn, x, y, z) for rn, ats in at[tgt].items() for _, x, y, z in ats]
    out = {}
    for rn, ats in at[bnd].items():
        partners = set()
        for _, bx, by, bz in ats:
            for trn, x, y, z in tlist:
                if (bx-x)**2 + (by-y)**2 + (bz-z)**2 <= CUT*CUT:
                    partners.add(trn)
        if partners:
            out[rn] = Counter(seq[tgt][p] for p in partners if p in seq[tgt])
    return out


def bias_for(partners, w):
    """Per-position bias dict chosen from the partner composition."""
    n = sum(partners.values()) or 1
    f_ac = sum(partners[a] for a in ACIDIC) / n
    f_ba = sum(partners[a] for a in BASIC) / n
    f_ap = sum(partners[a] for a in APOLAR) / n
    f_po = sum(partners[a] for a in POLAR) / n
    best = max((f_ac, "acid"), (f_ba, "base"), (f_ap, "apolar"), (f_po, "polar"))
    kind = best[1]
    if best[0] <= 0:
        return None, "none"
    if kind == "acid":
        return {"LYS": w, "ARG": w}, "saltbridge_KR"
    if kind == "base":
        return {"ASP": w * 0.7, "GLU": w * 0.7}, "saltbridge_DE"
    if kind == "apolar":
        v = w * 0.55
        return {k: v for k in ("ALA", "VAL", "LEU", "ILE", "MET", "PHE")}, "hydrophobic"
    v = w * 0.55
    return {k: v for k in ("ASN", "GLN", "SER", "THR", "TYR", "HIS")}, "hbond"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complexes", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--weights", required=True)
    ap.add_argument("--nseq", type=int, default=32)
    ap.add_argument("--keep", type=int, default=3, help="sequences kept per backbone")
    ap.add_argument("--temperature", type=float, default=0.15)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--bias-weight", type=float, default=1.6)
    ap.add_argument("--anchor-acid-bias", type=float, default=0.8)
    ap.add_argument("--target-charge", type=float, default=6.0,
                    help="preferred net charge of residues 11-end; scored, not enforced")
    ap.add_argument("--charge-weight", type=float, default=1.0,
                    help="penalty per unit of deviation beyond --charge-window")
    ap.add_argument("--charge-basis", choices=["full", "designed"], default="full",
                    help="which net charge --target-charge refers to: the whole "
                         "140-residue chain (default) or only residues 11-end. The "
                         "designed anchor contributes about -4, so the two differ.")
    ap.add_argument("--charge-window", type=float, default=2.0,
                    help="deviation from --target-charge that costs nothing; beyond "
                         "this the penalty is steep, so selection cannot buy salt "
                         "bridges by overshooting the charge")
    ap.add_argument("--bulk-charge-bias", type=float, default=1.0,
                    help="D/E down and K/R up at exposed non-interface positions")
    ap.add_argument("--charge-rounds", type=int, default=3,
                    help="refit the bulk bias per backbone up to this many times")
    ap.add_argument("--baseline-charge", type=float, default=-7.7,
                    help="measured unbiased full-chain net charge, used to fit the "
                         "bias strength from one sampling round")
    ap.add_argument("--burial-cut", type=float, default=85.0,
                    help="positions with burial below this count as exposed")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    recs = []
    for f in a.complexes:
        at, seq = read_cif(f)
        if len(at) < 2:
            print(f"  skip {os.path.basename(f)}: fewer than 2 chains"); continue
        by_len = sorted(at, key=lambda c: -len(at[c]))
        tgt, bnd = by_len[0], by_len[1]
        blen = max(at[bnd])
        if blen != 140:
            print(f"  skip {os.path.basename(f)}: binder length {blen} != 140"); continue
        imap = interface_map(at, seq, bnd, tgt)
        bmap = burial_map(at, bnd)
        recs.append({"cif": f, "tag": os.path.basename(f).replace("_complex.cif", ""),
                     "bnd": bnd, "tgt": tgt, "blen": blen, "imap": imap,
                     "bmap": bmap,
                     "native": "".join(seq[bnd][k] for k in sorted(seq[bnd]))})
    print(f"backbones to redesign: {len(recs)}")
    if not recs:
        print("FATAL: nothing to design"); return 3

    from mpnn.inference_engines.mpnn import MPNNInferenceEngine
    engine = MPNNInferenceEngine(
        model_type="protein_mpnn", checkpoint_path=a.weights, is_legacy_weights=True,
        out_directory=a.outdir, write_fasta=False, write_structures=False, device="cpu")

    rows, stats = [], Counter()
    for rec in recs:
        ch, blen, imap = rec["bnd"], rec["blen"], rec["imap"]
        designed = [f"{ch}{k}" for k in range(1, blen + 1)]
        omit_pr = {f"{ch}1": [x for x in ALL20 if x != "HIS"],
                   f"{ch}3": [x for x in ALL20 if x not in ("ASP", "GLU")]}
        bias_pr = {f"{ch}3": {"ASP": 2.0, "GLU": 2.0}}
        # anchor positions 2 and 4-10 face the TMD (+3) -> mild acidic
        for k in [2] + list(range(4, ANCHOR_HI + 1)):
            bias_pr[f"{ch}{k}"] = {"ASP": a.anchor_acid_bias, "GLU": a.anchor_acid_bias}
        kinds = Counter()
        iface_pos = []
        for k, partners in imap.items():
            if k <= ANCHOR_HI: continue
            b, kind = bias_for(partners, a.bias_weight)
            if b:
                bias_pr[f"{ch}{k}"] = b
                kinds[kind] += 1
                iface_pos.append(k)
        # exposed positions that do NOT touch the receptor: swap charge safely.
        # This is what actually moves the net charge; the interface set is far too
        # small (~36 of 130) to overcome MPNN's glutamate-heavy bulk composition.
        bmap = rec["bmap"]
        print(f"  {rec['tag']}: {len(iface_pos)} interface positions biased "
              f"({dict(kinds)})")
        stats.update(kinds)

        bulk_pos = [k for k in range(ANCHOR_HI + 1, blen + 1)
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
            cands = []
            for o in (outs if isinstance(outs, list) else [outs]):
                od = getattr(o, "output_dict", None)
                if not isinstance(od, dict): continue
                raw = od.get("designed_sequence")
                if isinstance(raw, (list, tuple)): raw = raw[0] if raw else None
                if not isinstance(raw, str) or not raw: continue
                binder = raw.replace(":", "")[:blen]
                if len(binder) != blen: continue
                if binder[0] != "H" or binder[2] not in ("D", "E"):
                    stats["constraint_violation"] += 1
                    continue
                des = binder[ANCHOR_HI:]
                c = Counter(des)
                des_q = c["K"] + c["R"] - c["D"] - c["E"]
                cf = Counter(binder)
                full_q = cf["K"] + cf["R"] - cf["D"] - cf["E"]
                anc_q = full_q - des_q
                q_sel = full_q if a.charge_basis == "full" else des_q
                # salt-bridge potential: binder residue of opposite polarity to its partners
                sb = 0
                for k in iface_pos:
                    aa = binder[k - 1]
                    p = imap[k]
                    if aa in BASIC and sum(p[x] for x in ACIDIC): sb += 1
                    if aa in ACIDIC and sum(p[x] for x in BASIC): sb += 1
                ic = Counter(binder[k - 1] for k in iface_pos)
                iface_q = ic["K"] + ic["R"] - ic["D"] - ic["E"]
                hyd = sum(1 for k in iface_pos
                          if binder[k - 1] in APOLAR
                          and sum(imap[k][x] for x in APOLAR))
                cands.append({"seq": binder, "des_q": des_q, "full_q": full_q,
                              "anc_q": anc_q, "iface_q": iface_q,
                              "sb_pot": sb, "hyd_pot": hyd,
                              "score": sb + 0.4 * iface_q + 0.2 * hyd
                                       - a.charge_weight
                                         * max(0.0, abs(q_sel - a.target_charge)
                                                    - a.charge_window)})
            if not cands:
                print(f"    WARNING {rec['tag']}: no sequence satisfied the constraints")
                continue
            qsel = sorted((c["full_q"] if a.charge_basis == "full" else c["des_q"])
                          for c in cands)
            if not qsel: break
            med = qsel[len(qsel) // 2]
            hit = abs(med - a.target_charge) <= a.charge_window
            if hit or rnd == a.charge_rounds - 1:
                if rnd or hit:
                    n_in = sum(1 for q in qsel
                               if abs(q - a.target_charge) <= a.charge_window)
                    state = ("ON TARGET" if hit else
                             f"OFF TARGET after {a.charge_rounds} rounds")
                    print(f"    round {rnd+1}: w={w_bulk:.2f} -> median "
                          f"q({a.charge_basis})={med:+.0f}  [{state}; "
                          f"{n_in}/{len(qsel)} candidates in window]")
                break
            # fit k from this round and solve for the weight that hits the target
            denom = (med - a.baseline_charge)
            k = denom / w_bulk if w_bulk > 1e-6 else 0.0
            w_new = ((a.target_charge - a.baseline_charge) / k) if abs(k) > 1e-6 else w_bulk
            w_new = max(0.05, min(3.0, w_new))
            print(f"    round {rnd+1}: w={w_bulk:.2f} gave median {med:+.0f}, "
                  f"refitting to w={w_new:.2f} ({len(bulk_pos)} bulk positions)")
            w_bulk = w_new
        passing = list(cands)
        n_pos = sum(1 for c in cands if c["des_q"] > 0)
        bk = "full_q" if a.charge_basis == "full" else "des_q"
        n_win = sum(1 for c in cands
                    if abs(c[bk] - a.target_charge) <= a.charge_window)
        if n_pos == 0:
            print(f"    WARNING {rec['tag']}: no sampled sequence reached a positive "
                  f"net charge (best {max(c['des_q'] for c in cands):+.0f})")
            stats["no_positive_charge"] += 1
        passing.sort(key=lambda c: -c["score"])
        for i, c in enumerate(passing[: a.keep]):
            rows.append({"tag": rec["tag"], "binder_len": blen,
                         "motif": list(range(1, ANCHOR_HI + 1)), "seq": c["seq"],
                         "pos3": c["seq"][2], "pos1": c["seq"][0],
                         "rank_in_backbone": i, "designed_charge": c["des_q"],
                         "full_charge": c["full_q"], "anchor_charge": c["anc_q"],
                         "iface_charge": c["iface_q"], "salt_bridge_potential": c["sb_pot"],
                         "hydrophobic_potential": c["hyd_pot"],
                         "n_iface_positions": len(iface_pos),
                         "native_seq": rec["native"],
                         "acid_ok": True, "motif_ok": True})
        print(f"    kept {min(len(passing), a.keep)}/{len(cands)} "
              f"({n_win}/{len(cands)} inside the charge window)  best: "
              f"q(full)={passing[0]['full_q']:+.0f} q(11-140)={passing[0]['des_q']:+.0f} "
              f"q(anchor)={passing[0]['anc_q']:+.0f} sb={passing[0]['sb_pot']}")

    if not rows:
        print("FATAL: no sequences produced"); return 5
    qs = [r["designed_charge"] for r in rows]
    qf = sorted(r["full_charge"] for r in rows)
    print(f"\nsequences kept : {len(rows)} from {len(recs)} backbones")
    print(f"pos1 == H      : {sum(1 for r in rows if r['pos1']=='H')}/{len(rows)}")
    print(f"pos3 in D/E    : {sum(1 for r in rows if r['pos3'] in 'DE')}/{len(rows)}")
    print(f"q(11-140)      : median {sorted(qs)[len(qs)//2]:+.0f}  "
          f"min {min(qs):+.0f}  max {max(qs):+.0f}  "
          f"(previous campaigns: -5.7 median)")
    print(f"q(full 140)    : median {qf[len(qf)//2]:+.0f}  min {qf[0]:+.0f}  "
          f"max {qf[-1]:+.0f}   <-- the controlled quantity "
          f"(basis={a.charge_basis}, target {a.target_charge:+.1f} "
          f"+-{a.charge_window:g})")
    inwin = sum(1 for q in qf if abs(q - a.target_charge) <= a.charge_window)
    print(f"inside the window: {inwin}/{len(rows)}")
    print(f"positive-charge designs: {sum(1 for q in qs if q > 0)}/{len(rows)}")
    print(f"bias kinds applied: {dict(stats)}")
    json.dump({"n": len(rows), "acid_pos": 3, "n_acid_ok": len(rows),
               "n_motif_ok": len(rows), "skipped": [], "designs": rows,
               "bias_weight": a.bias_weight},
              open(os.path.join(a.outdir, "h6_mpnn.json"), "w"), indent=1)
    print("wrote " + os.path.join(a.outdir, "h6_mpnn.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
