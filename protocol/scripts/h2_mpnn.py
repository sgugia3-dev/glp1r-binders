#!/usr/bin/env python
"""Stage 2 -- ProteinMPNN on the clean RFD3 backbones.

Two requirements, both enforced rather than nudged:

 1. Binder position 3 must be ASP or GLU. A logit `bias` cannot guarantee that, so the
    other 18 amino acids are passed in `omit_per_residue` for that position, and the
    result is ASSERTED afterwards.
 2. The 10 fixed peptide-motif residues must survive unchanged. `--fixed_residues` is
    mutually exclusive with the chain selectors, so every DESIGNABLE position is
    enumerated in `designed_residues` instead, and the motif is verified byte-identical
    to the input afterwards.

Motif location differs per design: it is read from each design's own
`diffused_index_map` (B<i> -> A<j>), never assumed.
"""
import argparse, glob, gzip, json, os, shutil, sys

ACIDS = ["ASP", "GLU"]
ALL20 = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
         "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
THREE2ONE = dict(zip(ALL20, "ARNDCQEGHILKMFPSTWYV"))


def binder_seq_from_cif(path, chain="A"):
    """1-based {resnum: one-letter} for one chain of an mmCIF."""
    hdr, rows, inl = [], [], False
    op = gzip.open if path.endswith(".gz") else open
    for line in op(path, "rt"):
        if line.startswith("_atom_site."):
            hdr.append(line.strip().split(".")[1]); inl = True; continue
        if inl and line.startswith(("ATOM", "HETATM")): rows.append(line.split())
        elif inl and line.startswith("#"): break
    i = {h: k for k, h in enumerate(hdr)}
    ca, cr, cn = i["auth_asym_id"], i["auth_seq_id"], i["label_comp_id"]
    out = {}
    for r in rows:
        if len(r) <= max(ca, cr, cn) or r[ca] != chain: continue
        try: n = int(r[cr])
        except ValueError: continue
        out[n] = THREE2ONE.get(r[cn], "X")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, help="RFD3 output dir")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--weights", required=True)
    ap.add_argument("--nseq", type=int, default=16)
    ap.add_argument("--temperature", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=37)
    ap.add_argument("--acid-pos", type=int, default=3, help="binder position forced to D/E")
    ap.add_argument("--max-ca-dev", type=float, default=2.0)
    ap.add_argument("--max-breaks", type=int, default=1)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    work = os.path.join(a.outdir, "_cif"); os.makedirs(work, exist_ok=True)

    # ---- free screen: keep only backbones that actually held the motif
    recs, skipped = [], []
    for j in sorted(glob.glob(os.path.join(a.indir, "*.json"))):
        d = json.load(open(j)); m = d["metrics"]
        tag = os.path.basename(j)[:-5]
        if m["max_ca_deviation"] > a.max_ca_dev or m["n_chainbreaks"] > a.max_breaks:
            skipped.append((tag, round(m["max_ca_deviation"], 2), m["n_chainbreaks"])); continue
        gz = os.path.join(a.indir, tag + ".cif.gz")
        if not os.path.exists(gz):
            skipped.append((tag, "no cif", 0)); continue
        cif = os.path.join(work, tag + ".cif")
        with gzip.open(gz, "rt") as fi, open(cif, "w") as fo: shutil.copyfileobj(fi, fo)
        motif = sorted(int(v[1:]) for k, v in d["diffused_index_map"].items()
                       if k.startswith("B") and v.startswith("A"))
        recs.append({"tag": tag, "cif": cif, "motif": motif})

    print(f"backbones in   : {len(recs) + len(skipped)}")
    print(f"passed screen  : {len(recs)}  (max_ca_dev<={a.max_ca_dev}, breaks<={a.max_breaks})")
    for t, dv, br in skipped: print(f"  skipped {t}: dev={dv} breaks={br}")
    if not recs:
        print("FATAL: no backbone passed the screen -- nothing to design")
        return 3

    from mpnn.inference_engines.mpnn import MPNNInferenceEngine
    engine = MPNNInferenceEngine(
        model_type="protein_mpnn", checkpoint_path=a.weights, is_legacy_weights=True,
        out_directory=a.outdir, write_fasta=False, write_structures=False, device="cpu")

    rows, fails = [], 0
    for rec in recs:
        seq_in = binder_seq_from_cif(rec["cif"], "A")
        blen = max(seq_in)
        motif = set(rec["motif"])
        designed = [f"A{k}" for k in range(1, blen + 1) if k not in motif]
        # With C-terminal-only extension the motif occupies binder positions 1-10, so the
        # requested acid position can fall INSIDE the fixed motif. That is not an error:
        # the position simply cannot be designed, so the constraint has to be satisfied by
        # the native residue instead. Verify it rather than silently forcing nothing.
        if a.acid_pos in motif:
            native = seq_in.get(a.acid_pos)
            acid_by_motif = native in ("D", "E")
            print(f"  note: position {a.acid_pos} is inside the fixed motif "
                  f"{min(motif)}-{max(motif)}; native residue is {native} "
                  f"({'already an acid' if acid_by_motif else 'NOT an acid'})")
            if not acid_by_motif:
                print(f"FATAL {rec['tag']}: position {a.acid_pos} is fixed as {native}, "
                      f"so the D/E requirement cannot be met without moving the motif")
                return 4
            omit_pr, bias_pr = None, None
        else:
            # force D/E at the requested position by omitting everything else there
            omit_pr = {f"A{a.acid_pos}": [x for x in ALL20 if x not in ACIDS]}
            bias_pr = {f"A{a.acid_pos}": {"ASP": 2.0, "GLU": 2.0}}
        inp = {"bias": None, "bias_per_residue": bias_pr,
               "omit": ["CYS"], "omit_per_residue": omit_pr,
               "structure_path": rec["cif"], "name": rec["tag"],
               "seed": a.seed, "batch_size": a.nseq, "number_of_batches": 1,
               "temperature": a.temperature, "designed_residues": designed}
        outs = engine.run(input_dicts=[inp], atom_arrays=None)
        # engine.run returns a LIST of mpnn.utils.inference.MPNNInferenceOutput, one per
        # sampled sequence, each carrying .output_dict["designed_sequence"] as a single
        # STRING for the whole complex (binder first, chain A). Verified by probe 1088560.
        got = 0
        for o in (outs if isinstance(outs, list) else [outs]):
            od = getattr(o, "output_dict", None)
            if not isinstance(od, dict): continue
            raw = od.get("designed_sequence")
            if isinstance(raw, (list, tuple)): raw = raw[0] if raw else None
            if not isinstance(raw, str) or not raw: continue
            seq = raw.replace(":", "")
            if True:
                binder = seq[:blen]
                aa3 = binder[a.acid_pos - 1]
                acid_src = "motif" if a.acid_pos in motif else "designed"
                motif_ok = all(binder[k - 1] == seq_in[k] for k in sorted(motif))
                if aa3 not in ("D", "E"): fails += 1
                rows.append({"tag": rec["tag"], "binder_len": blen,
                             "motif": sorted(motif), "seq": binder,
                             "pos%d" % a.acid_pos: aa3,
                             "acid_ok": aa3 in ("D", "E"), "acid_source": acid_src,
                             "motif_ok": motif_ok})
                got += 1
        print(f"  {rec['tag']}: {got} seqs, len {blen}, motif {min(motif)}-{max(motif)}")

    nacid = sum(1 for r in rows if r["acid_ok"])
    nmot = sum(1 for r in rows if r["motif_ok"])
    print(f"\nsequences      : {len(rows)}")
    print(f"pos{a.acid_pos} is D/E   : {nacid}/{len(rows)}")
    print(f"motif preserved: {nmot}/{len(rows)}")
    json.dump({"n": len(rows), "acid_pos": a.acid_pos,
               "n_acid_ok": nacid, "n_motif_ok": nmot,
               "skipped": skipped, "designs": rows},
              open(os.path.join(a.outdir, "h2_mpnn.json"), "w"), indent=1)
    if not rows: print("FATAL: engine returned no sequences"); return 5
    if nacid != len(rows):
        print(f"FATAL: {len(rows)-nacid} sequences violate the D/E constraint"); return 6
    if nmot != len(rows):
        print(f"FATAL: {len(rows)-nmot} sequences altered the fixed motif"); return 7
    print("wrote " + os.path.join(a.outdir, "h2_mpnn.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
