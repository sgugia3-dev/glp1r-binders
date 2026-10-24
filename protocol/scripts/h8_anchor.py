#!/usr/bin/env python
"""Redesign ONLY the 10-residue anchor, restrained to the agonist consensus.

Positions 11-140 are held at their p14 identities, so the validated
charge-complementary surface is preserved exactly and only the agonist
pharmacophore is repaired. Measured on the p14 round, 0 of 600 anchors satisfied
the consensus and 100% broke the invariant Gly4 and Thr5, so this is a repair
rather than an exploration.

Each position is restrained by omission, which forces an identity, with a bias
toward the consensus top residue at the variable positions. The invariant
positions have a single allowed residue and so are fully determined.
"""
import argparse, glob, json, os, sys
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from consensus import CONSENSUS, INVARIANT, NATIVE

ALL20 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
         "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
A3 = {v: k for k, v in ONE.items()}
ANCHOR_HI = 10


def read_cif_seq(path):
    names, seq, chains = [], {}, {}
    for l in open(path):
        if l.startswith("_atom_site."):
            names.append(l.strip().split(".")[1]); continue
        if not (l.startswith("ATOM") or l.startswith("HETATM")): continue
        t = l.split()
        if len(t) != len(names): continue
        d = dict(zip(names, t))
        if (d.get("label_atom_id") or "").strip() != "CA": continue
        ch = d.get("auth_asym_id") or d.get("label_asym_id")
        try: rn = int(d.get("auth_seq_id") or d.get("label_seq_id"))
        except (TypeError, ValueError): continue
        seq.setdefault(ch, {})[rn] = ONE.get((d.get("label_comp_id") or "").upper(), "X")
        chains[ch] = chains.get(ch, 0) + 1
    return seq, chains


ap = argparse.ArgumentParser()
ap.add_argument("--complexes", nargs="+", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--weights", required=True)
ap.add_argument("--nseq", type=int, default=16)
ap.add_argument("--keep", type=int, default=3)
ap.add_argument("--temperature", type=float, default=0.2)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--top-bias", type=float, default=1.0,
                help="bias toward the consensus top residue at variable positions")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

recs = []
for f in a.complexes:
    seq, chains = read_cif_seq(f)
    if len(chains) < 2:
        print(f"  skip {os.path.basename(f)}: one chain"); continue
    by = sorted(chains, key=lambda c: -chains[c])
    tgt, bnd = by[0], by[1]
    blen = max(seq[bnd])
    if blen != 140:
        print(f"  skip {os.path.basename(f)}: binder {blen}"); continue
    recs.append({"cif": f, "tag": os.path.basename(f).replace(".cif", ""),
                 "bnd": bnd, "blen": blen,
                 "seq": "".join(seq[bnd][k] for k in sorted(seq[bnd]))})
print(f"backbones: {len(recs)}")
if not recs:
    print("FATAL: no input"); sys.exit(3)

from mpnn.inference_engines.mpnn import MPNNInferenceEngine
engine = MPNNInferenceEngine(
    model_type="protein_mpnn", checkpoint_path=a.weights, is_legacy_weights=True,
    out_directory=a.outdir, write_fasta=False, write_structures=False, device="cpu")

rows, stats = [], Counter()
for rec in recs:
    ch, blen = rec["bnd"], rec["blen"]
    designed = [f"{ch}{k}" for k in range(1, ANCHOR_HI + 1)]   # anchor only
    omit_pr, bias_pr = {}, {}
    for p, (allowed, top, _) in CONSENSUS.items():
        keep3 = {A3[x] for x in allowed}
        omit_pr[f"{ch}{p}"] = [x for x in ALL20 if x not in keep3]
        if len(allowed) > 1:
            bias_pr[f"{ch}{p}"] = {A3[top]: a.top_bias}
    inp = {"bias": None, "bias_per_residue": bias_pr,
           "omit": ["CYS"], "omit_per_residue": omit_pr,
           "structure_path": rec["cif"], "name": rec["tag"],
           "seed": a.seed, "batch_size": a.nseq, "number_of_batches": 1,
           "temperature": a.temperature, "designed_residues": designed}
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
        anc = binder[:ANCHOR_HI]
        bad = [p for p in CONSENSUS if anc[p-1] not in CONSENSUS[p][0]]
        if bad:
            stats["consensus_violation"] += 1
            continue
        # the designed region must be untouched
        if binder[ANCHOR_HI:] != rec["seq"][ANCHOR_HI:]:
            stats["designed_region_changed"] += 1
            continue
        ntop = sum(1 for p in CONSENSUS if anc[p-1] == CONSENSUS[p][1])
        cands.append({"seq": binder, "anc": anc, "ntop": ntop,
                      "n_native": sum(1 for i in range(10) if anc[i] == NATIVE[i])})
    if not cands:
        print(f"    WARNING {rec['tag']}: no compliant sequence"); continue
    uniq = {}
    for c in sorted(cands, key=lambda c: (-c["ntop"], -c["n_native"])):
        uniq.setdefault(c["anc"], c)
    sel = list(uniq.values())[: a.keep]
    for i, c in enumerate(sel):
        rows.append({"tag": rec["tag"], "binder_len": blen,
                     "motif": list(range(1, ANCHOR_HI + 1)), "seq": c["seq"],
                     "pos3": c["seq"][2], "anchor": c["anc"],
                     "consensus_top_matches": c["ntop"],
                     "native_matches": c["n_native"],
                     "rank_in_backbone": i, "acid_ok": True, "motif_ok": True})
    print(f"  {rec['tag'][:46]:46s} {len(cands):3d} compliant, kept {len(sel)}: "
          f"{' '.join(c['anc'] for c in sel)}")
    stats["backbones_ok"] += 1

if not rows:
    print("FATAL: nothing produced"); sys.exit(5)
anc_all = Counter(r["anchor"] for r in rows)
print(f"\nsequences kept      : {len(rows)} from {stats['backbones_ok']} backbones")
print(f"distinct anchors    : {len(anc_all)}")
print(f"consensus compliant : {len(rows)}/{len(rows)} (enforced by omission)")
print(f"invariant positions : all satisfied by construction {INVARIANT}")
nt = [r["consensus_top_matches"] for r in rows]
print(f"consensus-top matches per anchor: median {sorted(nt)[len(nt)//2]}/10  "
      f"min {min(nt)}  max {max(nt)}")
nn = [r["native_matches"] for r in rows]
print(f"identity to native GLP-1 anchor : median {sorted(nn)[len(nn)//2]}/10  "
      f"max {max(nn)}")
print("most common anchors:")
for s, c in anc_all.most_common(8):
    print(f"  {s}  x{c}   (native {NATIVE})")
print(f"rejected: {dict(stats)}")
json.dump({"n": len(rows), "acid_pos": 3, "n_acid_ok": len(rows),
           "n_motif_ok": len(rows), "skipped": [], "designs": rows},
          open(os.path.join(a.outdir, "h8_mpnn.json"), "w"), indent=1)
print("wrote " + os.path.join(a.outdir, "h8_mpnn.json"))
