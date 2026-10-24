#!/usr/bin/env python
"""Stage 3a -- build Boltz YAMLs for the designed binders against the receptor.

Receptor MSA: pass --receptor-msa with the cached AlphaFold3 a3m. Boltz 2.2.1 reads that
file directly (verified, job 1088573: 2/2 manifest records, 0 failures) and it is not
optional in practice -- the same binder scored ipTM 0.841 with the MSA and 0.316 without.
The designed binder always stays `msa: empty`: it is de novo and has no homologues.

The binder is written FIRST (chain A) and the receptor second (chain B); the analysis
step still verifies the assignment by length and sequence identity rather than trusting
the order, because Boltz renumbers and relabels chains.

"""
import argparse, json, os, sys

ap = argparse.ArgumentParser()
ap.add_argument("--mpnn-json", required=True)
ap.add_argument("--af3-request", required=True)
ap.add_argument("--out-dir", required=True)
ap.add_argument("--per-backbone", type=int, default=3)
ap.add_argument("--manifest", required=True, help="written ONE LEVEL ABOVE out-dir")
ap.add_argument("--receptor-msa", default="empty",
                help="path to the cached receptor a3m, or 'empty'. Measured 2026-10-18: "
                     "the AF3 a3m works directly and lifts ipTM 0.316 -> 0.841 on the same "
                     "binder, so 'empty' is for controls only.")
a = ap.parse_args()

req = json.load(open(a.af3_request))
seqs = req[0]["sequences"]
receptor = seqs[0]["proteinChain"]["sequence"]
peptide = seqs[1]["proteinChain"]["sequence"]
print(f"receptor {len(receptor)} aa, reference peptide {len(peptide)} aa")

d = json.load(open(a.mpnn_json))
by = {}
for r in d["designs"]:
    by.setdefault(r["tag"], []).append(r)

os.makedirs(a.out_dir, exist_ok=True)
man = []
for tag, rows in sorted(by.items()):
    for i, r in enumerate(rows[: a.per_backbone]):
        name = f"{tag}__s{i}"
        p = os.path.join(a.out_dir, name + ".yaml")
        with open(p, "w") as fh:
            fh.write("version: 1\nsequences:\n")
            fh.write("  - protein:\n      id: A\n")
            fh.write(f"      sequence: {r['seq']}\n      msa: empty\n")
            fh.write("  - protein:\n      id: B\n")
            fh.write(f"      sequence: {receptor}\n      msa: {a.receptor_msa}\n")
        man.append({"name": name, "tag": tag, "binder_seq": r["seq"],
                    "binder_len": r["binder_len"], "motif": r["motif"],
                    "acid_pos3": r.get("pos3")})

# the manifest must NOT sit in out_dir: `boltz predict` parses EVERY file in the
# directory it is given and dies on a stray .json
os.makedirs(os.path.dirname(a.manifest), exist_ok=True)
json.dump({"n": len(man), "receptor_len": len(receptor),
           "reference_peptide": peptide, "yamls": man},
          open(a.manifest, "w"), indent=1)
print(f"wrote {len(man)} yamls to {a.out_dir}")
print(f"manifest {a.manifest}")
