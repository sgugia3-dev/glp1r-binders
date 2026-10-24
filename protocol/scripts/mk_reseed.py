#!/usr/bin/env python
"""Build Boltz inputs to re-fold the top p15 designs at a new seed.

Inputs are byte-identical to the originals -- same binder sequence, same receptor,
same cached MSA -- so the only difference between runs is the Boltz seed. That is
what makes the comparison a seed-replicate test rather than a re-design.

The manifest keeps the ORIGINAL design name so every downstream table joins
straight back to the seed-1 results.
"""
import csv, glob, json, os, re, sys

A, CPLX, TOPN, SHARD, NSH, OUTDIR, MANIFEST = (
    sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]),
    int(sys.argv[5]), sys.argv[6], sys.argv[7])

ranked = list(csv.DictReader(open(os.path.join(A, "p15_RANKED.csv"))))
ranked.sort(key=lambda r: -float(r["score"]))
top = ranked[:TOPN]
print(f"{len(ranked)} ranked designs -> top {len(top)}")

# binder sequences come from the original run manifests
seqs = {}
for m in sorted(glob.glob(os.path.join(CPLX, "sh*", "manifest.json"))):
    sh = os.path.basename(os.path.dirname(m)).replace("sh", "")
    d = json.load(open(m))
    for y in d["yamls"]:
        seqs[(sh, y["name"])] = y
print(f"manifest entries: {len(seqs)}")

req = json.load(open(os.environ["AF3_REQUEST"]))
receptor = req[0]["sequences"][0]["proteinChain"]["sequence"]
msa = os.environ["RECEPTOR_MSA"]

mine = [r for i, r in enumerate(top) if i % NSH == SHARD]
os.makedirs(OUTDIR, exist_ok=True)
man, missing = [], 0
for r in mine:
    k = (r["shard"], r["name"])
    y = seqs.get(k)
    if not y:
        missing += 1; continue
    name = r["name"]
    with open(os.path.join(OUTDIR, name + ".yaml"), "w") as fh:
        fh.write("version: 1\nsequences:\n")
        fh.write("  - protein:\n      id: A\n")
        fh.write(f"      sequence: {y['binder_seq']}\n      msa: empty\n")
        fh.write("  - protein:\n      id: B\n")
        fh.write(f"      sequence: {receptor}\n      msa: {msa}\n")
    man.append({"name": name, "tag": y.get("tag", name),
                "binder_seq": y["binder_seq"], "binder_len": y["binder_len"],
                "motif": y["motif"], "acid_pos3": y.get("acid_pos3"),
                "orig_shard": r["shard"], "orig_score": r["score"]})
os.makedirs(os.path.dirname(MANIFEST), exist_ok=True)
json.dump({"n": len(man), "receptor_len": len(receptor),
           "reference_peptide": req[0]["sequences"][1]["proteinChain"]["sequence"],
           "yamls": man}, open(MANIFEST, "w"), indent=1)
print(f"shard {SHARD}/{NSH}: wrote {len(man)} yamls ({missing} without a manifest entry)")
