#!/usr/bin/env python
"""Write the ten queried designs out as a self-describing deliverable."""
import csv, json, os, shutil, sys

SRC, OUT = sys.argv[1], sys.argv[2]
d = json.load(open(os.path.join(SRC, "data.json")))
des = d["designs"]
REP = {"s1": "seed1", "s2": "seed2", "s3": "seed3"}
for sub in REP.values():
    os.makedirs(os.path.join(OUT, "structures", sub), exist_ok=True)

bundles = {k: json.load(open(os.path.join(SRC, f"q_{k}.json"))) for k in REP}
n = 0
for x in des:
    base = f"{x['rank']:02d}_{x['parent']}_{x['source']}"
    for k, sub in REP.items():
        txt = bundles[k].get(x["tag"])
        if not txt: continue
        p = os.path.join(OUT, "structures", sub, base + f"_{sub}.pdb")
        open(p, "w").write(txt); n += 1

with open(os.path.join(OUT, "sequences.fasta"), "w") as fh:
    for x in des:
        fh.write(f">{x['rank']:02d}_{x['parent']}_{x['source']} "
                 f"anchor={x['anchor']} mean_ipsae={x.get('mean_ipsae','')} "
                 f"gate={x.get('n_seeds_gate','')}/3 len={len(x['seq'])}\n")
        for i in range(0, len(x["seq"]), 60):
            fh.write(x["seq"][i:i+60] + "\n")

cols = ["rank", "parent", "source", "name", "anchor", "n_seeds_gate", "n_seeds_full",
        "n_hotspots_total", "mean_ipsae", "range_ipsae", "mean_plddt",
        "mean_motif_rmsd", "range_motif_rmsd"]
per = ["motif", "recept", "hits", "sb", "hc", "ipsae", "plddt", "iptm"]
for k in REP:
    cols += [f"{REP[k]}_{m}" for m in per]
cols += ["sequence"]
with open(os.path.join(OUT, "metrics.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(cols)
    for x in des:
        row = [x["rank"], x["parent"], x["source"], x["name"], x["anchor"],
               x.get("n_seeds_gate"), x.get("n_seeds_full"), x.get("n_total"),
               x.get("mean_ipsae"), x.get("range_ipsae"), x.get("mean_plddt"),
               x.get("mean_motif"), x.get("range_motif")]
        for k in REP:
            s = x.get(k, {})
            row += [s.get(m, "") for m in per]
        row.append(x["seq"])
        w.writerow(row)

m = d["meta"]
open(os.path.join(OUT, "README.txt"), "w").write(f"""Ten GLP-1R binder candidates
============================
The ten sequences supplied, matched by exact sequence against every design this
project produced, with all three independent predictions of each.

CONTENTS
  sequences.fasta        the ten binder sequences, 140 aa each
  metrics.csv            every metric, per replicate, one row per candidate
  structures/seed1/      prediction at replicate 1 (run {m['replicates'][0]})
  structures/seed2/      prediction at replicate 2 (run {m['replicates'][1]})
  structures/seed3/      prediction at replicate 3 (run {m['replicates'][2]})

STRUCTURE FILES
  PDB. Chain A = binder, residues 1-140 (1-10 is the anchor).
  Chain B = GLP1R, residues 1-409. Both chains start at 1.
  The B-factor column carries per-residue pLDDT on a 0-100 scale.

WHAT THE REPLICATES ARE
  Three Boltz predictions from byte-identical input, differing only in seed.
  Disagreement between them is prediction variance, not a difference in design.
  {m['n_gate3']} of {m['n_designs']} hold the acceptance gate in all three.
  {m['n_full3']} of {m['n_designs']} keep full hotspot contact in all three.

  Gate: anchor RMSD < 2.0 A AND >= 3 hotspots contacted AND receptor RMSD < 5.0 A.

COHORT FIGURES (median across the ten)
  mean ipSAE          {m['median_mean_ipsae']}   (spread across replicates {m['median_range_ipsae']})
  mean binder pLDDT   {m['median_mean_plddt']}
  mean anchor RMSD    {m['median_mean_motif']} A  (spread {m['median_range_motif']} A)

CAVEAT
  Hotspot set differs by lineage: p3 uses A86-A89, p10 and p12 use A27/A42/A86,
  p13 uses A27/A42/A86/A95. The `n_hotspots_total` column records which applies,
  so a design showing 3 of 3 is not weaker than one showing 4 of 4.
""")
print(f"wrote {n} structure files for {len(des)} candidates")
for f in ("sequences.fasta", "metrics.csv", "README.txt"):
    print(f"  {f}: {os.path.getsize(os.path.join(OUT,f))} bytes")
