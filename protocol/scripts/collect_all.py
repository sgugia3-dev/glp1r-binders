#!/usr/bin/env python
"""Collect every predicted structure from the p14 round into one tree.

Named so the files are self-describing without the metrics table: source lineage,
band membership, gate status and score are all in the filename.
"""
import csv, glob, os, re, shutil, sys

A, ROOT, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
meta = {}
for r in csv.DictReader(open(os.path.join(A, "p14_band_q4to9.csv"))):
    meta[(r["shard"], r["name"])] = r
ranked = {}
for r in csv.DictReader(open(os.path.join(A, "p14_RANKED.csv"))):
    ranked[(r["shard"], r["name"])] = r
epi = {}
for f in sorted(glob.glob(os.path.join(A, "p14_1092297_sh*_epitope.csv"))):
    sh = re.search(r"_sh(\d+)_", f).group(1)
    for r in csv.DictReader(open(f)):
        epi[(sh, r["name"])] = r

for sub in ("all_complexes", "band_q4to9", "gated_full"):
    os.makedirs(os.path.join(OUT, sub), exist_ok=True)

n = nb = ng = 0
index = []
for d in sorted(glob.glob(os.path.join(ROOT, "sh*", "boltz_results_yaml",
                                       "predictions", "*"))):
    if not os.path.isdir(d): continue
    nm = os.path.basename(d)
    sh = re.search(r"/sh(\d+)/", d).group(1)
    cif = glob.glob(os.path.join(d, "*_model_0.cif"))
    if not cif: continue
    k = (sh, nm)
    e = epi.get(k, {})
    rk = ranked.get(k)
    src = e.get("source", "unk")
    hits = e.get("n_hotspots_hit", "?")
    tot = e.get("n_hotspots_total", "?")
    score = rk["score"] if rk else ""
    tag = f"{src}_sh{sh}_{nm}_hits{hits}of{tot}"
    if score: tag += f"_score{float(score):.3f}"
    shutil.copy(cif[0], os.path.join(OUT, "all_complexes", tag + ".cif"))
    n += 1
    inband = k in meta
    if inband:
        shutil.copy(cif[0], os.path.join(OUT, "band_q4to9", tag + ".cif")); nb += 1
    full = str(hits) == str(tot)
    if rk and full:
        shutil.copy(cif[0], os.path.join(OUT, "gated_full", tag + ".cif")); ng += 1
    index.append({"file": tag + ".cif", "shard": sh, "name": nm, "source": src,
                  "hits": hits, "n_hotspots": tot, "gated": bool(rk),
                  "in_band_q4to9": inband, "full_epitope": full,
                  "score": score,
                  "q_designed": meta.get(k, {}).get("q_des", ""),
                  "motif_rmsd": meta.get(k, {}).get("motif",
                                 rk["motif_rmsd"] if rk else ""),
                  "salt_bridges": meta.get(k, {}).get("sb", ""),
                  "ipsae": rk["ipsae"] if rk else "",
                  "plddt_binder": rk["plddt_binder"] if rk else ""})

cols = list(index[0]) if index else []
with open(os.path.join(OUT, "INDEX.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(index)

open(os.path.join(OUT, "README.txt"), "w").write(f"""GLP-1R p14 predicted structures (campaign p14_1092297)
======================================================
Boltz-2 complexes, mmCIF. Chain A = binder residues 1-140
(1-10 is the redesigned anchor, 11-140 the designed region),
chain B = GLP1R 1-409, renumbered 141-549 in the file.

  all_complexes/   {n} every refold in the round
  band_q4to9/      {nb} subset with designed-region charge q(11-140) +4..+9
  gated_full/      {ng} passed the gate AND contact every hotspot of their source
  INDEX.csv        one row per structure: source, hits, score, charge,
                   motif RMSD, salt bridges, ipSAE, binder pLDDT

Filenames carry source lineage, shard, design name, hotspot hits and score.
Source lineages and their hotspot sets:
  p3   A86-89          p10  A27/A42/A86
  p12  A27/A42/A86     p13  A27/A42/A86/A95
""")
print(f"all_complexes {n}   band_q4to9 {nb}   gated_full {ng}")
print(f"wrote {OUT}")
