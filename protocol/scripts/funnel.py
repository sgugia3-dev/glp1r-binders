#!/usr/bin/env python
"""Survival funnel for every GLP-1R campaign.

Counts BACKBONES as well as folds. ProteinMPNN produced several sequences per
diffused backbone (the __sN suffix), so each backbone was refolded 2-3 times and
a fold count overstates the number of distinct candidates.

Stages, in order of increasing stringency:
  refolded        every Boltz prediction that was measured
  motif < 2 A     anchor placed correctly against the AF3 peptide
  full epitope    every hotspot in that campaign's set contacted (3/3 or 4/4)
  GATE            motif < 2 A AND >= 3 hotspots AND receptor RMSD < 5 A
  ELITE           GATE plus full epitope AND motif < 1 A
  packaged        distinct backbones in the delivered *_best directory
"""
import csv, glob, os, re, sys
from collections import defaultdict

A = sys.argv[1]
CAMPAIGNS = [("p2_1088883", "P2  A24/25/27/42", 4, "p2_best"),
             ("p3_1088996", "P3  A86-89", 4, "p3_best"),
             ("p10_1089071", "P10 A27/42/86 (smoke)", 3, "p10_best"),
             ("p12_1089085", "P12 A27/42/86", 3, "p12_best"),
             ("p13_1089119", "P13 A27/42/86/95", 4, "p13_best")]
G = "/lustre/home/gsgueglia/challenge2026/glp1r2026"


def load(run, pat):
    out = {}
    for f in sorted(glob.glob(os.path.join(A, f"{run}_sh*_{pat}.csv"))):
        sh = re.search(r"_sh(\d+)_", f).group(1)
        for r in csv.DictReader(open(f)):
            out[(sh, r["name"])] = r
    return out


def backbone(shard, name):
    """Strip the MPNN sequence index so sibling sequences collapse to one backbone."""
    return (shard, re.sub(r"__s\d+$", "", name))


print(f"{'campaign':24s} {'folds':>6s} {'bb':>5s} | {'motif<2A':>9s} {'bb':>5s} | "
      f"{'full':>6s} {'bb':>5s} | {'GATE':>6s} {'bb':>5s} | {'ELITE':>6s} {'bb':>5s} | "
      f"{'pkg':>4s}")
tot = defaultdict(int)
for run, lbl, nh, bd in CAMPAIGNS:
    dev, epi = load(run, "deviation"), load(run, "epitope")
    keys = sorted(set(dev) & set(epi))
    rows = []
    for k in keys:
        try:
            rows.append({"bb": backbone(*k),
                         "motif": float(dev[k]["motif_rmsd_vs_af3"]),
                         "recept": float(dev[k]["receptor_rmsd"]),
                         "hits": int(epi[k]["n_hotspots_hit"])})
        except (ValueError, KeyError):
            pass
    u2 = [r for r in rows if r["motif"] < 2.0]
    fu = [r for r in rows if r["hits"] == nh]
    ga = [r for r in rows if r["motif"] < 2.0 and r["hits"] >= 3 and r["recept"] < 5.0]
    el = [r for r in ga if r["hits"] == nh and r["motif"] < 1.0]
    nb = lambda s: len({r["bb"] for r in s})
    pkg = 0
    mf = os.path.join(G, bd, "metrics.csv")
    if os.path.exists(mf):
        pkg = len(list(csv.DictReader(open(mf))))
    print(f"{lbl:24s} {len(rows):6d} {nb(rows):5d} | {len(u2):9d} {nb(u2):5d} | "
          f"{len(fu):6d} {nb(fu):5d} | {len(ga):6d} {nb(ga):5d} | "
          f"{len(el):6d} {nb(el):5d} | {pkg:4d}")
    for k, v in (("folds", len(rows)), ("bb", nb(rows)), ("gate", len(ga)),
                 ("gate_bb", nb(ga)), ("elite", len(el)), ("elite_bb", nb(el)),
                 ("full", len(fu)), ("full_bb", nb(fu)), ("pkg", pkg)):
        tot[k] += v
print("-" * 118)
print(f"{'TOTAL':24s} {tot['folds']:6d} {tot['bb']:5d} | {'':9s} {'':5s} | "
      f"{tot['full']:6d} {tot['full_bb']:5d} | {tot['gate']:6d} {tot['gate_bb']:5d} | "
      f"{tot['elite']:6d} {tot['elite_bb']:5d} | {tot['pkg']:4d}")
print("\nbb = distinct diffused backbones; folds = Boltz predictions (2-3 sequences per backbone)")
print("'full' means every hotspot in that campaign's own set: 4/4 for P2/P3/P13, 3/3 for P10/P12")
