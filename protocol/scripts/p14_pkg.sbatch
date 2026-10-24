#!/bin/bash
#SBATCH -p htc
#SBATCH --qos=htc
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -J p14pkg
set -e
ulimit -u 8192
export OMP_NUM_THREADS=1
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate rfd3
G=$C26_ROOT/glp1r2026
rm -rf $G/p14_best
python $G/scripts/best_dir.py \
  --ranked $G/analysis/p14_RANKED.csv \
  --root $C26_SCRATCH/runs/glp1r/p14_1092297 \
  --out $G/p14_best \
  --prefix p14 --campaign p14_1092297 --analysis $G/analysis --order tier --hotspot-count 4 \
  --contig "redesigned from p3/p10/p12/p13_best backbones, anchor designed with His1 and Asp-Glu3" \
  --hotspots "p3 A86-89 / p10 and p12 A27-A42-A86 / p13 A27-A42-A86-A95" \
  --ipsae-note "charge-complementary round: the anchor was redesigned and the surface charge inverted"
