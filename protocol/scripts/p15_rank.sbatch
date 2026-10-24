#!/bin/bash
#SBATCH -p htc
#SBATCH --qos=htc
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -J p15rank
set -e
ulimit -u 8192
export OMP_NUM_THREADS=1
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate rfd3
G=$C26_ROOT/glp1r2026
A=$G/analysis
ROOT=$C26_SCRATCH/runs/glp1r/p15_1099477

echo "########## stage 1: verify the anchor repair + merge epitope tables"
python $G/scripts/p15_check.py $A p15_1099477

echo; echo "########## stage 2: ipSAE"
n=0; ok=0
for d in $ROOT/sh*/boltz_results_yaml/predictions/*/; do
  nm=$(basename $d)
  cif=$d/${nm}_model_0.cif; pae=$d/pae_${nm}_model_0.npz
  [ -f "$cif" ] && [ -f "$pae" ] || continue
  n=$((n+1))
  ls $d/*_15_15.txt >/dev/null 2>&1 && continue
  python $C26_ROOT/tools/ipsae/ipsae.py "$pae" "$cif" 15 15 >/dev/null && ok=$((ok+1))
done
echo "folds $n  newly scored $ok"

echo; echo "########## stage 3: rank"
python $G/scripts/rank_campaign.py $A $ROOT $A/p15_RANKED.csv p15_1099477

echo; echo "########## stage 4: package"
rm -rf $G/p15_best
python $G/scripts/best_dir.py --ranked $A/p15_RANKED.csv --root $ROOT \
  --out $G/p15_best --prefix p15 --campaign p15_1099477 --analysis $A \
  --order tier --hotspot-count 4 \
  --contig "anchor repaired under the 23-agonist consensus; residues 11-140 inherited from p14" \
  --hotspots "p3 A86-89 / p10 and p12 A27-A42-A86 / p13 A27-A42-A86-A95" \
  --ipsae-note "anchor restrained to the agonist consensus at all ten positions"
echo; echo "########## done"
