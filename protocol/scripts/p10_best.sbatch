#!/bin/bash
#SBATCH -p htc
#SBATCH --qos=htc
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -J p10best
set -e
ulimit -u 8192
export OMP_NUM_THREADS=1
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate rfd3
S=$C26_ROOT/glp1r2026/scripts
A=$C26_ROOT/glp1r2026/analysis
ROOT=$C26_SCRATCH/runs/glp1r/p10_1089071
TOOL=$C26_ROOT/tools/ipsae/ipsae.py

echo "=== stage 1: ipSAE over every P10 fold"
n=0; ok=0; skip=0
for d in $ROOT/sh*/boltz_results_yaml/predictions/*/; do
  nm=$(basename $d)
  cif=$d/${nm}_model_0.cif; pae=$d/pae_${nm}_model_0.npz
  [ -f "$cif" ] && [ -f "$pae" ] || continue
  n=$((n+1))
  if ls $d/*_15_15.txt >/dev/null 2>&1; then skip=$((skip+1)); continue; fi
  python $TOOL "$pae" "$cif" 15 15 >/dev/null && ok=$((ok+1))
done
echo "folds $n  newly scored $ok  already had $skip"

echo; echo "=== stage 2: rank"
python $S/rank_campaign.py $A $ROOT $A/p10_RANKED.csv p10_1089071

echo; echo "=== stage 3: best dir"
rm -rf $C26_ROOT/glp1r2026/p10_best
python $S/best_dir.py --ranked $A/p10_RANKED.csv --root $ROOT \
  --out $C26_ROOT/glp1r2026/p10_best \
  --prefix p10 --campaign p10_1089071 --analysis $A --order tier \
  --contig "B1-10,130-130,/0,A1-409" \
  --hotspots "A27:CD1, A42:OG1, A86:CD2" \
  --ipsae-note "see analysis/p10_RANKED.csv; the per-hotspot geometric check is what verifies engagement"
echo; echo "=== done"
