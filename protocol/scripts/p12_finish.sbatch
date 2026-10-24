#!/bin/bash
#SBATCH -p htc
#SBATCH --qos=htc
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 03:00:00
#SBATCH -J p12fin
set -e
ulimit -u 8192
export OMP_NUM_THREADS=1
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate rfd3
S=$C26_ROOT/glp1r2026/scripts
A=$C26_ROOT/glp1r2026/analysis
ROOT=$C26_SCRATCH/runs/glp1r/p12_1089085

echo "########## stage 1: aggregate + bridging test"
python $S/campaign_agg.py $A p12_1089085

echo; echo "########## stage 2: ipSAE over every fold"
n=0; ok=0; skip=0
for d in $ROOT/sh*/boltz_results_yaml/predictions/*/; do
  nm=$(basename $d)
  cif=$d/${nm}_model_0.cif; pae=$d/pae_${nm}_model_0.npz
  [ -f "$cif" ] && [ -f "$pae" ] || continue
  n=$((n+1))
  if ls $d/*_15_15.txt >/dev/null 2>&1; then skip=$((skip+1)); continue; fi
  python $C26_ROOT/tools/ipsae/ipsae.py "$pae" "$cif" 15 15 >/dev/null && ok=$((ok+1))
  [ $((n % 300)) -eq 0 ] && echo "  $n scored, ${SECONDS}s"
done
echo "folds $n  newly scored $ok  already had $skip"

echo; echo "########## stage 3: rank"
python $S/rank_campaign.py $A $ROOT $A/p12_RANKED.csv p12_1089085

echo; echo "########## stage 4: best dir"
rm -rf $C26_ROOT/glp1r2026/p12_best
python $S/best_dir.py --ranked $A/p12_RANKED.csv --root $ROOT \
  --out $C26_ROOT/glp1r2026/p12_best \
  --prefix p12 --campaign p12_1089085 --analysis $A --order tier \
  --contig "B1-10,130-130,/0,A1-409" \
  --hotspots "A27:CD1, A42:OG1, A86:CD2" \
  --ipsae-note "measured per hotspot-count bin in the rank log; the geometric per-hotspot check is what verifies engagement"
echo; echo "########## done"
