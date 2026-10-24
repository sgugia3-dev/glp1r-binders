#!/bin/bash
#SBATCH -p htc
#SBATCH --qos=htc
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -J p14fin
set -e
ulimit -u 8192
export OMP_NUM_THREADS=1
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate rfd3
S=$C26_ROOT/glp1r2026/scripts
A=$C26_ROOT/glp1r2026/analysis
PC=$A/physchem
ROOT=$C26_SCRATCH/runs/glp1r/p14_1092297
mkdir -p $PC

echo "########## stage 1: pre-redesign interaction baseline on the 200 inputs"
for D in p3_best p10_best p12_best p13_best; do
  python $S/h7_interactions.py --flat \
    --pred-dir $C26_ROOT/glp1r2026/$D/structures \
    --out-csv $PC/baseline_${D}_interactions.csv | tail -8
done

echo; echo "########## stage 2: aggregate the redesign round"
python $S/p14_agg.py $A p14_1092297

echo; echo "########## stage 3: ipSAE over every fold"
n=0; ok=0
for d in $ROOT/sh*/boltz_results_yaml/predictions/*/; do
  nm=$(basename $d)
  cif=$d/${nm}_model_0.cif; pae=$d/pae_${nm}_model_0.npz
  [ -f "$cif" ] && [ -f "$pae" ] || continue
  n=$((n+1))
  ls $d/*_15_15.txt >/dev/null 2>&1 && continue
  python $C26_ROOT/tools/ipsae/ipsae.py "$pae" "$cif" 15 15 >/dev/null && ok=$((ok+1))
  [ $((n % 200)) -eq 0 ] && echo "  $n scored, ${SECONDS}s"
done
echo "folds $n  newly scored $ok"

echo; echo "########## stage 4: rank"
python $S/rank_campaign.py $A $ROOT $A/p14_RANKED.csv p14_1092297

echo; echo "########## stage 5: best dir"
rm -rf $C26_ROOT/glp1r2026/p14_best
python $S/best_dir.py --ranked $A/p14_RANKED.csv --root $ROOT \
  --out $C26_ROOT/glp1r2026/p14_best \
  --prefix p14 --campaign p14_1092297 --analysis $A --order tier \
  --contig "redesigned from p3/p10/p12/p13_best backbones" \
  --hotspots "per source: p3 A86-89, p10/p12 A27/42/86, p13 A27/42/86/95" \
  --ipsae-note "this round redesigned the anchor and modulated charge; see p14_RANKED.csv"
echo; echo "########## done"
