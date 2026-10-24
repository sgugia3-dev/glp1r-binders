#!/bin/bash
#SBATCH -p htc
#SBATCH --qos=htc
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -J bestdirs
set -e
ulimit -u 8192
export OMP_NUM_THREADS=1
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate rfd3

S=$C26_ROOT/glp1r2026/scripts
A=$C26_ROOT/glp1r2026/analysis
CONTIG="B1-10,130-130,/0,A1-409"

echo "############ P3 Trp87 86-89 patch"
python "$S/best_dir.py" --ranked "$A/p3_RANKED.csv" \
  --root "$C26_SCRATCH/runs/glp1r/p3_1088996" \
  --out "$C26_ROOT/glp1r2026/p3_best" \
  --prefix p3 --campaign p3_1088996 --analysis "$A" --order tier \
  --contig "$CONTIG" \
  --hotspots "A86:CD2, A87:CD1, A88:CD1, A89:NE2" \
  --ipsae-note "median ipSAE 0.894 at 4/4 hotspots vs 0.870 at 0 hotspots, a 0.024 separation over 480 folds"

echo
echo "############ P2 24/25/27/42 patch, regenerated with the same tool"
python "$S/best_dir.py" --ranked "$A/p2_RANKED.csv" \
  --root "$C26_SCRATCH/runs/glp1r/p2_1088883" \
  --out "$C26_ROOT/glp1r2026/p2_best" \
  --prefix p2 --campaign p2_1088883 --analysis "$A" --order tier \
  --contig "$CONTIG" \
  --hotspots "A24:OE1,NE2  A25:NH1,NH2  A27:CD1  A42:OG1" \
  --ipsae-note "median ipSAE 0.882 at 4/4 hotspots vs 0.871 at 0 hotspots, a 0.011 separation over 1200 folds"
echo
echo "=== done"
