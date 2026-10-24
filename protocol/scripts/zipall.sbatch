#!/bin/bash
#SBATCH -p htc
#SBATCH --qos=htc
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH -J zipall
set -e
ulimit -u 8192
export OMP_NUM_THREADS=1
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate rfd3
G=$C26_ROOT/glp1r2026
rm -rf $G/p14_structures $G/p14_structures.zip $G/p14_band_structures.zip
python $G/scripts/collect_all.py $G/analysis $C26_SCRATCH/runs/glp1r/p14_1092297 $G/p14_structures
cd $G
echo "--- raw sizes"
du -sh p14_structures p14_structures/* | sed 's/^/  /'
echo "--- zipping (band + gated_full first, it is the smaller useful one)"
zip -qr9 p14_band_structures.zip p14_structures/band_q4to9 p14_structures/gated_full \
    p14_structures/INDEX.csv p14_structures/README.txt
zip -qr9 p14_structures.zip p14_structures
ls -lh p14_band_structures.zip p14_structures.zip | awk '{print "  ", $5, $9}'
