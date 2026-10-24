#!/usr/bin/env bash
#SBATCH --job-name=s7-length
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/s7_length_%j.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/s7_length_%j.out
#
# STRATEGY: binder LENGTH sweep. Everything fixed at the reference config
# (binder_ROG weight:1, guide_scale=10, quadratic decay, 4 hotspots, 6 designs);
# the ONLY knob is the designed binder length L in the contig: [A1-117/0 L-L].
# Four configs run SEQUENTIALLY in this single job (20-job cap, 8 agents).
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/anaconda1/etc/profile.d/conda.sh
conda activate SE3nv

G=/lustre/home/gsgueglia/challenge2026/glp1r2026
R=/lustre/home/gsgueglia/diff_test/RFdiffusion
PDB="$G/inputs/rfd2/glp1r_ecd_A1_117.pdb"
BASE=/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/s7_length
NDES=${NDES:-6}

echo "node $(hostname)  job $SLURM_JOB_ID  ndes=$NDES"
nvidia-smi --query-gpu=driver_version,name --format=csv,noheader | head -1
python -c "import torch;print('torch',torch.__version__,'cuda',torch.cuda.is_available())" || exit 1

cd "$R" || exit 2
for L in 50 70 90 110; do
  OUT="$BASE/L$L"
  for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
  echo "=== CONFIG L=$L START $(date +%T) ==="
  time python scripts/run_inference.py \
    inference.input_pdb="$PDB" \
    inference.output_prefix="$OUT/des" \
    inference.num_designs=$NDES \
    "contigmap.contigs=[A1-117/0 ${L}-${L}]" \
    'ppi.hotspot_res=[A24,A25,A27,A42]' \
    'potentials.guiding_potentials=["type:binder_ROG,weight:1"]' \
    potentials.guide_scale=10 potentials.guide_decay=quadratic
  rc=$?
  echo "=== CONFIG L=$L DONE rc=$rc designs=$(ls "$OUT"/des_*.pdb 2>/dev/null | wc -l) $(date +%T) ==="
done
echo "=== ALL CONFIGS DONE ==="
