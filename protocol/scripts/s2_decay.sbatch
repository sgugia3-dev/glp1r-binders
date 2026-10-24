#!/usr/bin/env bash
#SBATCH --job-name=s2-decay
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/s2_decay_%j.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/s2_decay_%j.out
#
# STRATEGY s2: vary potentials.guide_decay at fixed binder_ROG,weight:1, guide_scale=10,
# with the 4 hotspots. Everything else identical to p8_rfd2_pot.sbatch arm 1.
# manager.get_guide_scale: t counts DOWN from T, scale = (t/T)^n * guide_scale, so
#   constant = full pressure to the final step, linear/quadratic/cubic fade progressively
#   earlier. Hypothesis: a faster fade sets global shape early and leaves the late steps
#   free for the hotspot conditioning to settle the interface.
# 4 configs run SEQUENTIALLY in this ONE job (no array: 20-job cap, 8 agents).
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/anaconda1/etc/profile.d/conda.sh
conda activate SE3nv

G=/lustre/home/gsgueglia/challenge2026/glp1r2026
R=/lustre/home/gsgueglia/diff_test/RFdiffusion
PDB="$G/inputs/rfd2/glp1r_ecd_A1_117.pdb"
NDES=${NDES:-6}
BASE=/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/s2_decay

echo "job ${SLURM_JOB_ID} host $(hostname) start $(date -Is)"
nvidia-smi --query-gpu=driver_version --format=csv,noheader | head -1
python -c "import torch;print('torch',torch.__version__,'cuda avail',torch.cuda.is_available())" || exit 1
cd "$R" || exit 2

for DECAY in constant linear quadratic cubic; do
  OUT="$BASE/$DECAY"
  for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
  echo "=== CONFIG decay=$DECAY START $(date -Is) ==="
  time python scripts/run_inference.py \
    "inference.input_pdb=$PDB" \
    "inference.output_prefix=$OUT/des" \
    "inference.num_designs=$NDES" \
    'contigmap.contigs=[A1-117/0 70-70]' \
    'ppi.hotspot_res=[A24,A25,A27,A42]' \
    'potentials.guiding_potentials=["type:binder_ROG,weight:1"]' \
    potentials.guide_scale=10 "potentials.guide_decay=$DECAY"
  rc=$?
  echo "=== CONFIG decay=$DECAY rc=$rc designs=$(ls "$OUT"/des_*.pdb 2>/dev/null | wc -l) END $(date -Is) ==="
done
echo "=== ALL DONE $(date -Is) ==="
