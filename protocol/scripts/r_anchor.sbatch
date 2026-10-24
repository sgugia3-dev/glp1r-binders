#!/usr/bin/env bash
#SBATCH --job-name=r-anchor
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/r_anchor_%A_%a.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/r_anchor_%A_%a.out
#
# ANCHORED arms only: the peptide motif B1-10 is kept in every one, N-terminus free,
# extension C-terminal. The variables are extension LENGTH and WHICH PATCH is targeted.
# Contacts are scored from the DESIGNED region only, so a contact made by the fixed motif
# does not count.
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate rfd3
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
G="$C26_ROOT/glp1r2026"
ARMS=(r1_ecd_L40_90 r2_ecd_L90 r3_ecd_L130 r4_prox_L40_90 r5_prox_L90 r6_nohot_L40_90)
HOTS=(24,25,27,42 24,25,27,42 24,25,27,42 174,178,210,214,218 174,178,210,214,218 24,25,27,42)
S=${SLURM_ARRAY_TASK_ID:-0}
A=${ARMS[$S]}; HS=${HOTS[$S]}
OUT="$C26_SCRATCH/runs/glp1r/r_anchor_${SLURM_ARRAY_JOB_ID:-m}/$A"
for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
CKPT=$(ls "$C26_WEIGHTS/rfd3"/*.ckpt | head -1)
echo "arm $A  hotspots=$HS  node $(hostname)  out $OUT"
time rfd3 design out_dir="$OUT" inputs="$G/inputs/rfd3/$A.json" ckpt_path="$CKPT" \
  n_batches=1 diffusion_batch_size="${BS:-2}" seed=$((5150 + S)) \
  inference_sampler.step_scale=3 inference_sampler.gamma_0=0.2 \
  prevalidate_inputs=True skip_existing=True
rc=$?
echo "designs: $(find "$OUT" -name '*.cif.gz' | wc -l)"
echo -n "RESULT $A  "; python "$G/scripts/r_anchor_measure.py" "$OUT" "$HS"
echo "=== exit: $rc ==="
exit $rc
