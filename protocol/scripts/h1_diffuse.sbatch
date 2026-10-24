#!/usr/bin/env bash
#SBATCH --job-name=h1-glp-diff
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/h1_diff_%j.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/h1_diff_%j.out
#
# Smoke test, stage 1: motif-scaffold the peptide's first 10 residues into a 50-100 aa
# miniprotein, 5 input models x BATCH designs each.
# --cpus-per-task=2 deliberately: these nodes do not oversubscribe, and a large core
# request strands the other GPUs on the node.
set -o pipefail
ulimit -u 8192 2>/dev/null || true          # per-user nproc is 100 on compute nodes too
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate rfd3
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

G="$C26_ROOT/glp1r2026"
BATCH="${BATCH:-2}"
OUT="$C26_SCRATCH/runs/glp1r/h1_diff_${SLURM_JOB_ID}"
for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done

echo "node $(hostname)  out $OUT"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader | head -1
CKPT=$(ls "$C26_WEIGHTS/rfd3"/*.ckpt 2>/dev/null | head -1)
if [ -z "$CKPT" ]; then
  foundry install rfd3 --checkpoint-dir "$C26_WEIGHTS/rfd3" || true
  CKPT=$(ls "$C26_WEIGHTS/rfd3"/*.ckpt 2>/dev/null | head -1)
fi
[ -n "$CKPT" ] || { echo "FATAL: no rfd3 checkpoint"; echo "=== exit: 1 ==="; exit 1; }
echo "ckpt $CKPT"

# peak-VRAM sampler: the target is 409 aa, far larger than anything measured here,
# so OOM is the most likely failure mode and we want the number either way
( while :; do nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits >>"$OUT/.vram"; sleep 2; done ) &
VPID=$!; trap 'kill $VPID 2>/dev/null' EXIT

time rfd3 design \
  out_dir="$OUT" \
  inputs="${INPUTS:-$G/inputs/rfd3/glp1r_cterm.json}" \
  ckpt_path="$CKPT" \
  n_batches=1 \
  diffusion_batch_size="$BATCH" \
  inference_sampler.step_scale=3 \
  inference_sampler.gamma_0=0.2 \
  prevalidate_inputs=True \
  skip_existing=False
rc=$?

echo "designs produced: $(find "$OUT" -name '*.cif.gz' 2>/dev/null | wc -l) (expected $((BATCH*5)))"
echo "PEAK VRAM: $(sort -n "$OUT/.vram" 2>/dev/null | tail -1) MiB of 32768"
echo "$OUT" > "$G/logs/h1_diff_outdir.txt"
echo "=== exit: $rc ==="
exit $rc
