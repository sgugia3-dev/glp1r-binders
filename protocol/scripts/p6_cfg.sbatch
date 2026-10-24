#!/usr/bin/env bash
#SBATCH --job-name=p6-cfg
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p6_cfg_%A_%a.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p6_cfg_%A_%a.out
#
# CFG scale sweep, to find out whether the one guidance knob rfd3 DOES expose moves
# compactness. rfd3 has classifier-free guidance (inference_sampler.cfg_scale, default 2.0;
# applied as delta += (cfg_scale-1)*(delta - delta_uncond), inference_sampler.py:306-307).
# It is NOT a geometric potential: there is no monomer_ROG equivalent anywhere in the
# package. CFG sharpens adherence to the CONDITIONING (motif + hotspots).
#
# Hypothesis being tested: our conditioning is what stretches the chain (motif in the TM
# core, hotspots 34-40 A away in the ECD), so LOWER cfg_scale should let the chain fold
# more freely and compact, at the cost of on-target engagement. Higher should do the
# opposite. One array task per cfg value, same seed, so the only variable is cfg_scale.
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate rfd3
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

G="$C26_ROOT/glp1r2026"
CFGS=(1.0 1.5 2.0 3.0 5.0)
SH=${SLURM_ARRAY_TASK_ID:-0}
CFG=${CFGS[$SH]}
OUT="$C26_SCRATCH/runs/glp1r/p6_cfg_${SLURM_ARRAY_JOB_ID:-m}/cfg${CFG}"
for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
CKPT=$(ls "$C26_WEIGHTS/rfd3"/*.ckpt | head -1)
echo "shard $SH  cfg_scale=$CFG  out $OUT  node $(hostname)"

time rfd3 design out_dir="$OUT" inputs="$G/inputs/rfd3/glp1r_ecd4.json" ckpt_path="$CKPT" \
  n_batches=1 diffusion_batch_size=2 seed=4242 \
  inference_sampler.step_scale=3 inference_sampler.gamma_0=0.2 \
  inference_sampler.cfg_scale="$CFG" \
  prevalidate_inputs=True skip_existing=False
rc=$?
echo "cfg=$CFG designs=$(find "$OUT" -name '*.cif.gz' | wc -l) rc=$rc"
echo "=== exit: $rc ==="
exit $rc
