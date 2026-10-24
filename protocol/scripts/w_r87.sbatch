#!/usr/bin/env bash
#SBATCH --job-name=w-r87
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/w_r87_%A_%a.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/w_r87_%A_%a.out
#
# Second ECD region, around Trp87. Same r3 protocol otherwise: anchor B1-10 with a free
# N-terminus, 130-residue C-terminal extension, origin inferred from the hotspots.
# Arm 0 is the single residue as requested; arm 1 adds its well-exposed neighbours,
# because a single hotspot previously gave rg_ratio 2.7 with zero contacts.
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate rfd3
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
G="$C26_ROOT/glp1r2026"
ARMS=(w1_r87_only w2_r86_89)
HOTS=(87 86,87,88,89)
S=${SLURM_ARRAY_TASK_ID:-0}
A=${ARMS[$S]}; HS=${HOTS[$S]}
OUT="$C26_SCRATCH/runs/glp1r/w_r87_${SLURM_ARRAY_JOB_ID:-m}/$A"
for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
CKPT=$(ls "$C26_WEIGHTS/rfd3"/*.ckpt | head -1)
echo "arm $A  hotspots=$HS  node $(hostname)"
time rfd3 design out_dir="$OUT" inputs="$G/inputs/rfd3/$A.json" ckpt_path="$CKPT" \
  n_batches=1 diffusion_batch_size=2 seed=$((8700 + S)) \
  inference_sampler.step_scale=3 inference_sampler.gamma_0=0.2 \
  prevalidate_inputs=True skip_existing=True
rc=$?
echo "designs: $(find "$OUT" -name '*.cif.gz' | wc -l)"
echo -n "RESULT $A  "; python "$G/scripts/r_seg_measure.py" "$OUT" "$HS" 30
echo "=== exit: $rc ==="
exit $rc
