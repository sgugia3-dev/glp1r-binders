#!/usr/bin/env bash
#SBATCH --job-name=q1-free
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/q1_free_%A_%a.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/q1_free_%A_%a.out
#
# UNANCHORED binder to the hotspot patch: no peptide motif, no C-terminal extension.
# Task 0 = rfd3 (no geometric potential available) -- the missing comparator: does rfd3
#          produce compact binders once the motif/hotspot conflict is removed?
# Task 1 = classic RFdiffusion + binder_ROG at the settings that looked least bad so far
#          (guide_scale 10, min_dist 12, weight 1), with n raised for statistical power
#          since every arm so far is only 6 designs.
set -o pipefail
ulimit -u 8192 2>/dev/null || true
G=/lustre/home/gsgueglia/challenge2026/glp1r2026
S=${SLURM_ARRAY_TASK_ID:-0}
OUT=/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/q1_free_${SLURM_ARRAY_JOB_ID:-m}/t$S
for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
echo "task $S on $(hostname)  out $OUT"

if [ "$S" = "0" ]; then
  source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
  conda activate rfd3
  export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
  CKPT=$(ls "$C26_WEIGHTS/rfd3"/*.ckpt | head -1)
  time rfd3 design out_dir="$OUT" inputs="$G/inputs/rfd3/glp1r_free.json" ckpt_path="$CKPT" \
    n_batches=1 diffusion_batch_size=12 seed=909 \
    inference_sampler.step_scale=3 inference_sampler.gamma_0=0.2 \
    prevalidate_inputs=True skip_existing=False
  rc=$?
  echo "rfd3 unanchored: $(find "$OUT" -name '*.cif.gz' | wc -l) designs rc=$rc"
else
  source /lustre/home/gsgueglia/anaconda1/etc/profile.d/conda.sh
  conda activate SE3nv
  cd /lustre/home/gsgueglia/diff_test/RFdiffusion || exit 2
  time python scripts/run_inference.py \
    inference.input_pdb="$G/inputs/rfd2/glp1r_ecd_A1_117.pdb" \
    inference.output_prefix="$OUT/des" \
    inference.num_designs=24 \
    'contigmap.contigs=[A1-117/0 70-70]' \
    'ppi.hotspot_res=[A24,A25,A27,A42]' \
    'potentials.guiding_potentials=["type:binder_ROG,weight:1,min_dist:12"]' \
    potentials.guide_scale=10 potentials.guide_decay=quadratic
  rc=$?
  echo "classic unanchored: $(ls "$OUT"/des_*.pdb 2>/dev/null | wc -l) designs rc=$rc"
  python "$G/scripts/pot_measure.py" "$OUT" || true
fi
echo "=== exit: $rc ==="
exit $rc
