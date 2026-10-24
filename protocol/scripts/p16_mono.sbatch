#!/usr/bin/env bash
#SBATCH --job-name=p16-mono
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p16_mono_%j.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p16_mono_%j.out
#
# Fold every p15 binder ALONE, no receptor chain.
#
# This asks whether the designed chain holds its fold by itself or is only
# ordered because the receptor scaffolds it. A 140-residue monomer folds far
# faster than the 549-residue complex, so the whole set fits in one job.
#
# Identical binder sequences are folded once; the manifest maps each monomer back
# to every complex design that shares it.
set -o pipefail
ulimit -u 8192 2>/dev/null || true
T0=$SECONDS
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
G="$C26_ROOT/glp1r2026"
RUN=${SLURM_JOB_ID:-manual}
OUT="$C26_SCRATCH/runs/glp1r/p16_mono_${RUN}"
CPLX="$C26_SCRATCH/runs/glp1r/p15_1099477"
mkdir -p "$OUT"
echo "monomer run $RUN on $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | head -1

conda activate boltz2
export OMP_NUM_THREADS=1
Y="$OUT/yaml"
python "$G/scripts/mk_monomer_yamls.py" "$CPLX" "$Y" "$OUT/manifest.json" \
  || { echo "=== exit: 1 ==="; exit 1; }
NY=$(ls "$Y"/*.yaml | wc -l); echo "monomer yamls: $NY"

SHMC="/dev/shm/boltzcache_${SLURM_JOB_ID}"
cleanup(){ rm -rf "$SHMC" 2>/dev/null || true; }
trap cleanup EXIT
mkdir -p "$SHMC"
cp "$C26_WEIGHTS/boltz/boltz2_conf.ckpt" "$SHMC/" || { echo "FATAL ckpt"; exit 1; }
for x in mols mols.tar boltz2_aff.ckpt; do
  [ -e "$C26_WEIGHTS/boltz/$x" ] && ln -sfn "$C26_WEIGHTS/boltz/$x" "$SHMC/$x"
done
export NUMBA_CACHE_DIR="$OUT/numba"; mkdir -p "$NUMBA_CACHE_DIR"

tB=$SECONDS
boltz predict "$Y" --out_dir "$OUT" --cache "$SHMC" --model boltz2 \
  --override --accelerator gpu --devices 1 \
  --num_workers 1 --preprocessing-threads 1 --diffusion_samples 1 --seed 1
echo "FOLD rc=$? elapsed=$((SECONDS-tB))s"

P=$(find "$OUT" -maxdepth 3 -type d -name predictions | head -1)
NS=$(find "$P" -name '*_model_0.cif' 2>/dev/null | wc -l)
echo "monomer structures: $NS of $NY"
if [ "$NS" -gt 0 ]; then
  conda activate rfd3
  python "$G/scripts/measure_monomer.py" "$P" "$CPLX" "$OUT/manifest.json" \
    "$G/analysis/p16_monomer_vs_bound.csv"
fi
echo "TOTAL elapsed=$((SECONDS-T0))s"
echo "=== exit: 0 ==="
