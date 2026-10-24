#!/usr/bin/env bash
#SBATCH --job-name=p18-seed
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p18_seed%a_%A_%a.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p18_seed%a_%A_%a.out
#
# Seed-replicate re-fold of the top p15 designs.
#
# Inputs are identical to the seed-1 run; only Boltz's --seed differs. The EGFR
# submission in challenge2026/submission required every criterion to hold across
# independent folds because only 2 of 9 restrained passes survived re-folding,
# and the GLP-1R set has so far been validated at one seed only.
#
# Required: SEED (the Boltz seed) and TOPN. Array index selects the shard.
set -o pipefail
ulimit -u 8192 2>/dev/null || true
T0=$SECONDS
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
G="$C26_ROOT/glp1r2026"
SEED=${SEED:?set SEED}
TOPN=${TOPN:-150}
NSH=${NSHARDS:-10}
SH=${SLURM_ARRAY_TASK_ID:-0}
RUN="${SLURM_ARRAY_JOB_ID:-manual}"
TAG="p18s${SEED}"
OUT="$C26_SCRATCH/runs/glp1r/${TAG}_${RUN}/sh${SH}"
mkdir -p "$OUT"
echo "seed $SEED  shard $SH/$NSH  top $TOPN  on $(hostname)"

conda activate boltz2
export OMP_NUM_THREADS=1
export AF3_REQUEST="$G/fold_glp1r_glp_job_request.json"
export RECEPTOR_MSA="$G/inputs/msa/glp1r_receptor_A.a3m"
Y="$OUT/yaml"
python "$G/scripts/mk_reseed.py" "$G/analysis" \
  "$C26_SCRATCH/runs/glp1r/p15_1099477" "$TOPN" "$SH" "$NSH" \
  "$Y" "$OUT/manifest.json" || { echo "=== exit: 1 ==="; exit 1; }
NY=$(ls "$Y"/*.yaml 2>/dev/null | wc -l)
echo "yamls=$NY"
[ "$NY" -gt 0 ] || { echo "nothing for this shard"; echo "=== exit: 0 ==="; exit 0; }

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
  --write_full_pae --override --accelerator gpu --devices 1 \
  --num_workers 1 --preprocessing-threads 1 --diffusion_samples 1 --seed "$SEED"
echo "FOLD rc=$? elapsed=$((SECONDS-tB))s"

P=$(find "$OUT" -maxdepth 3 -type d -name predictions | head -1)
NS=$(find "$P" -name '*_model_0.cif' 2>/dev/null | wc -l)
echo "structures=$NS of $NY"
if [ "$NS" -gt 0 ]; then
  conda activate rfd3
  A="$G/analysis"
  python "$G/scripts/h3_measure.py" --pred-dir "$P" --manifest "$OUT/manifest.json" \
    --af3-model "$G/fold_glp1r_glp_model_0.cif" \
    --out-csv "$A/${TAG}_${RUN}_sh${SH}_deviation.csv" || true
  python "$G/scripts/h7_interactions.py" --pred-dir "$P" \
    --out-csv "$A/${TAG}_${RUN}_sh${SH}_interactions.csv" || true
  for SET in "p3:86,87,88,89" "p10:27,42,86" "p12:27,42,86" "p13:27,42,86,95"; do
    PFX="${SET%%:*}"; HOTS="${SET##*:}"
    if ls "$P" | grep -q "^${PFX}_"; then
      python "$G/scripts/h5_epitope.py" --pred-dir "$P" --manifest "$OUT/manifest.json" \
        --hotspots "$HOTS" --name-prefix "${PFX}_" \
        --dev-csv "$A/${TAG}_${RUN}_sh${SH}_deviation.csv" \
        --out-csv "$A/${TAG}_${RUN}_sh${SH}_epitope-${PFX}.csv" || true
    fi
  done
  for d in "$P"/*/; do
    nm=$(basename "$d"); cif="$d/${nm}_model_0.cif"; pae="$d/pae_${nm}_model_0.npz"
    [ -f "$cif" ] && [ -f "$pae" ] || continue
    ls "$d"/*_15_15.txt >/dev/null 2>&1 && continue
    python "$C26_ROOT/tools/ipsae/ipsae.py" "$pae" "$cif" 15 15 >/dev/null || true
  done
  echo "ipsae files: $(ls "$P"/*/*_15_15.txt 2>/dev/null | wc -l)"
fi
echo "TOTAL elapsed=$((SECONDS-T0))s"
echo "=== exit: 0 ==="
