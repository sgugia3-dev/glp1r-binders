#!/usr/bin/env bash
#SBATCH --job-name=h3-glp-fold
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/h3_fold_%j.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/h3_fold_%j.out
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate boltz2

G="$C26_ROOT/glp1r2026"
MP="$(cat "$G/logs/h2_mpnn_outdir.txt")/h2_mpnn.json"
OUT="$C26_SCRATCH/runs/glp1r/h3_fold_${SLURM_JOB_ID}"
YAML="$OUT/yaml"
for i in 1 2 3 4 5; do mkdir -p "$YAML" && break || { echo "mkdir retry $i"; sleep 20; }; done

echo "node $(hostname)  mpnn $MP  out $OUT"
python "$G/scripts/h3_make_yamls.py" \
  --mpnn-json "$MP" \
  --af3-request "$G/fold_glp1r_glp_job_request.json" \
  --out-dir "$YAML" \
  --per-backbone "${PERBB:-3}" \
  --receptor-msa "${RMSA:-$G/inputs/msa/glp1r_receptor_A.a3m}" \
  --manifest "$OUT/manifest.json" || { echo "=== exit: 1 ==="; exit 1; }
N=$(ls "$YAML"/*.yaml | wc -l); echo "yamls: $N"
[ "$N" -gt 0 ] || { echo "FATAL: no yamls"; echo "=== exit: 1 ==="; exit 1; }

# stage the checkpoint into node RAM: reading it off Lustre from many slots caused
# ESHUTDOWN and hung logins before. One read per node, cleaned up on exit.
SHMC="/dev/shm/boltzcache_${SLURM_JOB_ID}"
cleanup(){ rm -rf "$SHMC" 2>/dev/null || true; }
trap cleanup EXIT
mkdir -p "$SHMC"
cp "$C26_WEIGHTS/boltz/boltz2_conf.ckpt" "$SHMC/" || { echo "FATAL: ckpt stage"; echo "=== exit: 1 ==="; exit 1; }
for x in mols mols.tar boltz2_aff.ckpt; do
  [ -e "$C26_WEIGHTS/boltz/$x" ] && ln -sfn "$C26_WEIGHTS/boltz/$x" "$SHMC/$x"
done
export NUMBA_CACHE_DIR="$OUT/numba"; mkdir -p "$NUMBA_CACHE_DIR"

t0=$SECONDS
boltz predict "$YAML" --out_dir "$OUT" --cache "$SHMC" --model boltz2 \
  --write_full_pae --override --accelerator gpu --devices 1 \
  --num_workers 2 --preprocessing-threads 2 --diffusion_samples 1 --seed 1
rc=$?
dt=$((SECONDS-t0))
echo "boltz rc=$rc wall ${dt}s for $N folds ($((dt/(N>0?N:1)))s each)"

nf=$(grep -c "Failed to process" "$G/logs/h3_fold_${SLURM_JOB_ID}.out" 2>/dev/null || true)
echo "failed-to-process records: $nf"
MAN=$(find "$OUT" -name manifest.json -path "*boltz*" | head -1)
P=$(find "$OUT" -maxdepth 3 -type d -name predictions | head -1)
echo "predictions dir: $P  structures: $(find "$P" -name '*_model_0.cif' 2>/dev/null | wc -l)"

conda activate rfd3 2>/dev/null || true
python "$G/scripts/h3_measure.py" \
  --pred-dir "$P" --manifest "$OUT/manifest.json" \
  --af3-model "$G/fold_glp1r_glp_model_0.cif" \
  --out-csv "$G/analysis/h3_deviation_${SLURM_JOB_ID}.csv"
mrc=$?
echo "$OUT" > "$G/logs/h3_fold_outdir.txt"
echo "=== exit: $((rc||mrc)) ==="
exit $((rc||mrc))
