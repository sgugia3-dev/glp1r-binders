#!/usr/bin/env bash
#SBATCH --job-name=h2-glp-mpnn
#SBATCH --partition=htc
#SBATCH --qos=htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/h2_mpnn_%j.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/h2_mpnn_%j.out
#
# CPU on purpose: htc is a separate QOS from gpus, so design costs no GPU slot.
set -o pipefail
ulimit -u 8192 2>/dev/null || true
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate rfd3          # the mpnn engine ships inside the rfd3 env

G="$C26_ROOT/glp1r2026"
IN="$(cat "$G/logs/h1_diff_outdir.txt")"
OUT="$C26_SCRATCH/runs/glp1r/h2_mpnn_${SLURM_JOB_ID}"
for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done

W="$C26_WEIGHTS/proteinmpnn/soluble_v_48_020.pt"
[ -f "$W" ] || W="$C26_ROOT/tools/mosaic/src/boltz_binder_design/proteinmpnn/weights/soluble_v_48_020.pt"
[ -f "$W" ] || { echo "FATAL: no MPNN weights"; echo "=== exit: 1 ==="; exit 1; }
echo "node $(hostname)  in $IN  out $OUT  weights $W"

python "$G/scripts/h2_mpnn.py" \
  --indir "$IN" --outdir "$OUT" --weights "$W" \
  --nseq "${NSEQ:-16}" --temperature 0.1 --seed 37 \
  --acid-pos 3 --max-ca-dev 2.0 --max-breaks 1
rc=$?
echo "$OUT" > "$G/logs/h2_mpnn_outdir.txt"
echo "=== exit: $rc ==="
exit $rc
