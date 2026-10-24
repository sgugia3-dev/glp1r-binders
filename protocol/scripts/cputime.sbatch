#!/usr/bin/env bash
#SBATCH --partition=htc
#SBATCH --qos=htc
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/cputime_%j.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/cputime_%j.out
#
# How long does ONE 549-residue co-fold take on CPU?
#
# Decides whether the htc partition beats the gpus partition for the seed-replicate
# work. Both QOSes cap a user at 20 concurrent jobs, and only 2 GPUs are currently
# reachable, so 20 CPU jobs win if a fold takes less than ~18 min
# (20 jobs / T_cpu  >  2 GPUs / 110 s).
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
G="$C26_ROOT/glp1r2026"
NT=${SLURM_CPUS_PER_TASK:-16}
OUT="$C26_SCRATCH/runs/glp1r/cputime_${SLURM_JOB_ID}"
mkdir -p "$OUT/yaml"
echo "cores=$NT host=$(hostname)"

conda activate boltz2
export OMP_NUM_THREADS=$NT MKL_NUM_THREADS=$NT OPENBLAS_NUM_THREADS=$NT
export TORCH_NUM_THREADS=$NT
python -c "import torch; torch.set_num_threads($NT); print('torch threads', torch.get_num_threads())"

# one real input, copied from the seed-1 run so the system size is identical
SRC=$(ls $C26_SCRATCH/runs/glp1r/p15_1099477/sh0/yaml/*.yaml | head -1)
cp "$SRC" "$OUT/yaml/"
echo "input: $(basename $SRC)"
grep -c "sequence:" "$OUT/yaml/$(basename $SRC)" | sed "s/^/  chains: /"

SHMC="/dev/shm/bc_${SLURM_JOB_ID}"
trap 'rm -rf "$SHMC"' EXIT
mkdir -p "$SHMC"
cp "$C26_WEIGHTS/boltz/boltz2_conf.ckpt" "$SHMC/"
for x in mols mols.tar; do
  [ -e "$C26_WEIGHTS/boltz/$x" ] && ln -sfn "$C26_WEIGHTS/boltz/$x" "$SHMC/$x"
done
export NUMBA_CACHE_DIR="$OUT/numba"; mkdir -p "$NUMBA_CACHE_DIR"

t0=$SECONDS
boltz predict "$OUT/yaml" --out_dir "$OUT" --cache "$SHMC" --model boltz2 \
  --override --accelerator cpu --devices 1 \
  --num_workers 1 --preprocessing-threads 2 --diffusion_samples 1 --seed 2
rc=$?
T=$((SECONDS-t0))
N=$(find "$OUT" -name "*_model_0.cif" 2>/dev/null | wc -l)
echo "RESULT cores=$NT rc=$rc structures=$N seconds=$T"
awk -v t=$T -v c=$NT 'BEGIN{
  printf "  per fold: %.1f min\n", t/60
  printf "  20 concurrent htc jobs -> %.2f folds/hour\n", 20*3600/t
  printf "  2 reachable GPUs at 110 s -> %.2f folds/hour\n", 2*3600/110
  printf "  300 folds on htc: %.1f h   on 2 GPUs: %.1f h\n", 300*t/20/3600, 300*110/2/3600
  if (t < 1100) print "  VERDICT: htc is faster"; else print "  VERDICT: gpus is faster"
}'
echo "=== exit: $rc ==="
