#!/usr/bin/env bash
#SBATCH --job-name=s8-hotcount
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/s8_hotcount/s8_%j.out
#SBATCH --error=/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/s8_hotcount/s8_%j.out
#
# STRATEGY: hotspot-count. Does imposing FEWER hotspots let the binder compact while
# still contacting? Everything else fixed: binder_ROG weight:1, guide_scale=10, quadratic.
# Four configs run SEQUENTIALLY in this one job (20-job cap, 8 agents concurrent).
#   hot4      [A24,A25,A27,A42]  <- reproduces the known hot+ROG_w1 arm (rg 1.46, 0/4)
#   hot2      [A24,A42]          <- the two extremes of the patch, ~12 A apart
#   hot1_A42  [A42]              <- single anchor
#   hot1_A24  [A24]              <- single anchor, other end
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/anaconda1/etc/profile.d/conda.sh
conda activate SE3nv

G=/lustre/home/gsgueglia/challenge2026/glp1r2026
R=/lustre/home/gsgueglia/diff_test/RFdiffusion
BASE=/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/s8_hotcount
NDES=${NDES:-6}
BLEN=${BLEN:-70}
PDB="$G/inputs/rfd2/glp1r_ecd_A1_117.pdb"

echo "node $(hostname)  binder_len=$BLEN  n=$NDES  job=$SLURM_JOB_ID"
nvidia-smi --query-gpu=driver_version --format=csv,noheader | head -1
python -c "import torch;print('torch',torch.__version__,'cuda avail',torch.cuda.is_available())" || exit 1
cd "$R" || exit 2

run_cfg () {
  local LBL="$1" HOT="$2"
  local OUT="$BASE/$LBL"
  for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
  echo "=============== CONFIG $LBL  hotspots=$HOT ==============="
  local ARGS=(
    "inference.input_pdb=$PDB"
    "inference.output_prefix=$OUT/des"
    "inference.num_designs=$NDES"
    "contigmap.contigs=[A1-117/0 ${BLEN}-${BLEN}]"
    "ppi.hotspot_res=$HOT"
    'potentials.guiding_potentials=["type:binder_ROG,weight:1"]'
    "potentials.guide_scale=10"
    "potentials.guide_decay=quadratic"
  )
  echo "FLAGS[$LBL]: ${ARGS[*]}"
  time python scripts/run_inference.py "${ARGS[@]}"
  local rc=$?
  echo "config $LBL designs=$(ls "$OUT"/des_*.pdb 2>/dev/null | wc -l) rc=$rc"
  return 0
}

run_cfg hot4     '[A24,A25,A27,A42]'
run_cfg hot2     '[A24,A42]'
run_cfg hot1_A42 '[A42]'
run_cfg hot1_A24 '[A24]'

echo "################ MEASUREMENT ################"
for LBL in hot4 hot2 hot1_A42 hot1_A24; do
  echo "---- MEASURE $LBL (default hotspot set 24,25,27,42) ----"
  python "$G/scripts/pot_measure.py" "$BASE/$LBL"
done
echo "################ SUPPLEMENTARY: measured against each config's OWN imposed subset ################"
echo "---- $LBL=hot2 own subset 24,42 ----";    python "$G/scripts/pot_measure.py" "$BASE/hot2" 24,42
echo "---- hot1_A42 own subset 42 ----";        python "$G/scripts/pot_measure.py" "$BASE/hot1_A42" 42
echo "---- hot1_A24 own subset 24 ----";        python "$G/scripts/pot_measure.py" "$BASE/hot1_A24" 24
echo "=== exit: 0 ==="
exit 0
