#!/usr/bin/env bash
#SBATCH --job-name=s1-gscale
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/s1_guidescale/s1_%j.out
#SBATCH --error=/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/s1_guidescale/s1_%j.out
#
# STRATEGY s1: guide_scale sweep at FIXED binder_ROG weight:1 with the 4 hotspots.
# Hypothesis: weight 5 overshot (rg_ratio 2.01) because the guidance STEP is too large,
# not because the potential is wrong. guide_scale multiplies the step (manager.py:180
# get_guide_scale -> quadratic: t^2/T^2 * guide_scale), so scale is the knob to sweep.
# Four configs run SEQUENTIALLY in this one job (20-job cap, 8 agents concurrent).
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/anaconda1/etc/profile.d/conda.sh
conda activate SE3nv

G=/lustre/home/gsgueglia/challenge2026/glp1r2026
R=/lustre/home/gsgueglia/diff_test/RFdiffusion
BASE=/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/s1_guidescale
NDES=${NDES:-6}
BLEN=${BLEN:-70}
PDB="$G/inputs/rfd2/glp1r_ecd_A1_117.pdb"
HOT='[A24,A25,A27,A42]'

echo "node $(hostname)  binder_len=$BLEN  n=$NDES per config"
nvidia-smi --query-gpu=driver_version,name --format=csv,noheader | head -1
python -c "import torch;print('torch',torch.__version__,'cuda avail',torch.cuda.is_available())" || exit 1
cd "$R" || exit 2

rc_all=0
for GS in 2 5 20 40; do
  OUT="$BASE/gs$GS"
  for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
  echo "=== CONFIG gs$GS : binder_ROG weight:1, guide_scale=$GS, guide_decay=quadratic, hotspots=$HOT ==="
  time python scripts/run_inference.py \
    "inference.input_pdb=$PDB" \
    "inference.output_prefix=$OUT/des" \
    "inference.num_designs=$NDES" \
    "contigmap.contigs=[A1-117/0 ${BLEN}-${BLEN}]" \
    "ppi.hotspot_res=$HOT" \
    'potentials.guiding_potentials=["type:binder_ROG,weight:1"]' \
    "potentials.guide_scale=$GS" \
    potentials.guide_decay=quadratic
  rc=$?
  n=$(ls "$OUT"/des_*.pdb 2>/dev/null | wc -l)
  echo "=== CONFIG gs$GS done: designs=$n rc=$rc ==="
  [ "$rc" -ne 0 ] && rc_all=$rc
done

echo "=== ALL CONFIGS DONE rc_all=$rc_all ==="
for GS in 2 5 20 40; do
  echo "--- MEASURE gs$GS ---"
  python "$G/scripts/pot_measure.py" "$BASE/gs$GS"
done
echo "=== exit: $rc_all ==="
exit $rc_all
