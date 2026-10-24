#!/usr/bin/env bash
#SBATCH --job-name=s4-iface
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/s4_iface_%j.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/s4_iface_%j.out
#
# STRATEGY 4: express the contact objective as a guiding POTENTIAL
# (interface_ncontacts, defaults r_0=8 d_0=6) instead of relying only on soft
# hotspot conditioning, combined with binder_ROG for compactness.
# Hotspots stay set in every arm. Four configs run SEQUENTIALLY in this one job.
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/anaconda1/etc/profile.d/conda.sh
conda activate SE3nv

G=/lustre/home/gsgueglia/challenge2026/glp1r2026
R=/lustre/home/gsgueglia/diff_test/RFdiffusion
BASE=/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/s4_ifacecontacts
NDES=${NDES:-6}
BLEN=${BLEN:-70}
PDB="$G/inputs/rfd2/glp1r_ecd_A1_117.pdb"
HOT='[A24,A25,A27,A42]'

echo "node $(hostname)  binder_len=$BLEN  n=$NDES"
nvidia-smi --query-gpu=driver_version --format=csv,noheader | head -1
python -c "import torch;print('torch',torch.__version__,'cuda avail',torch.cuda.is_available())" || exit 1

run_cfg () {
  local LBL="$1"; local POTLIST="$2"
  local OUT="$BASE/$LBL"
  for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
  echo "=== CONFIG $LBL ==="
  echo "FLAGS: potentials.guiding_potentials=$POTLIST potentials.guide_scale=10 potentials.guide_decay=quadratic ppi.hotspot_res=$HOT"
  cd "$R" || return 2
  time python scripts/run_inference.py \
    "inference.input_pdb=$PDB" \
    "inference.output_prefix=$OUT/des" \
    "inference.num_designs=$NDES" \
    "contigmap.contigs=[A1-117/0 ${BLEN}-${BLEN}]" \
    "ppi.hotspot_res=$HOT" \
    "potentials.guiding_potentials=$POTLIST" \
    potentials.guide_scale=10 potentials.guide_decay=quadratic
  local rc=$?
  echo "CONFIG $LBL designs=$(ls "$OUT"/des_*.pdb 2>/dev/null | wc -l) rc=$rc"
  return 0
}

run_cfg rogw1_ifw1  '["type:binder_ROG,weight:1","type:interface_ncontacts,weight:1"]'
run_cfg rogw1_ifw3  '["type:binder_ROG,weight:1","type:interface_ncontacts,weight:3"]'
run_cfg rogw1_ifw10 '["type:binder_ROG,weight:1","type:interface_ncontacts,weight:10"]'
run_cfg ifw3_noROG  '["type:interface_ncontacts,weight:3"]'

echo "=== MEASUREMENTS ==="
for L in rogw1_ifw1 rogw1_ifw3 rogw1_ifw10 ifw3_noROG; do
  echo "--- $L"
  python "$G/scripts/pot_measure.py" "$BASE/$L" 2>&1 | tail -12
done
echo "=== ALL DONE ==="
