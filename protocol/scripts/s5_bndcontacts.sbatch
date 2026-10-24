#!/usr/bin/env bash
#SBATCH --job-name=s5-bndcon
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/s5_bndcon_%j.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/s5_bndcon_%j.out
#
# STRATEGY 5: binder_ncontacts as the compactness driver instead of binder_ROG.
# binder_ncontacts (potentials.py:118) maximises the smooth count of intra-binder Ca-Ca
# contacts, so it rewards the binder folding onto itself rather than penalising the
# distance of every residue from the binder centroid. Hypothesis: a contact-count reward
# saturates once the chain is packed and so stops fighting the interface, whereas
# binder_ROG keeps pulling and at weight 5 drags the binder away from the target.
#
# Four configs, run SEQUENTIALLY in this one job (20-job cap, 8 agents):
#   c1  hotspots + binder_ncontacts w1
#   c2  hotspots + binder_ncontacts w5
#   c3  hotspots + binder_ncontacts w1 + binder_ROG w1   (do they stack?)
#   c4  binder_ncontacts w5, NO hotspots                 (control: effect of the potential alone)
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/anaconda1/etc/profile.d/conda.sh
conda activate SE3nv

G=/lustre/home/gsgueglia/challenge2026/glp1r2026
R=/lustre/home/gsgueglia/diff_test/RFdiffusion
BASE=/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/s5_bndcontacts
NDES=${NDES:-6}
BLEN=${BLEN:-70}
PDB="$G/inputs/rfd2/glp1r_ecd_A1_117.pdb"
HOTS='[A24,A25,A27,A42]'

echo "host $(hostname)  job ${SLURM_JOB_ID}  binder_len=$BLEN  n=$NDES"
nvidia-smi --query-gpu=driver_version,name --format=csv,noheader | head -1
python -c "import torch;print('torch',torch.__version__,'cuda avail',torch.cuda.is_available())" || exit 1
cd "$R" || exit 2

run_cfg () {
  local LBL="$1" HOT="$2" POTLIST="$3"
  local OUT="$BASE/$LBL"
  for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
  local ARGS=(
    "inference.input_pdb=$PDB"
    "inference.output_prefix=$OUT/des"
    "inference.num_designs=$NDES"
    "contigmap.contigs=[A1-117/0 ${BLEN}-${BLEN}]"
    "potentials.guiding_potentials=$POTLIST"
    "potentials.guide_scale=10"
    "potentials.guide_decay=quadratic"
  )
  [ -n "$HOT" ] && ARGS+=( "ppi.hotspot_res=$HOT" )
  echo "##### CONFIG $LBL"
  echo "##### FLAGS: ${ARGS[*]}"
  time python scripts/run_inference.py "${ARGS[@]}"
  local rc=$?
  echo "##### $LBL designs=$(ls "$OUT"/des_*.pdb 2>/dev/null | wc -l) rc=$rc"
}

run_cfg "c1_bndcon_w1"        "$HOTS" '["type:binder_ncontacts,weight:1"]'
run_cfg "c2_bndcon_w5"        "$HOTS" '["type:binder_ncontacts,weight:5"]'
run_cfg "c3_bndcon_w1_rog_w1" "$HOTS" '["type:binder_ncontacts,weight:1","type:binder_ROG,weight:1"]'
run_cfg "c4_bndcon_w5_noHot"  ""      '["type:binder_ncontacts,weight:5"]'

echo "=== ALL CONFIGS DONE ==="
for d in c1_bndcon_w1 c2_bndcon_w5 c3_bndcon_w1_rog_w1 c4_bndcon_w5_noHot; do
  echo "--- $d : $(ls $BASE/$d/des_*.pdb 2>/dev/null | wc -l) pdbs"
done
exit 0
