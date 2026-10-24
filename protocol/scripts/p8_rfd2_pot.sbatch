#!/usr/bin/env bash
#SBATCH --job-name=p8-rfd2pot
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p8_rfd2_%A_%a.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p8_rfd2_%A_%a.out
#
# Do classic RFdiffusion's guiding potentials enforce hotspot binding AND a compact binder
# at the same time? rfd3 has no geometric potential; this checkout does
# (rfdiffusion/potentials/potentials.py: monomer_ROG, binder_ROG, binder_ncontacts,
# interface_ncontacts, ...). binder_ROG is the right one here because the system is a
# binder + target complex, and it clamps at min_dist so it stops pulling once compact.
#
# Four arms, one per array task, everything else identical:
#   0  hotspots, no potential            <- baseline, shows the untreated Rg
#   1  hotspots + binder_ROG weight 1
#   2  hotspots + binder_ROG weight 5
#   3  binder_ROG weight 5, NO hotspots  <- control: does the potential do anything at all
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/anaconda1/etc/profile.d/conda.sh
conda activate SE3nv

G=/lustre/home/gsgueglia/challenge2026/glp1r2026
R=/lustre/home/gsgueglia/diff_test/RFdiffusion
S=${SLURM_ARRAY_TASK_ID:-0}
NDES=${NDES:-6}
BLEN=${BLEN:-70}
PDB="$G/inputs/rfd2/glp1r_ecd_A1_117.pdb"
OUT="/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/p8_rfd2_${SLURM_ARRAY_JOB_ID:-m}/arm$S"
for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done

case $S in
  0) POT=""                   ; HOT='[A24,A25,A27,A42]' ; LBL="hotspots_only" ;;
  1) POT="binder_ROG,weight:1"; HOT='[A24,A25,A27,A42]' ; LBL="hot+ROG_w1" ;;
  2) POT="binder_ROG,weight:5"; HOT='[A24,A25,A27,A42]' ; LBL="hot+ROG_w5" ;;
  3) POT="binder_ROG,weight:5"; HOT=''                  ; LBL="ROG_w5_noHot" ;;
esac
echo "arm $S ($LBL)  node $(hostname)  binder_len=$BLEN  n=$NDES"
nvidia-smi --query-gpu=driver_version --format=csv,noheader | head -1
python -c "import torch;print('torch',torch.__version__,'cuda avail',torch.cuda.is_available())" || exit 1

ARGS=(
  "inference.input_pdb=$PDB"
  "inference.output_prefix=$OUT/des"
  "inference.num_designs=$NDES"
  "contigmap.contigs=[A1-117/0 ${BLEN}-${BLEN}]"
)
[ -n "$HOT" ] && ARGS+=( "ppi.hotspot_res=$HOT" )
if [ -n "$POT" ]; then
  ARGS+=( "potentials.guiding_potentials=[\"type:$POT\"]"
          "potentials.guide_scale=10" "potentials.guide_decay=quadratic" )
fi
echo "ARGS: ${ARGS[*]}"
cd "$R" || exit 2
time python scripts/run_inference.py "${ARGS[@]}"
rc=$?
echo "arm $S ($LBL) designs=$(ls "$OUT"/des_*.pdb 2>/dev/null | wc -l) rc=$rc"
echo "=== exit: $rc ==="
exit $rc
