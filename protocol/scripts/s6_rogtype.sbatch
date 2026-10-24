#!/usr/bin/env bash
#SBATCH --job-name=s6-rogtype
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/s6_rogtype_%j.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/s6_rogtype_%j.out
#
# STRATEGY s6: monomer_ROG vs binder_ROG at matched weight, with the 4 hotspots.
# Source fact (rfdiffusion/potentials/potentials.py): monomer_ROG.compute uses xyz[:,1]
# (ALL L residues) and monomer_ROG is NOT in require_binderlen (line 469), so it is
# constructed with weight only and measures the Rg of the WHOLE 187-res complex, not the
# 70-res binder. Combined with Ca_grads[diffusion_mask]=0 (inference/utils.py:397) the
# gradient lands only on the binder, so monomer_ROG here = "drag binder toward the
# binder+target centroid". Four configs, SEQUENTIAL in one job.
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/anaconda1/etc/profile.d/conda.sh
conda activate SE3nv

G=/lustre/home/gsgueglia/challenge2026/glp1r2026
R=/lustre/home/gsgueglia/diff_test/RFdiffusion
PDB="$G/inputs/rfd2/glp1r_ecd_A1_117.pdb"
BASE=/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/s6_rogtype
NDES=${NDES:-6}
BLEN=70
HOT='[A24,A25,A27,A42]'

echo "node $(hostname)  job ${SLURM_JOB_ID}  binder_len=$BLEN  n=$NDES"
nvidia-smi --query-gpu=driver_version --format=csv,noheader | head -1
python -c "import torch;print('torch',torch.__version__,'cuda avail',torch.cuda.is_available())" || exit 1
cd "$R" || exit 2

run_cfg () {
  local NAME="$1" POT="$2" GS="$3"
  local OUT="$BASE/$NAME"
  for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
  echo "##### CONFIG $NAME  pot=$POT  guide_scale=$GS  out=$OUT"
  local ARGS=(
    "inference.input_pdb=$PDB"
    "inference.output_prefix=$OUT/des"
    "inference.num_designs=$NDES"
    "contigmap.contigs=[A1-117/0 ${BLEN}-${BLEN}]"
    "ppi.hotspot_res=$HOT"
    "potentials.guiding_potentials=[\"type:$POT\"]"
    "potentials.guide_scale=$GS"
    "potentials.guide_decay=quadratic"
  )
  echo "FLAGS $NAME: ${ARGS[*]}"
  time python scripts/run_inference.py "${ARGS[@]}"
  local rc=$?
  echo "##### CONFIG $NAME done designs=$(ls "$OUT"/des_*.pdb 2>/dev/null | wc -l) rc=$rc"
}

run_cfg c1_mono_w1        "monomer_ROG,weight:1" 10
run_cfg c2_mono_w5        "monomer_ROG,weight:5" 10
run_cfg c3_binder_w1_ref  "binder_ROG,weight:1"  10
run_cfg c4_mono_w1_gs5    "monomer_ROG,weight:1" 5

echo "===== MEASUREMENTS ====="
for N in c1_mono_w1 c2_mono_w5 c3_binder_w1_ref c4_mono_w1_gs5; do
  echo "--- $N"
  python "$G/scripts/pot_measure.py" "$BASE/$N" || echo "measure failed $N"
done
echo "=== ALL DONE ==="
