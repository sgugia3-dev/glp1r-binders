#!/usr/bin/env bash
#SBATCH --job-name=s3-mindist
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/s3_mindist_%j.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/s3_mindist_%j.out
#
# STRATEGY rog-mindist: vary binder_ROG min_dist at FIXED weight 1, guide_scale 10,
# quadratic decay, 4 hotspots. binder_ROG clamps the per-residue centroid distance with
# max(min_dist, d), so residues already inside min_dist contribute a constant and get ZERO
# gradient. min_dist therefore sets how compact "compact enough" is:
#   min_dist 8  -> keeps pulling until every CA is within 8 A of the centroid (most compact)
#   min_dist 20 -> essentially every CA is already clamped -> potential inert (expect baseline)
# Four configs, SEQUENTIAL in one job. Everything else identical.
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/anaconda1/etc/profile.d/conda.sh
conda activate SE3nv

G=/lustre/home/gsgueglia/challenge2026/glp1r2026
R=/lustre/home/gsgueglia/diff_test/RFdiffusion
BASE=/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/s3_mindist
NDES=${NDES:-6}
BLEN=${BLEN:-70}
PDB="$G/inputs/rfd2/glp1r_ecd_A1_117.pdb"
HOT='[A24,A25,A27,A42]'

echo "host $(hostname)  job ${SLURM_JOB_ID}  binder_len=$BLEN  n=$NDES per config"
nvidia-smi --query-gpu=driver_version,name --format=csv,noheader | head -1
python -c "import torch;print('torch',torch.__version__,'cuda avail',torch.cuda.is_available())" || exit 1
cd "$R" || exit 2

overall=0
for MD in 8 12 15 20; do
  LBL=$(printf "md%02d" "$MD")
  OUT="$BASE/$LBL"
  for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
  POT="type:binder_ROG,weight:1,min_dist:${MD}"
  echo "=========================================================="
  echo "### CONFIG $LBL  potential=[$POT]  guide_scale=10 quadratic  hotspots=$HOT"
  echo "### outdir=$OUT  start $(date -Is)"
  time python scripts/run_inference.py \
    "inference.input_pdb=$PDB" \
    "inference.output_prefix=$OUT/des" \
    "inference.num_designs=$NDES" \
    "contigmap.contigs=[A1-117/0 ${BLEN}-${BLEN}]" \
    "ppi.hotspot_res=$HOT" \
    "potentials.guiding_potentials=[\"$POT\"]" \
    potentials.guide_scale=10 potentials.guide_decay=quadratic
  rc=$?
  n=$(ls "$OUT"/des_*.pdb 2>/dev/null | wc -l)
  echo "### CONFIG $LBL done rc=$rc designs=$n  end $(date -Is)"
  [ "$rc" -ne 0 ] && overall=$rc
done

echo "=========================================================="
for MD in 8 12 15 20; do
  LBL=$(printf "md%02d" "$MD")
  echo "--- MEASURE $LBL ---"
  python "$G/scripts/pot_measure.py" "$BASE/$LBL" 2>&1
done
echo "=== exit: $overall ==="
exit $overall
