#!/usr/bin/env bash
#SBATCH --job-name=p15-anc
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p15_anc_%A_%a.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p15_anc_%A_%a.out
#
# Anchor repair under the agonist consensus, on the best low-charge p14 binders.
#
# Positions 11-140 are held fixed, so the charge-complementary surface that round
# produced is preserved and only the 10-residue agonist pharmacophore changes.
# Measured justification: 0/600 p14 anchors satisfied the consensus and 100% broke
# the invariant Gly4 and Thr5 (scripts/consensus.py).
#
# Restraints per position, enforced by omission:
#   1 H/Y/F   2 G/A/S   3 E/D   4 G*  5 T*  6 F*  7 T/I   8 S*  9 D*  10 Y/V/L/K
#   (* invariant at 100% across the 23 verified agonists)
set -o pipefail
ulimit -u 8192 2>/dev/null || true
T0=$SECONDS
DEADLINE_S=${DEADLINE_S:-12600}
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
G="$C26_ROOT/glp1r2026"
SH=${SLURM_ARRAY_TASK_ID:-0}
NSH=${NSHARDS:-10}
RUN="${SLURM_ARRAY_JOB_ID:-manual}"
OUT="$C26_SCRATCH/runs/glp1r/p15_${RUN}/sh${SH}"
for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
NSEQ=${NSEQ:-16}; KEEP=${KEEP:-3}
echo "shard $SH/$NSH on $(hostname)"

conda activate rfd3
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LIST="$OUT/inputs.txt"
python "$G/scripts/pick_inputs.py" "$G/analysis" \
  "$C26_SCRATCH/runs/glp1r/p14_1092297" "$LIST"
mapfile -t ALL < "$LIST"
MINE=()
for ((i=SH; i<${#ALL[@]}; i+=NSH)); do MINE+=("${ALL[$i]}"); done
echo "this shard: ${#MINE[@]} of ${#ALL[@]} backbones"
[ ${#MINE[@]} -gt 0 ] || { echo "nothing"; echo "=== exit: 0 ==="; exit 0; }

W="$C26_WEIGHTS/proteinmpnn/soluble_v_48_020.pt"
[ -f "$W" ] || W="$C26_ROOT/tools/mosaic/src/boltz_binder_design/proteinmpnn/weights/soluble_v_48_020.pt"
MP="$OUT/mpnn"; mkdir -p "$MP"
tA=$SECONDS
python "$G/scripts/h8_anchor.py" --complexes "${MINE[@]}" --outdir "$MP" \
  --weights "$W" --nseq "$NSEQ" --keep "$KEEP" --temperature 0.2 \
  --seed $((900 + SH)) --top-bias 1.0
rc=$?
echo "STAGE_A rc=$rc elapsed=$((SECONDS-tA))s"
[ $rc -eq 0 ] || { echo "FATAL mpnn"; echo "=== exit: $rc ==="; exit $rc; }

conda activate boltz2
Y="$OUT/yaml"; mkdir -p "$Y"
python "$G/scripts/h3_make_yamls.py" --mpnn-json "$MP/h8_mpnn.json" \
  --af3-request "$G/fold_glp1r_glp_job_request.json" --out-dir "$Y" \
  --per-backbone "$KEEP" \
  --receptor-msa "$G/inputs/msa/glp1r_receptor_A.a3m" \
  --manifest "$OUT/manifest.json" || { echo "=== exit: 1 ==="; exit 1; }
NY=$(ls "$Y"/*.yaml | wc -l); echo "yamls=$NY"

SHMC="/dev/shm/boltzcache_${SLURM_JOB_ID}"
cleanup(){ rm -rf "$SHMC" 2>/dev/null || true; }
trap cleanup EXIT
mkdir -p "$SHMC"
cp "$C26_WEIGHTS/boltz/boltz2_conf.ckpt" "$SHMC/" || { echo "FATAL ckpt"; echo "=== exit: 1 ==="; exit 1; }
for x in mols mols.tar boltz2_aff.ckpt; do
  [ -e "$C26_WEIGHTS/boltz/$x" ] && ln -sfn "$C26_WEIGHTS/boltz/$x" "$SHMC/$x"
done
export NUMBA_CACHE_DIR="$OUT/numba"; mkdir -p "$NUMBA_CACHE_DIR"
LEFT=$((DEADLINE_S - (SECONDS - T0)))
echo "STAGE_B budget_left=${LEFT}s for $NY folds"
if [ "$LEFT" -lt 300 ]; then echo "WARN: skipping folds"; else
  tB=$SECONDS
  timeout "${LEFT}s" boltz predict "$Y" --out_dir "$OUT" --cache "$SHMC" --model boltz2 \
    --write_full_pae --override --accelerator gpu --devices 1 \
    --num_workers 1 --preprocessing-threads 1 --diffusion_samples 1 --seed 1
  echo "STAGE_B rc=$? elapsed=$((SECONDS-tB))s"
fi

P=$(find "$OUT" -maxdepth 3 -type d -name predictions | head -1)
NS=$(find "$P" -name '*_model_0.cif' 2>/dev/null | wc -l)
echo "STAGE_C structures=$NS of $NY"
if [ "$NS" -gt 0 ]; then
  conda activate rfd3
  A="$G/analysis"
  python "$G/scripts/h3_measure.py" --pred-dir "$P" --manifest "$OUT/manifest.json" \
    --af3-model "$G/fold_glp1r_glp_model_0.cif" \
    --out-csv "$A/p15_${RUN}_sh${SH}_deviation.csv" || true
  python "$G/scripts/h7_interactions.py" --pred-dir "$P" \
    --out-csv "$A/p15_${RUN}_sh${SH}_interactions.csv" || true
  for SET in "p3:86,87,88,89" "p10:27,42,86" "p12:27,42,86" "p13:27,42,86,95"; do
    PFX="${SET%%:*}"; HOTS="${SET##*:}"
    if ls "$P" | grep -q "^${PFX}_"; then
      python "$G/scripts/h5_epitope.py" --pred-dir "$P" --manifest "$OUT/manifest.json" \
        --hotspots "$HOTS" --name-prefix "${PFX}_" \
        --dev-csv "$A/p15_${RUN}_sh${SH}_deviation.csv" \
        --out-csv "$A/p15_${RUN}_sh${SH}_epitope-${PFX}.csv" || true
    fi
  done
fi
echo "TOTAL elapsed=$((SECONDS-T0))s"
echo "=== exit: 0 ==="
exit 0
