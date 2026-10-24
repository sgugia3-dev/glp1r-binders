#!/usr/bin/env bash
#SBATCH --job-name=p14-chg
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p14_chg_%A_%a.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p14_chg_%A_%a.out
#
# Charge-complementary redesign round on the validated backbones from P3, P10,
# P12 and P13 (P2 excluded as instructed). 200 input complexes.
#
# What changes vs every previous round:
#   - the 10-residue anchor is DESIGNED, not fixed: position 1 forced to HIS and
#     position 3 to ASP/GLU by omission; the other eight are free
#   - residues 11-140 get a per-position bias derived from their measured
#     receptor partners, so charge complementarity, salt bridges, H-bonds and
#     hydrophobic contacts are all targeted at the positions that actually touch
#     the receptor, leaving interior positions unbiased so the fold survives
#
# Grounding (scripts/contact_charge.py over the 200 inputs):
#   binder 11-140 contacts ECD: 44 res, net -3.0  -> designed region biased BASIC
#   binder 1-10   contacts TMD: 27 res, net +3.0  -> anchor biased ACIDIC
#
# Epitope scoring uses each backbone's OWN hotspot set, taken from its filename
# prefix, because the inputs come from four campaigns with three different sets.
set -o pipefail
ulimit -u 8192 2>/dev/null || true
T0=$SECONDS
DEADLINE_S=${DEADLINE_S:-12600}

source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
G="$C26_ROOT/glp1r2026"
SH=${SLURM_ARRAY_TASK_ID:-0}
NSH=${NSHARDS:-20}
RUN="${SLURM_ARRAY_JOB_ID:-manual}"
OUT="$C26_SCRATCH/runs/glp1r/p14_${RUN}/sh${SH}"
for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
NSEQ=${NSEQ:-32}; KEEP=${KEEP:-3}
echo "shard $SH/$NSH on $(hostname)  nseq=$NSEQ keep=$KEEP"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | head -1

# ---------- shard the 200 input complexes deterministically ----------
mapfile -t ALL < <(ls $G/p3_best/structures/*.cif $G/p10_best/structures/*.cif \
                      $G/p12_best/structures/*.cif $G/p13_best/structures/*.cif 2>/dev/null | sort)
echo "total input complexes: ${#ALL[@]}"
MINE=()
for ((i=SH; i<${#ALL[@]}; i+=NSH)); do MINE+=("${ALL[$i]}"); done
echo "this shard designs ${#MINE[@]}"
[ ${#MINE[@]} -gt 0 ] || { echo "nothing for this shard"; echo "=== exit: 0 ==="; exit 0; }

# ---------- A: charge-aware MPNN (CPU) ----------
conda activate rfd3
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W="$C26_WEIGHTS/proteinmpnn/soluble_v_48_020.pt"
[ -f "$W" ] || W="$C26_ROOT/tools/mosaic/src/boltz_binder_design/proteinmpnn/weights/soluble_v_48_020.pt"
MP="$OUT/mpnn"; mkdir -p "$MP"
tA=$SECONDS
python "$G/scripts/h6_mpnn_charge.py" --complexes "${MINE[@]}" --outdir "$MP" \
  --weights "$W" --nseq "$NSEQ" --keep "$KEEP" --temperature 0.15 \
  --seed $((500 + SH)) --bias-weight 1.6 --target-charge 5.5 --charge-window 1.5 --charge-weight 1.2 --charge-basis full --bulk-charge-bias 0.85 --charge-rounds 3
rc=$?
echo "STAGE_A rc=$rc elapsed=$((SECONDS-tA))s"
[ $rc -eq 0 ] || { echo "FATAL mpnn"; echo "=== exit: $rc ==="; exit $rc; }

# ---------- B: fold (GPU) ----------
conda activate boltz2
Y="$OUT/yaml"; mkdir -p "$Y"
python "$G/scripts/h3_make_yamls.py" --mpnn-json "$MP/h6_mpnn.json" \
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
if [ "$LEFT" -lt 300 ]; then
  echo "WARN: under 5 min left, skipping the fold stage"
else
  tB=$SECONDS
  timeout "${LEFT}s" boltz predict "$Y" --out_dir "$OUT" --cache "$SHMC" --model boltz2 \
    --write_full_pae --override --accelerator gpu --devices 1 \
    --num_workers 2 --preprocessing-threads 2 --diffusion_samples 1 --seed 1
  echo "STAGE_B rc=$? elapsed=$((SECONDS-tB))s (124 = deadline guard, partial is fine)"
fi

# ---------- C: measure ----------
P=$(find "$OUT" -maxdepth 3 -type d -name predictions | head -1)
NS=$(find "$P" -name '*_model_0.cif' 2>/dev/null | wc -l)
echo "STAGE_C structures=$NS of $NY"
if [ "$NS" -gt 0 ]; then
  conda activate rfd3
  A="$G/analysis"
  python "$G/scripts/h3_measure.py" --pred-dir "$P" --manifest "$OUT/manifest.json" \
    --af3-model "$G/fold_glp1r_glp_model_0.cif" \
    --out-csv "$A/p14_${RUN}_sh${SH}_deviation.csv" || true
  python "$G/scripts/h7_interactions.py" --pred-dir "$P" \
    --out-csv "$A/p14_${RUN}_sh${SH}_interactions.csv" || true
  # each backbone carries its own hotspot set; run the epitope check per group
  for SET in "p3:86,87,88,89" "p10:27,42,86" "p12:27,42,86" "p13:27,42,86,95"; do
    PFX="${SET%%:*}"; HOTS="${SET##*:}"
    if ls "$P" | grep -q "^${PFX}_"; then
      python "$G/scripts/h5_epitope.py" --pred-dir "$P" --manifest "$OUT/manifest.json" \
        --hotspots "$HOTS" --name-prefix "${PFX}_" \
        --dev-csv "$A/p14_${RUN}_sh${SH}_deviation.csv" \
        --out-csv "$A/p14_${RUN}_sh${SH}_epitope-${PFX}.csv" || true
    fi
  done
fi
echo "$OUT" >> "$G/logs/p14_${RUN}_outdirs.txt"
echo "TOTAL elapsed=$((SECONDS-T0))s"
echo "=== exit: 0 ==="
exit 0
