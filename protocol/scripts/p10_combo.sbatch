#!/usr/bin/env bash
#SBATCH --job-name=p10-combo
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=05:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p10_combo_%A_%a.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p10_combo_%A_%a.out
#
# SMOKE TEST: the P2/P3 protocol exactly, with hotspots A27 + A42 + A86 - one
# patch drawn from each face. A27:CD1 and A42:OG1 come from P2, A86:CD2 from P3,
# so the atoms are the ones those campaigns actually steered on.
#
# Geometry measured before launching (scripts/hot_geom.py):
#   A27-A42  6.1 A, A42-A86 11.9 A, A27-A86 16.5 A -> max span 16.5 A
#   references: P2 24/25/27/42 span 11.6 A, P3 86-89 span 8.5 A
#   anchor to A27 37.4 A, to A42 37.2 A, to A86 47.1 A
# So A27+A42 are effectively one site and A86 a second site ~12-16 A away. The
# open question is whether one binder bridges both, or trades one for the other.
#
# Everything else - contig B1-10,130-130,/0,A1-409, infer_ori_strategy=hotspots,
# is_non_loopy, step_scale, gamma_0, MPNN settings, Boltz settings - is identical
# to p3_all.sbatch, so the hit rate is directly comparable to P2's 64% at 4/4
# and P3's 23% at 4/4 (note: 3 hotspots here, not 4, so compare the FRACTION).
#
# All three stages in ONE job so the GPU allocation is never released between them:
#   A rfd3 diffusion (GPU) -> B motif screen -> C ProteinMPNN (CPU, GPU idles) ->
#   D Boltz refold (GPU)   -> E deviation + epitope measurement
# The GPU sits idle during stage C. That is the deliberate cost of not freeing the
# partition; on these numbers it is ~10 min of a ~3 h job.
#
# Each array task is an independent shard with its own rfd3 seed, so shards do not
# duplicate designs (`seed` is a real config field: rfd3/engine.py:80).
#
# DEADLINE GUARD: folding stops accepting new work at DEADLINE_S and the job still
# measures whatever finished, rather than being killed mid-write at the wall limit.
set -o pipefail
ulimit -u 8192 2>/dev/null || true
ulimit -c 0 2>/dev/null || true
T0=$SECONDS
DEADLINE_S=${DEADLINE_S:-15300}          # 4 h 15 m of a 5 h wall

source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
G="$C26_ROOT/glp1r2026"
SH=${SLURM_ARRAY_TASK_ID:-0}
RUN="${SLURM_ARRAY_JOB_ID:-manual}"
OUT="$C26_SCRATCH/runs/glp1r/p10_${RUN}/sh${SH}"
for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
NB=${NBATCH:-1}; BS=${BATCH:-4}; NSEQ=${NSEQ:-16}; PERBB=${PERBB:-3}
echo "shard $SH on $(hostname)  out $OUT  designs=$((NB*BS*5))  nseq=$NSEQ  per_bb=$PERBB"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | head -1

# ---------- A: diffusion ----------
conda activate rfd3
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
CKPT=$(ls "$C26_WEIGHTS/rfd3"/*.ckpt 2>/dev/null | head -1)
[ -n "$CKPT" ] || { echo "FATAL: no rfd3 ckpt"; echo "=== exit: 1 ==="; exit 1; }
DIFF="$OUT/diff"; mkdir -p "$DIFF"
tA=$SECONDS
rc=1
for att in 1 2 3; do
  rfd3 design out_dir="$DIFF" inputs="$G/inputs/rfd3/w3_r27_42_86.json" ckpt_path="$CKPT" \
    n_batches="$NB" diffusion_batch_size="$BS" seed=$((9000 + SH)) \
    inference_sampler.step_scale=3 inference_sampler.gamma_0=0.2 \
    prevalidate_inputs=True skip_existing=True
  rc=$?
  [ $rc -eq 0 ] && break
  echo "diffusion attempt $att rc=$rc; have $(find "$DIFF" -name '*.cif.gz' | wc -l)"; sleep 30
done
ND=$(find "$DIFF" -name '*.cif.gz' | wc -l)
echo "STAGE_A designs=$ND rc=$rc elapsed=$((SECONDS-tA))s"
[ "$ND" -gt 0 ] || { echo "FATAL: diffusion produced nothing"; echo "=== exit: 1 ==="; exit 1; }

# ---------- B+C: screen and design (CPU) ----------
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W="$C26_WEIGHTS/proteinmpnn/soluble_v_48_020.pt"
[ -f "$W" ] || W="$C26_ROOT/tools/mosaic/src/boltz_binder_design/proteinmpnn/weights/soluble_v_48_020.pt"
MP="$OUT/mpnn"; mkdir -p "$MP"
tC=$SECONDS
python "$G/scripts/h2_mpnn.py" --indir "$DIFF" --outdir "$MP" --weights "$W" \
  --nseq "$NSEQ" --temperature 0.1 --seed $((71 + SH)) \
  --acid-pos 3 --max-ca-dev 2.0 --max-breaks 1
mrc=$?
echo "STAGE_C rc=$mrc elapsed=$((SECONDS-tC))s"
[ $mrc -eq 0 ] || { echo "FATAL: mpnn stage failed"; echo "=== exit: $mrc ==="; exit $mrc; }

# ---------- D: fold (GPU, same allocation) ----------
conda activate boltz2
Y="$OUT/yaml"; mkdir -p "$Y"
python "$G/scripts/h3_make_yamls.py" --mpnn-json "$MP/h2_mpnn.json" \
  --af3-request "$G/fold_glp1r_glp_job_request.json" --out-dir "$Y" \
  --per-backbone "$PERBB" \
  --receptor-msa "$G/inputs/msa/glp1r_receptor_A.a3m" \
  --manifest "$OUT/manifest.json" || { echo "=== exit: 1 ==="; exit 1; }
NY=$(ls "$Y"/*.yaml | wc -l); echo "yamls=$NY"

SHMC="/dev/shm/boltzcache_${SLURM_JOB_ID}"
cleanup(){ rm -rf "$SHMC" 2>/dev/null || true; }
trap cleanup EXIT
mkdir -p "$SHMC"
cp "$C26_WEIGHTS/boltz/boltz2_conf.ckpt" "$SHMC/" || { echo "FATAL ckpt stage"; echo "=== exit: 1 ==="; exit 1; }
for x in mols mols.tar boltz2_aff.ckpt; do
  [ -e "$C26_WEIGHTS/boltz/$x" ] && ln -sfn "$C26_WEIGHTS/boltz/$x" "$SHMC/$x"
done
export NUMBA_CACHE_DIR="$OUT/numba"; mkdir -p "$NUMBA_CACHE_DIR"

LEFT=$((DEADLINE_S - (SECONDS - T0)))
echo "STAGE_D budget_left=${LEFT}s for $NY folds"
if [ "$LEFT" -lt 300 ]; then
  echo "WARN: under 5 min of budget left, skipping the fold stage"
else
  tD=$SECONDS
  timeout "${LEFT}s" boltz predict "$Y" --out_dir "$OUT" --cache "$SHMC" --model boltz2 \
    --write_full_pae --override --accelerator gpu --devices 1 \
    --num_workers 2 --preprocessing-threads 2 --diffusion_samples 1 --seed 1
  brc=$?
  echo "STAGE_D rc=$brc elapsed=$((SECONDS-tD))s (124 = hit the deadline guard, partial is fine)"
fi

# ---------- E: measure whatever folded ----------
P=$(find "$OUT" -maxdepth 3 -type d -name predictions | head -1)
NS=$(find "$P" -name '*_model_0.cif' 2>/dev/null | wc -l)
echo "STAGE_E structures=$NS of $NY"
if [ "$NS" -gt 0 ]; then
  conda activate rfd3
  python "$G/scripts/h3_measure.py" --pred-dir "$P" --manifest "$OUT/manifest.json" \
    --af3-model "$G/fold_glp1r_glp_model_0.cif" \
    --out-csv "$G/analysis/p10_${RUN}_sh${SH}_deviation.csv" || true
  python "$G/scripts/r_seg_measure.py" "$DIFF" 27,42,86 30 > "$G/analysis/p10_${RUN}_sh${SH}_segrg.txt" 2>&1 || true
  echo -n "SEGRG-combo sh${SH}  "; tail -1 "$G/analysis/p10_${RUN}_sh${SH}_segrg.txt"
  python "$G/scripts/h5_epitope.py" --pred-dir "$P" --manifest "$OUT/manifest.json" \
    --hotspots 27,42,86 \
    --dev-csv "$G/analysis/p10_${RUN}_sh${SH}_deviation.csv" \
    --out-csv "$G/analysis/p10_${RUN}_sh${SH}_epitope.csv" || true
fi
echo "$OUT" >> "$G/logs/p10_${RUN}_outdirs.txt"
echo "TOTAL elapsed=$((SECONDS-T0))s"
echo "=== exit: 0 ==="
exit 0
