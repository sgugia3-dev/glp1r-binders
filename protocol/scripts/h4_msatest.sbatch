#!/usr/bin/env bash
#SBATCH --job-name=h4-msatest
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/h4_msatest_%j.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/h4_msatest_%j.out
#
# Does Boltz 2.2.1 accept the AlphaFold3 .a3m directly as `msa:`?
# One YAML, receptor -> cached a3m, designed binder -> msa: empty.
# Success is NOT the exit code: a bad MSA yields manifest {"records": []} at exit 0.
# So the record count and the structure file are both asserted.
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate boltz2

G="$C26_ROOT/glp1r2026"
A3M="$G/inputs/msa/glp1r_receptor_A.a3m"
MP="$(cat "$G/logs/h2_mpnn_outdir.txt")/h2_mpnn.json"
OUT="$C26_SCRATCH/runs/glp1r/h4_msatest_${SLURM_JOB_ID}"
Y="$OUT/yaml"
for i in 1 2 3 4 5; do mkdir -p "$Y" && break || { echo "mkdir retry $i"; sleep 20; }; done
[ -s "$A3M" ] || { echo "FATAL: no a3m"; echo "=== exit: 1 ==="; exit 1; }

python - "$MP" "$G/fold_glp1r_glp_job_request.json" "$A3M" "$Y" <<'PY'
import json, sys
mp, req, a3m, out = sys.argv[1:5]
rec = json.load(open(req))[0]["sequences"][0]["proteinChain"]["sequence"]
d = json.load(open(mp))["designs"][0]
for tag, msa in (("withmsa", a3m), ("nomsa", "empty")):
    with open(f"{out}/{tag}.yaml", "w") as fh:
        fh.write("version: 1\nsequences:\n")
        fh.write(f"  - protein:\n      id: A\n      sequence: {d['seq']}\n      msa: empty\n")
        fh.write(f"  - protein:\n      id: B\n      sequence: {rec}\n      msa: {msa}\n")
print("binder len", d["binder_len"], "receptor len", len(rec))
PY
echo "yamls: $(ls $Y/*.yaml | wc -l)"; grep -H "msa:" $Y/*.yaml | sed 's/:.*msa:/  msa:/'

SHMC="/dev/shm/boltzcache_${SLURM_JOB_ID}"
cleanup(){ rm -rf "$SHMC" 2>/dev/null || true; }
trap cleanup EXIT
mkdir -p "$SHMC"
cp "$C26_WEIGHTS/boltz/boltz2_conf.ckpt" "$SHMC/" || { echo "FATAL ckpt"; echo "=== exit: 1 ==="; exit 1; }
for x in mols mols.tar boltz2_aff.ckpt; do
  [ -e "$C26_WEIGHTS/boltz/$x" ] && ln -sfn "$C26_WEIGHTS/boltz/$x" "$SHMC/$x"
done
export NUMBA_CACHE_DIR="$OUT/numba"; mkdir -p "$NUMBA_CACHE_DIR"

t0=$SECONDS
boltz predict "$Y" --out_dir "$OUT" --cache "$SHMC" --model boltz2 \
  --write_full_pae --override --accelerator gpu --devices 1 \
  --num_workers 2 --preprocessing-threads 2 --diffusion_samples 1 --seed 1
rc=$?
echo "boltz rc=$rc wall $((SECONDS-t0))s"

echo "--- the checks that actually decide it ---"
MAN=$(find "$OUT" -name manifest.json | head -1)
echo "manifest: $MAN"
[ -n "$MAN" ] && python -c "
import json,sys; m=json.load(open('$MAN')); r=m.get('records',[])
print('manifest records:', len(r), '(expected 2)')
sys.exit(0 if len(r)==2 else 9)"
mrc=$?
P=$(find "$OUT" -maxdepth 3 -type d -name predictions | head -1)
echo "structures: $(find "$P" -name '*_model_0.cif' 2>/dev/null | wc -l) (expected 2)"
for t in withmsa nomsa; do
  c=$(find "$P" -path "*$t*" -name '*_model_0.cif' 2>/dev/null | head -1)
  if [ -n "$c" ]; then
    echo "  $t: OK $(basename $c)"
    python -c "
import json,glob
g=glob.glob('$P/$t/confidence_*_model_0.json')
if g:
    d=json.load(open(g[0]))
    print('    ipTM %.3f  pTM %.3f  plddt %.3f' % (d.get('iptm',-1), d.get('ptm',-1), d.get('complex_plddt',-1)))
"
  else
    echo "  $t: MISSING"
  fi
done
echo "failed-to-process: $(grep -c 'Failed to process' "$G/logs/h4_msatest_${SLURM_JOB_ID}.out" 2>/dev/null || echo 0)"
echo "=== exit: $((rc||mrc)) ==="
exit $((rc||mrc))
