#!/usr/bin/env bash
#SBATCH --job-name=p4-ipsae
#SBATCH --partition=htc
#SBATCH --qos=htc
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p4_ipsae_%j.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p4_ipsae_%j.out
#
# ipSAE on every P1 fold. The project rule is explicit: ipSAE is the ranking metric,
# computed from the Boltz PAE + CIF at cutoffs 15/15, and ipTM must NOT be substituted
# (Boltz's global ipTM is recorded as unusable on a multi-copy target).
# ipSAE writes its output NEXT TO ITS INPUT, not the cwd, and silencing it hides both the
# errors and the fact it worked elsewhere -- so stderr is kept.
set -o pipefail
ulimit -u 8192 2>/dev/null || true
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
source /lustre/home/gsgueglia/challenge2026/scripts/env.sh
conda activate rfd3

G="$C26_ROOT/glp1r2026"
ROOT="$C26_SCRATCH/runs/glp1r/p1_1088733"
TOOL="$C26_ROOT/tools/ipsae/ipsae.py"
[ -f "$TOOL" ] || { echo "FATAL: no ipsae.py at $TOOL"; echo "=== exit: 1 ==="; exit 1; }

n=0; ok=0; miss=0; fail=0
for d in "$ROOT"/sh*/boltz_results_yaml/predictions/*/; do
  nm=$(basename "$d")
  cif="$d/${nm}_model_0.cif"; pae="$d/pae_${nm}_model_0.npz"
  if [ ! -f "$cif" ] || [ ! -f "$pae" ]; then miss=$((miss+1)); continue; fi
  n=$((n+1))
  # skip if already scored, so this is resumable
  if ls "$d"/*_15_15.txt >/dev/null 2>&1; then ok=$((ok+1)); continue; fi
  if python "$TOOL" "$pae" "$cif" 15 15 >/dev/null; then ok=$((ok+1)); else fail=$((fail+1)); fi
  [ $((n % 200)) -eq 0 ] && echo "  progress: $n scored, $fail failed, $((SECONDS))s"
done
echo "folds seen: $n   scored ok: $ok   failed: $fail   missing cif/pae: $miss"
echo "ipsae output files: $(find "$ROOT" -name '*_15_15.txt' | wc -l)"
echo "=== exit: 0 ==="
exit 0
