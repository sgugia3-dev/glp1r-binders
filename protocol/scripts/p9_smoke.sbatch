#!/usr/bin/env bash
#SBATCH --job-name=p9-smoke
#SBATCH --partition=gpus
#SBATCH --qos=gpus
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p9_smoke_%A_%a.out
#SBATCH --error=/lustre/home/gsgueglia/challenge2026/glp1r2026/logs/p9_smoke_%A_%a.out
#
# Smoke test: can classic RFdiffusion run the ANCHORED architecture (140-residue
# binder + 409-residue receptor = 549 total) on a V100, and does a contact
# potential fire correctly on the 86-89 patch?
#
# Everything before now used rfd3, which has NO geometric potential. Potentials
# only exist in this checkout. Earlier contact-potential arms (s4/s5/p8) used an
# ECD-only target and an UNANCHORED binder, so neither the size nor the anchor
# has ever been tested with a potential.
#
#   0  hotspots only                      <- baseline, isolates the stack change
#   1  hotspots + interface_ncontacts     <- global: contacts to ALL 409 residues
#   2  hotspots + patch_ncontacts 86-89   <- restricted to the patch
set -o pipefail
ulimit -u 8192 2>/dev/null || true
source /lustre/home/gsgueglia/anaconda1/etc/profile.d/conda.sh
conda activate SE3nv
G=/lustre/home/gsgueglia/challenge2026/glp1r2026
R=/lustre/home/gsgueglia/diff_test/RFdiffusion
S=${SLURM_ARRAY_TASK_ID:-0}
NDES=${NDES:-1}
PDB="$G/inputs/rfd2/glp1r_full_pep10.pdb"
OUT="/lustre/scratch/gsgueglia/challenge2026/runs/glp1r/p9_smoke_${SLURM_ARRAY_JOB_ID:-m}/arm$S"
for i in 1 2 3 4 5; do mkdir -p "$OUT" && break || { echo "mkdir retry $i"; sleep 20; }; done
HOT='[A86,A87,A88,A89]'
case $S in
  0) POT=""                                                      ; LBL="hotspots_only" ;;
  1) POT="interface_ncontacts,weight:1"                          ; LBL="iface_ncontacts_w1" ;;
  2) POT="patch_ncontacts,weight:1,patch_lo:86,patch_hi:89"      ; LBL="patch_ncontacts_w1" ;;
esac
echo "arm $S ($LBL)  node $(hostname)"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader | head -1
python -c "import torch;print('torch',torch.__version__,'cuda',torch.cuda.is_available())" || exit 1
ARGS=(
  "inference.input_pdb=$PDB"
  "inference.output_prefix=$OUT/des"
  "inference.num_designs=$NDES"
  "contigmap.contigs=[B1-10/130-130/0 A1-409]"
  "ppi.hotspot_res=$HOT"
)
if [ -n "$POT" ]; then
  ARGS+=( "potentials.guiding_potentials=[\"type:$POT\"]"
          "potentials.guide_scale=10" "potentials.guide_decay=quadratic" )
fi
echo "ARGS: ${ARGS[*]}"
cd "$R" || exit 2
time python scripts/run_inference.py "${ARGS[@]}" 2>&1 | tail -60
rc=${PIPESTATUS[0]}
n=$(ls "$OUT"/des_*.pdb 2>/dev/null | wc -l)
echo "arm $S ($LBL) designs=$n rc=$rc"
if [ "$n" -gt 0 ]; then
  python - <<PY
import glob
f=sorted(glob.glob("$OUT/des_*.pdb"))[0]
ch={}
for l in open(f):
    if l.startswith("ATOM"): ch.setdefault(l[21],set()).add(int(l[22:26]))
print("  output chains:", {c:(len(v),min(v),max(v)) for c,v in sorted(ch.items())})
AA3={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
     "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
     "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
s=[AA3.get(l[17:20].strip(),"X") for l in open(f) if l.startswith("ATOM") and l[12:16].strip()=="CA"]
print("  first 20 res of the output:", "".join(s[:20]))
PY
fi
echo "=== exit: $rc ==="
exit $rc
