import csv, glob, os, warnings
from collections import Counter
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PPBuilder, MMCIFParser
D = "/lustre/home/gsgueglia/challenge2026/glp1r2026/p14_best"
pp_, cp, ppb = PDBParser(QUIET=True), MMCIFParser(QUIET=True), PPBuilder()
pdbs = sorted(glob.glob(D + "/binder_only/*.pdb"))
cifs = sorted(glob.glob(D + "/structures/*.cif"))
rows = list(csv.DictReader(open(D + "/metrics.csv")))
print(f"metrics {len(rows)}  pdbs {len(pdbs)}  cifs {len(cifs)}")
lens, bad = {}, 0
for f in pdbs:
    s = pp_.get_structure("x", f); ch = list(s[0])
    r = [x for x in ch[0] if x.id[0] == " "]
    lens[len(r)] = lens.get(len(r), 0) + 1
    if len(ch) != 1: bad += 1
    bf = [a.get_bfactor() for a in ch[0].get_atoms()]
    if not (0 < max(bf) <= 100): bad += 1
print(f"binder residue counts {lens}  problems {bad}")
cb = [f for f in cifs
      if sorted(len([r for r in c if r.id[0] == " "]) for c in cp.get_structure("x", f)[0]) != [140, 409]]
print(f"complexes not 140+409: {len(cb)}")
seqs, k = {}, None
for line in open(D + "/candidates.fasta"):
    if line.startswith(">"): k = line[1:].split()[0]; seqs[k] = ""
    else: seqs[k] += line.strip()
mm = sum(1 for f in pdbs
         if seqs.get(os.path.basename(f).replace("_binder.pdb", ""), "")
         != "".join(str(x.get_sequence()) for x in ppb.build_peptides(pp_.get_structure("x", f))))
print(f"FASTA vs PDB mismatches {mm} | metrics vs FASTA "
      f"{sum(1 for r in rows if r['seq'] != seqs.get(r['tag'], ''))}")
print(f"CONSTRAINTS  pos1==H {sum(1 for r in rows if r['seq'][0]=='H')}/{len(rows)}   "
      f"pos3 in D/E {sum(1 for r in rows if r['seq'][2] in 'DE')}/{len(rows)}")
print(f"distinct anchors: {len({r['seq'][:10] for r in rows})}")
q = [Counter(r["seq"])["K"] + Counter(r["seq"])["R"]
     - Counter(r["seq"])["D"] - Counter(r["seq"])["E"] for r in rows]
import statistics as st
print(f"q(full 140): median {st.median(q):+.0f}  min {min(q):+.0f}  max {max(q):+.0f}")
print(f"hotspot counts: {dict(Counter(int(r['hotspots']) for r in rows))}")
print(f"source mix: {dict(Counter(r['name'].split('_')[0] for r in rows))}")
rd = open(D + "/README.txt").read()
print(f"README: own campaign {'p14_1092297' in rd} | says 600 refolds {'600 refolds' in rd}")
