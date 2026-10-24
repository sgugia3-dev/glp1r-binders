import csv, glob, os, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PPBuilder, MMCIFParser
G = "/lustre/home/gsgueglia/challenge2026/glp1r2026"
pp_, cp, ppb = PDBParser(QUIET=True), MMCIFParser(QUIET=True), PPBuilder()
for d in ("p2_best", "p3_best"):
    D = os.path.join(G, d)
    print("=" * 56); print(d)
    pdbs = sorted(glob.glob(D + "/binder_only/*.pdb"))
    cifs = sorted(glob.glob(D + "/structures/*.cif"))
    rows = list(csv.DictReader(open(D + "/metrics.csv")))
    print(f"  metrics rows {len(rows)}  pdbs {len(pdbs)}  cifs {len(cifs)}")
    bad = []
    lens = {}
    for f in pdbs:
        try:
            s = pp_.get_structure("x", f); ch = list(s[0])
            res = [r for r in ch[0] if r.id[0] == " "]
            lens[len(res)] = lens.get(len(res), 0) + 1
            if len(ch) != 1: bad.append((f, "chains=%d" % len(ch)))
            bf = [at.get_bfactor() for at in ch[0].get_atoms()]
            if not (0 < max(bf) <= 100): bad.append((f, "bfactor %.1f-%.1f" % (min(bf), max(bf))))
        except Exception as e: bad.append((f, repr(e)[:60]))
    print(f"  binder residue counts {lens}   problems {len(bad)}")
    for f, m in bad[:3]: print("    ", os.path.basename(f), m)
    cb = []
    for f in cifs:
        try:
            s = cp.get_structure("x", f)
            n = sorted(len([r for r in c if r.id[0] == " "]) for c in s[0])
            if n != [140, 409]: cb.append((os.path.basename(f), n))
        except Exception as e: cb.append((os.path.basename(f), repr(e)[:50]))
    print(f"  complexes not 140+409 in 2 chains: {len(cb)}")
    for m in cb[:3]: print("    ", m)
    # FASTA vs PDB vs metrics
    seqs, k = {}, None
    for line in open(D + "/candidates.fasta"):
        if line.startswith(">"): k = line[1:].split()[0]; seqs[k] = ""
        else: seqs[k] += line.strip()
    mism = 0
    for f in pdbs:
        tag = os.path.basename(f).replace("_binder.pdb", "")
        got = "".join(str(x.get_sequence()) for x in ppb.build_peptides(pp_.get_structure("x", f)))
        if seqs.get(tag, "") != got: mism += 1
    mm2 = sum(1 for r in rows if r["seq"] != seqs.get(r["tag"], ""))
    print(f"  FASTA vs PDB mismatches {mism}/{len(pdbs)}   metrics vs FASTA {mm2}/{len(rows)}")
    h = [int(r["hotspots"]) for r in rows]
    print(f"  hotspot counts in the set: {{{', '.join(f'{v}:{h.count(v)}' for v in sorted(set(h), reverse=True))}}}")
    anch = set(r["seq"][:10] for r in rows)
    print(f"  anchor(s): {sorted(anch)}")
    # README must not mention the other campaign
    rd = open(D + "/README.txt").read()
    other = "p3_1088996" if d == "p2_best" else "p2_1088883"
    wrong_hot = "A86" if d == "p2_best" else "A24"
    print(f"  README names its own campaign: {('p2_1088883' if d=='p2_best' else 'p3_1088996') in rd}")
    print(f"  README leaks the other campaign id: {other in rd}")
    print(f"  README leaks the other hotspots: {wrong_hot in rd}")
