#!/usr/bin/env python
"""Build the anchored input PDB for classic RFdiffusion.

Needs the FULL receptor, not the ECD-only file already in inputs/rfd2: the
10-residue anchor is GLP-1 7-16, which binds the TMD core, so against an
ECD-only target the anchor would sit in empty space.

Writes chain A = receptor 1-409 and chain B = peptide 1-10, both renumbered
from 1, so the contig  [B1-10/130-130/0 A1-409]  maps exactly.
"""
import sys

src, out = sys.argv[1], sys.argv[2]
KEEP_PEP = int(sys.argv[3]) if len(sys.argv) > 3 else 10


def read_cif(path):
    names, rows = [], []
    for line in open(path):
        if line.startswith("_atom_site."):
            names.append(line.strip().split(".")[1]); continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")): continue
        t = line.split()
        if len(t) == len(names): rows.append(dict(zip(names, t)))
    return rows


rows = read_cif(src)
ch = lambda r: r.get("auth_asym_id") or r.get("label_asym_id")
rn = lambda r: r.get("auth_seq_id") or r.get("label_seq_id")

chains = {}
for r in rows:
    chains.setdefault(ch(r), []).append(r)
sizes = {c: len({rn(r) for r in v}) for c, v in chains.items()}
print("chains in the source:", sizes)

rec = max(sizes, key=lambda c: sizes[c])
pep = min(sizes, key=lambda c: sizes[c])
if rec == pep:
    sys.exit("could not tell receptor from peptide")
print(f"receptor chain {rec} ({sizes[rec]} res) -> A ; peptide chain {pep} "
      f"({sizes[pep]} res) -> B, keeping 1-{KEEP_PEP}")

n = 0
lines = []
for src_ch, dst_ch, limit in ((rec, "A", None), (pep, "B", KEEP_PEP)):
    order, seen = [], set()
    for r in chains[src_ch]:
        k = rn(r)
        if k not in seen: seen.add(k); order.append(k)
    remap = {k: i + 1 for i, k in enumerate(order)}
    for r in chains[src_ch]:
        new = remap[rn(r)]
        if limit and new > limit: continue
        el = (r.get("type_symbol") or "").strip()
        if el == "H": continue
        an = (r.get("label_atom_id") or "").strip().strip('"')
        nm = an.ljust(4) if len(an) >= 4 or len(el) > 1 else " " + an.ljust(3)
        n += 1
        lines.append("ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
                     % (n, nm, (r.get("label_comp_id") or "UNK")[:3], dst_ch, new,
                        float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]),
                        1.0, 0.0, el[:2].rjust(2)))
    lines.append("TER\n")
with open(out, "w") as fh:
    fh.writelines(lines); fh.write("END\n")

# verify what we wrote
cnt = {}
for l in lines:
    if l.startswith("ATOM"):
        cnt.setdefault(l[21], set()).add(int(l[22:26]))
print("written:", {c: (len(v), min(v), max(v)) for c, v in cnt.items()}, f"atoms={n}")
AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
seq = {}
for l in lines:
    if l.startswith("ATOM") and l[12:16].strip() == "CA":
        seq.setdefault(l[21], []).append(AA3.get(l[17:20].strip(), "X"))
for c, v in seq.items():
    print(f"  chain {c}: {''.join(v)[:60]}{'...' if len(v) > 60 else ''}")
print(f"wrote {out}")
