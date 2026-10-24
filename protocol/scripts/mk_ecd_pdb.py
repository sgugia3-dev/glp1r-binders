#!/usr/bin/env python
"""Extract receptor chain A residues 1-117 (the ECD) from the AF3 model as PDB.

Classic RFdiffusion takes a PDB, not mmCIF. The ECD alone is used because the hotspot
patch (24/25/27/42) is entirely within it, it is the only part a soluble binder can reach,
and it is ~3.5x cheaper than the full 409-residue receptor.
"""
import sys
cif, out, lo, hi = sys.argv[1], sys.argv[2], 1, 117
hdr, rows, inl = [], [], False
for line in open(cif):
    if line.startswith("_atom_site."):
        hdr.append(line.strip().split(".")[1]); inl = True; continue
    if inl and line.startswith(("ATOM", "HETATM")): rows.append(line.split())
    elif inl and line.startswith("#"): break
i = {h: k for k, h in enumerate(hdr)}
n = 0
with open(out, "w") as fh:
    for r in rows:
        if r[i["auth_asym_id"]] != "A": continue
        try: resi = int(r[i["auth_seq_id"]])
        except ValueError: continue
        if not (lo <= resi <= hi): continue
        el = r[i["type_symbol"]] if "type_symbol" in i else r[i["label_atom_id"]][0]
        if el == "H": continue
        n += 1
        atom = r[i["label_atom_id"]]
        fh.write("ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n" % (
            n, (" " + atom) if len(atom) < 4 else atom, r[i["label_comp_id"]], resi,
            float(r[i["Cartn_x"]]), float(r[i["Cartn_y"]]), float(r[i["Cartn_z"]]), el))
    fh.write("TER\nEND\n")
print(f"wrote {out}: {n} atoms, chain A residues {lo}-{hi}")
