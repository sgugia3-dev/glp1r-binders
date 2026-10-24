#!/usr/bin/env python
"""Per-atom SASA at the requested hotspot residues, plus reach from the motif.

Two things decide whether a hotspot spec is usable:
  1. EXPOSURE. rfd3 `select_hotspots` values are ATOM NAMES. Naming a buried atom
     silently wastes the run (recorded: 4 arms differing only in the atom spec gave
     0/32 designs within 5 A, the atoms having 0.0-0.5 A^2 exposure). SASA is computed
     on the RECEPTOR ALONE, since that is the surface a binder sees.
  2. REACH. The motif's coordinates are fixed by the input model, so the designed
     extension has to physically span from the motif to the hotspots.
"""
import sys
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

cif, targets = sys.argv[1], sys.argv[2]
want = set()
for part in targets.split(","):
    if "-" in part:
        a, b = part.split("-"); want.update(range(int(a), int(b) + 1))
    else:
        want.add(int(part))

f = pdbx.CIFFile.read(cif)
arr = pdbx.get_structure(f, model=1)
arr = arr[struc.filter_amino_acids(arr)]
rec = arr[arr.chain_id == "A"]
pep = arr[arr.chain_id == "B"]
print(f"receptor {struc.get_residue_count(rec)} res, peptide {struc.get_residue_count(pep)} res")

# SASA on the receptor alone
sasa = struc.sasa(rec, vdw_radii="Single", point_number=200)
print(f"\n=== per-atom SASA at residues {sorted(want)} (receptor alone) ===")
print(f"{'res':>6s} {'aa':>4s} {'atom':>6s} {'SASA':>7s}  {'verdict':<10s}")
best = {}
for i in range(rec.array_length()):
    r = int(rec.res_id[i])
    if r not in want: continue
    nm = rec.atom_name[i]
    s = float(sasa[i]) if not np.isnan(sasa[i]) else 0.0
    if nm in ("N", "C", "O"): continue          # backbone, not useful as a hotspot
    v = "EXPOSED" if s >= 10 else ("partial" if s >= 3 else "buried")
    print(f"{r:6d} {rec.res_name[i]:>4s} {nm:>6s} {s:7.2f}  {v:<10s}")
    key = r
    if s > best.get(key, (0, None))[0]: best[key] = (s, nm, rec.res_name[i])

print("\n=== most exposed side-chain atom per residue ===")
usable = {}
for r in sorted(best):
    s, nm, aa = best[r]
    ok = s >= 10
    print(f"  A{r} ({aa}): {nm} at {s:.1f} A^2  {'USABLE' if ok else 'TOO BURIED'}")
    if ok: usable[r] = nm

# reach: motif B1-10 to the hotspot residues
mca = pep[(pep.atom_name == "CA") & (pep.res_id <= 10)]
hca = rec[(rec.atom_name == "CA") & np.isin(rec.res_id, sorted(want))]
if mca.array_length() and hca.array_length():
    d = np.linalg.norm(mca.coord[:, None, :] - hca.coord[None, :, :], axis=-1)
    print(f"\n=== reach from the fixed motif (B1-10 CA) to the hotspots ===")
    print(f"  min {d.min():.1f} A   max {d.max():.1f} A")
    # the C-terminal motif residue is where the extension starts
    c = pep[(pep.atom_name == "CA") & (pep.res_id == 10)]
    if c.array_length():
        dc = np.linalg.norm(c.coord[0] - hca.coord, axis=-1)
        print(f"  from B10 (extension start): min {dc.min():.1f} A  max {dc.max():.1f} A")
        print(f"  a 40-90 residue helix spans ~{40*1.5:.0f}-{90*1.5:.0f} A end to end")
print("\nHOTSPOT_SPEC " + " ".join(f"A{r}:{nm}" for r, nm in sorted(usable.items())))
