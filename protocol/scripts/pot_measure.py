#!/usr/bin/env python
"""Shared measurement for every compactness/hotspot experiment. Do not modify.

Usage:  python pot_measure.py <dir of des_*.pdb> [hotspot list, default 24,25,27,42]

Prints one TSV line per design plus a SUMMARY line. Definitions, fixed so every arm is
comparable:
  rg_ratio  = binder Rg / (2.2 * N^0.38), the globular expectation. ~1.0 folded, >1.3 extended.
  hits      = how many hotspot residues have a binder heavy atom within 5.0 A
  min_dist  = closest binder heavy atom to any hotspot residue
Classic RFdiffusion renumbers output continuously (binder chain A 1..L, target chain B
L+1..), so the target offset is DERIVED per design, never assumed.
"""
import glob, math, os, statistics as st, sys

CUT = 5.0


def read(p):
    ca, heavy = {}, {}
    for line in open(p):
        if not line.startswith("ATOM"): continue
        nm, ch = line[12:16].strip(), line[21]
        try: resi = int(line[22:26])
        except ValueError: continue
        if nm.startswith("H"): continue
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        heavy.setdefault(ch, {}).setdefault(resi, []).append(xyz)
        if nm == "CA": ca.setdefault(ch, {})[resi] = xyz
    return ca, heavy


def main():
    d = sys.argv[1]
    hot = [int(x) for x in (sys.argv[2] if len(sys.argv) > 2 else "24,25,27,42").split(",")]
    pdbs = sorted(glob.glob(os.path.join(d, "des_*.pdb")))
    print(f"# dir={d} hotspots={hot} cut={CUT}")
    print("design\tblen\trg\trg_exp\trg_ratio\thits\tmin_dist")
    R, H, D = [], [], []
    for p in pdbs:
        ca, heavy = read(p)
        if len(ca) < 2: continue
        tgt = max(ca, key=lambda c: len(ca[c]))
        bnd = [c for c in ca if c != tgt][0]
        pts = list(ca[bnd].values()); L = len(pts)
        if L < 10: continue
        c = [sum(q[k] for q in pts) / L for k in range(3)]
        rg = math.sqrt(sum(sum((q[k] - c[k]) ** 2 for k in range(3)) for q in pts) / L)
        exp = 2.2 * L ** 0.38
        off = min(heavy[tgt]) - 1
        bh = [x for v in heavy[bnd].values() for x in v]
        ds = [min(math.dist(a, b) for a in bh for b in heavy[tgt][h + off])
              for h in hot if (h + off) in heavy[tgt]]
        hits = sum(1 for x in ds if x <= CUT) if ds else -1
        md = min(ds) if ds else -1
        R.append(rg / exp); H.append(hits); D.append(md)
        print(f"{os.path.basename(p)}\t{L}\t{rg:.2f}\t{exp:.2f}\t{rg/exp:.3f}\t{hits}\t{md:.2f}")
    if not R:
        print("SUMMARY\tn=0\tNOTHING MEASURED"); return 1
    print(f"SUMMARY\tn={len(R)}\trg_ratio_med={st.median(R):.3f}\t"
          f"rg_ratio_min={min(R):.3f}\thits_med={st.median(H):.1f}\t"
          f"hits_max={max(H)}\tmin_dist_med={st.median(D):.2f}\t"
          f"n_with_any_contact={sum(1 for x in H if x >= 1)}\t"
          f"n_compact_and_contacting={sum(1 for r, h in zip(R, H) if r <= 1.3 and h >= 1)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
