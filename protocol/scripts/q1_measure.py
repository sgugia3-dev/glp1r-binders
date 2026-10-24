#!/usr/bin/env python
"""Measure the rfd3 unanchored arm the same way as the classic arms.

rfd3 writes gzipped mmCIF + a JSON carrying radius_of_gyration. Output convention is
binder = chain A numbered from 1, target = chain B renumbered from 1, so an original
target residue r is simply r in chain B here (no offset) -- asserted by checking the
target chain length matches the input.
"""
import glob, gzip, json, math, os, statistics as st, sys

HOT = [24, 25, 27, 42]
CUT = 5.0
d = sys.argv[1]


def read_cifgz(p):
    hdr, rows, inl = [], [], False
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                hdr.append(line.strip().split(".")[1]); inl = True; continue
            if inl and line.startswith(("ATOM", "HETATM")): rows.append(line.split())
            elif inl and line.startswith("#"): break
    i = {h: k for k, h in enumerate(hdr)}
    ca, heavy = {}, {}
    for r in rows:
        if len(r) <= max(i.values()): continue
        nm = r[i["label_atom_id"]]
        if nm.startswith("H"): continue
        try: n = int(r[i["auth_seq_id"]])
        except ValueError: continue
        c = r[i["auth_asym_id"]]
        xyz = (float(r[i["Cartn_x"]]), float(r[i["Cartn_y"]]), float(r[i["Cartn_z"]]))
        heavy.setdefault(c, {}).setdefault(n, []).append(xyz)
        if nm == "CA": ca.setdefault(c, {})[n] = xyz
    return ca, heavy


R, H, D, L = [], [], [], []
for cif in sorted(glob.glob(os.path.join(d, "*.cif.gz"))):
    ca, heavy = read_cifgz(cif)
    if len(ca) < 2: continue
    tgt = max(ca, key=lambda c: len(ca[c]))
    bnd = [c for c in ca if c != tgt][0]
    pts = list(ca[bnd].values()); n = len(pts)
    if n < 10: continue
    c0 = [sum(q[k] for q in pts) / n for k in range(3)]
    rg = math.sqrt(sum(sum((q[k] - c0[k]) ** 2 for k in range(3)) for q in pts) / n)
    exp = 2.2 * n ** 0.38
    off = min(heavy[tgt]) - 1
    bh = [x for v in heavy[bnd].values() for x in v]
    ds = [min(math.dist(a, b) for a in bh for b in heavy[tgt][h + off])
          for h in HOT if (h + off) in heavy[tgt]]
    R.append(rg / exp); L.append(n)
    H.append(sum(1 for x in ds if x <= CUT) if ds else -1)
    D.append(min(ds) if ds else -1)

if not R:
    print("SUMMARY n=0 NOTHING MEASURED"); sys.exit(1)
print(f"n={len(R)}  binder_len med={st.median(L):.0f}")
print(f"rg_ratio_med={st.median(R):.3f}  rg_ratio_min={min(R):.3f}  rg_ratio_max={max(R):.3f}")
print(f"hits_med={st.median(H):.1f}  hits_max={max(H)}  min_dist_med={st.median(D):.2f}")
print(f"n_with_any_contact={sum(1 for x in H if x >= 1)}/{len(H)}")
print(f"n_compact_and_contacting={sum(1 for r, h in zip(R, H) if r <= 1.3 and h >= 1)}/{len(R)}")
print(f"n_compact(<=1.3)={sum(1 for r in R if r <= 1.3)}/{len(R)}")
