#!/usr/bin/env python
"""Anchored arms: is the binder compact, is the motif held, and does the DESIGNED part
contact the intended patch? All three, per design.

The motif is excluded from the contact calculation: a contact made by the fixed peptide
residues is not a design result. Motif positions come from each design's own
diffused_index_map.
"""
import glob, gzip, json, math, os, statistics as st, sys

CUT = 5.0
d, hots = sys.argv[1], [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else []


def read(p):
    hdr, rows, inl = [], [], False
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                hdr.append(line.strip().split(".")[1]); inl = True; continue
            if inl and line.startswith(("ATOM", "HETATM")): rows.append(line.split())
            elif inl and line.startswith("#"): break
    i = {h: k for h, k in zip(hdr, range(len(hdr)))}
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


R, H, D, HELD, L = [], [], [], [], []
for cif in sorted(glob.glob(os.path.join(d, "*.cif.gz"))):
    j = cif[:-7] + ".json"
    if not os.path.exists(j): continue
    meta = json.load(open(j))
    m = meta.get("metrics", {})
    im = meta.get("diffused_index_map", {})
    motif = {int(v[1:]) for k, v in im.items() if k.startswith("B") and v.startswith("A")}
    held = m.get("max_ca_deviation", 99) <= 2.0 and m.get("n_chainbreaks", 9) <= 1
    ca, heavy = read(cif)
    if len(ca) < 2: continue
    tgt = max(ca, key=lambda c: len(ca[c]))
    bnd = [c for c in ca if c != tgt][0]
    pts = list(ca[bnd].values()); n = len(pts)
    if n < 10: continue
    c0 = [sum(q[k] for q in pts) / n for k in range(3)]
    rg = math.sqrt(sum(sum((q[k] - c0[k]) ** 2 for k in range(3)) for q in pts) / n)
    ratio = rg / (2.2 * n ** 0.38)
    off = min(heavy[tgt]) - 1
    des = [x for k, v in heavy[bnd].items() if k not in motif for x in v]
    ds = [min(math.dist(a, b) for a in des for b in heavy[tgt][h + off])
          for h in hots if (h + off) in heavy[tgt]] if des else []
    R.append(ratio); L.append(n); HELD.append(held)
    H.append(sum(1 for x in ds if x <= CUT) if ds else -1)
    D.append(min(ds) if ds else 99.0)

if not R: print("n=0 NOTHING MEASURED"); sys.exit(1)
joint = sum(1 for r, h, k in zip(R, H, HELD) if r <= 1.3 and h >= 1 and k)
print(f"n={len(R)} blen_med={st.median(L):.0f} motif_held={sum(HELD)}/{len(HELD)} "
      f"rg_med={st.median(R):.3f} rg_min={min(R):.3f} "
      f"hits_med={st.median(H):.1f} hits_max={max(H)} mind_med={st.median(D):.2f} "
      f"any_contact={sum(1 for x in H if x>=1)} JOINT={joint}/{len(R)}")
