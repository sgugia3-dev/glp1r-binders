#!/usr/bin/env python
"""Re-analysis: compactness of the binder from a chosen residue onward.

Rationale: if the extension leaves the anchor as a reaching stalk and only then folds,
whole-chain Rg is dominated by the stalk and understates the distal domain. So Rg is
recomputed over binder residues [FROM..L] only, with the globular expectation taken from
the SEGMENT length (2.2 * n_seg^0.38), not the full chain.

Contacts are unchanged: still measured from the designed region (motif excluded).
Usage: r_seg_measure.py <dir> <hotspot,list> <FROM residue, 1-based>
"""
import glob, gzip, json, math, os, statistics as st, sys

CUT = 5.0
d = sys.argv[1]
hots = [int(x) for x in sys.argv[2].split(",")]
FROM = int(sys.argv[3]) if len(sys.argv) > 3 else 30


def read(p):
    hdr, rows, inl = [], [], False
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                hdr.append(line.strip().split(".")[1]); inl = True; continue
            if inl and line.startswith(("ATOM", "HETATM")): rows.append(line.split())
            elif inl and line.startswith("#"): break
    i = {h: n for n, h in enumerate(hdr)}
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


def rg_of(pts):
    n = len(pts)
    c = [sum(q[k] for q in pts) / n for k in range(3)]
    return math.sqrt(sum(sum((q[k] - c[k]) ** 2 for k in range(3)) for q in pts) / n)


FULL, SEG, H, D, HELD, SL = [], [], [], [], [], []
for cif in sorted(glob.glob(os.path.join(d, "*.cif.gz"))):
    j = cif[:-7] + ".json"
    if not os.path.exists(j): continue
    meta = json.load(open(j)); m = meta.get("metrics", {})
    im = meta.get("diffused_index_map", {})
    motif = {int(v[1:]) for k, v in im.items() if k.startswith("B") and v.startswith("A")}
    held = m.get("max_ca_deviation", 99) <= 2.0 and m.get("n_chainbreaks", 9) <= 1
    ca, heavy = read(cif)
    if len(ca) < 2: continue
    tgt = max(ca, key=lambda c: len(ca[c]))
    bnd = [c for c in ca if c != tgt][0]
    allpos = sorted(ca[bnd])
    if len(allpos) < 10: continue
    full = [ca[bnd][k] for k in allpos]
    segpos = [k for k in allpos if k >= FROM]
    if len(segpos) < 5: continue
    seg = [ca[bnd][k] for k in segpos]
    FULL.append(rg_of(full) / (2.2 * len(full) ** 0.38))
    SEG.append(rg_of(seg) / (2.2 * len(seg) ** 0.38))
    SL.append(len(seg)); HELD.append(held)
    off = min(heavy[tgt]) - 1
    des = [x for k, v in heavy[bnd].items() if k not in motif for x in v]
    ds = [min(math.dist(a, b) for a in des for b in heavy[tgt][h + off])
          for h in hots if (h + off) in heavy[tgt]] if des else []
    H.append(sum(1 for x in ds if x <= CUT) if ds else -1)
    D.append(min(ds) if ds else 99.0)

if not SEG: print("n=0"); sys.exit(1)
joint = sum(1 for r, h, k in zip(SEG, H, HELD) if r <= 1.3 and h >= 1 and k)
print(f"n={len(SEG)} seg_len_med={st.median(SL):.0f} "
      f"rg_FULL_med={st.median(FULL):.3f} "
      f"rg_SEG{FROM}_med={st.median(SEG):.3f} rg_SEG_min={min(SEG):.3f} "
      f"n_seg_compact={sum(1 for r in SEG if r <= 1.3)}/{len(SEG)} "
      f"hits_med={st.median(H):.1f} JOINT_SEG={joint}/{len(SEG)}")
