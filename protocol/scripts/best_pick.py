#!/usr/bin/env python
"""Rank the r3 arm per design on the joint objective and print the winner's details."""
import glob, gzip, json, math, os, sys
CUT = 5.0; HOT = [24, 25, 27, 42]; FROM = 30
d = sys.argv[1]

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

def rg(pts):
    n = len(pts); c = [sum(q[k] for q in pts)/n for k in range(3)]
    return math.sqrt(sum(sum((q[k]-c[k])**2 for k in range(3)) for q in pts)/n)

rows = []
for cif in sorted(glob.glob(os.path.join(d, "*.cif.gz"))):
    j = cif[:-7] + ".json"
    if not os.path.exists(j): continue
    meta = json.load(open(j)); m = meta["metrics"]
    im = meta["diffused_index_map"]
    motif = {int(v[1:]) for k, v in im.items() if k.startswith("B") and v.startswith("A")}
    ca, heavy = read(cif)
    tgt = max(ca, key=lambda c: len(ca[c])); bnd = [c for c in ca if c != tgt][0]
    pos = sorted(ca[bnd]); seg = [ca[bnd][k] for k in pos if k >= FROM]
    if len(seg) < 5: continue
    sr = rg(seg)/(2.2*len(seg)**0.38)
    off = min(heavy[tgt]) - 1
    des = [x for k, v in heavy[bnd].items() if k not in motif for x in v]
    ds = {h: min(math.dist(a, b) for a in des for b in heavy[tgt][h+off])
          for h in HOT if (h+off) in heavy[tgt]}
    rows.append({"f": os.path.basename(cif), "blen": len(pos), "seg": len(seg),
                 "seg_rg": round(sr, 3), "hits": sum(1 for v in ds.values() if v <= CUT),
                 "d": {k: round(v, 2) for k, v in ds.items()},
                 "dev": round(m["max_ca_deviation"], 2), "brk": m["n_chainbreaks"],
                 "helix": round(m.get("helix_fraction", 0), 2),
                 "ss": m.get("num_ss_elements", 0),
                 "clash": m.get("n_clashing.interresidue_clashes_w_sidechain", -1)})
rows.sort(key=lambda r: (-r["hits"], r["seg_rg"]))
print("%-44s %4s %4s %8s %5s %6s %4s %6s" % ("design","blen","seg","seg_rg","hits","dev","nSS","helix"))
for r in rows:
    print("%-44s %4d %4d %8.3f %5d %6.2f %4d %6.2f" % (
        r["f"].replace("r3_ecd_L130_r3_ecd_L130_","")[:44], r["blen"], r["seg"],
        r["seg_rg"], r["hits"], r["dev"], r["ss"], r["helix"]))
b = rows[0]
print("\nBEST: %s" % b["f"])
print("  binder %d res (motif 1-10 + %d designed), segment scored from res %d" % (b["blen"], b["blen"]-10, FROM))
print("  segment Rg ratio %.3f   motif deviation %.2f A   chainbreaks %d   sidechain clashes %d" % (b["seg_rg"], b["dev"], b["brk"], b["clash"]))
print("  helix fraction %.2f over %d SS elements" % (b["helix"], b["ss"]))
print("  hotspot distances from the DESIGNED region: %s" % b["d"])
