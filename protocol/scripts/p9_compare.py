#!/usr/bin/env python
"""Does binder_ROG give a compact binder AND keep hotspot contact? Plus the CFG sweep.

Compactness is reported as Rg / (2.2 * N^0.38), the globular expectation. ~1.0 is a
folded domain; >1.3 is extended.
"""
import glob, gzip, json, math, os, statistics as st, sys

HOT = [24, 25, 27, 42]


def pdb_chains(p):
    ca, heavy = {}, {}
    for line in open(p):
        if not line.startswith("ATOM"): continue
        nm = line[12:16].strip(); ch = line[21]; 
        try: resi = int(line[22:26])
        except ValueError: continue
        if nm.startswith("H"): continue
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        heavy.setdefault(ch, {}).setdefault(resi, []).append(xyz)
        if nm == "CA": ca.setdefault(ch, {})[resi] = xyz
    return ca, heavy


def rg(pts):
    c = [sum(p[k] for p in pts) / len(pts) for k in range(3)]
    return math.sqrt(sum(sum((p[k] - c[k]) ** 2 for k in range(3)) for p in pts) / len(pts))


def mind(a, b):
    return min(math.dist(p, q) for p in a for q in b)


print("=== classic RFdiffusion: binder_ROG potential ===")
print("%-16s %4s %7s %7s %7s %7s %7s" % ("arm", "n", "Rg", "exp", "ratio", "hits", "minD"))
for arm in sorted(glob.glob(os.path.join(sys.argv[1], "arm*"))):
    pdbs = sorted(glob.glob(os.path.join(arm, "des_*.pdb")))
    if not pdbs: continue
    rr, hh, dd = [], [], []
    for p in pdbs:
        ca, heavy = pdb_chains(p)
        if len(ca) < 2: continue
        tgt = max(ca, key=lambda c: len(ca[c]))
        bnd = [c for c in ca if c != tgt][0]
        pts = list(ca[bnd].values()); L = len(pts)
        if L < 10: continue
        rr.append(rg(pts) / (2.2 * L ** 0.38))
        bh = [x for v in heavy[bnd].values() for x in v]
        # classic RFdiffusion renumbers the output CONTINUOUSLY: binder is chain A 1..L and
        # the target follows as chain B L+1.., so an original target residue r appears at
        # r + offset where offset = min(target resnum) - 1. Derive it, never assume.
        off = min(heavy[tgt]) - 1
        ds = [mind(bh, heavy[tgt][h + off]) for h in HOT if (h + off) in heavy[tgt]]
        if ds:
            hh.append(sum(1 for d in ds if d <= 5.0)); dd.append(min(ds))
    if not rr: continue
    print("%-16s %4d %7.1f %7.1f %7.2f %7.1f %7.2f" % (
        os.path.basename(arm), len(rr),
        st.median([r * (2.2 * 70 ** 0.38) for r in rr]), 2.2 * 70 ** 0.38,
        st.median(rr), st.median(hh) if hh else -1, st.median(dd) if dd else -1))

print("\n=== rfd3: classifier-free guidance sweep ===")
print("%-10s %4s %7s %7s %6s %6s" % ("cfg_scale", "n", "Rg", "ratio", "helix", "motif_ok"))
for d in sorted(glob.glob(os.path.join(sys.argv[2], "cfg*")),
                key=lambda x: float(os.path.basename(x).replace("cfg", ""))):
    rs, hx, ok = [], [], 0
    for j in glob.glob(os.path.join(d, "*.json")):
        try: m = json.load(open(j))
        except Exception: continue
        mm = m.get("metrics", {}); sp = m.get("specification", {})
        sc = sp.get("extra", {}).get("sampled_contig", "")
        im = m.get("diffused_index_map", {})
        nb = len([v for k, v in im.items() if k.startswith("B") and v.startswith("A")])
        fl = sum(int(x[:-1]) for x in sc.split(",/0,")[0].split(",") if x.endswith("P"))
        L = fl + nb
        g = mm.get("radius_of_gyration")
        if not g or not L: continue
        rs.append((g, g / (2.2 * L ** 0.38)))
        hx.append(mm.get("helix_fraction", 0))
        if mm.get("max_ca_deviation", 99) <= 2.0 and mm.get("n_chainbreaks", 9) <= 1: ok += 1
    if not rs: continue
    print("%-10s %4d %7.1f %7.2f %6.2f %6s" % (
        os.path.basename(d).replace("cfg", ""), len(rs),
        st.median([a for a, _ in rs]), st.median([b for _, b in rs]),
        st.median(hx), f"{ok}/{len(rs)}"))
