#!/usr/bin/env python
"""Rank the P2 campaign on RMSD, ipSAE and pLDDT.

pLDDT choice: the complex is 409 receptor residues to ~140 binder residues, so
`complex_plddt` is three-quarters a constant and barely moves between designs.
This ranks on the BINDER's own mean pLDDT, read per-residue from the plddt npz,
and reports the motif and designed-region means separately.

Gate first, then rank. The epitope gate leads because ipSAE cannot see it: over
1800 folds earlier in this campaign, designs touching 4/4 hotspots had a median
ipSAE of 0.885 and designs touching NONE had 0.872.

Gate: motif RMSD < 2.0 A, >= 3 hotspots, receptor RMSD < 5.0 A
Rank: mean of three min-max normalised terms (motif RMSD inverted, ipSAE, binder pLDDT)
"""
import csv, glob, json, os, re, statistics as st, sys
import numpy as np

analysis, root, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]


def cif_chain_residues(path):
    """Ordered list of (chain, resnum) over CA atoms, in file order."""
    names, rows, inloop = [], [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                names.append(line.strip().split(".")[1]); inloop = True; continue
            if inloop and (line.startswith("ATOM") or line.startswith("HETATM")):
                t = line.split()
                if len(t) != len(names): continue
                d = dict(zip(names, t))
                if d.get("label_atom_id") != "CA": continue
                ch = d.get("auth_asym_id") or d.get("label_asym_id")
                rows.append((ch, d.get("auth_seq_id") or d.get("label_seq_id")))
            elif inloop and rows and line.startswith("#"):
                break
    return rows


def parse_ipsae(p):
    rows, hdr = {}, None
    for line in open(p):
        t = line.split()
        if not t: continue
        if t[0] == "Chn1": hdr = t; continue
        if hdr and len(t) >= len(hdr) - 1:
            d = dict(zip(hdr, t))
            if d.get("Type"): rows[d["Type"]] = d
    return rows


def load(pat):
    out = {}
    for f in sorted(glob.glob(os.path.join(analysis, pat))):
        m = re.search(r"_sh(\d+)_", f)
        sh = m.group(1) if m else "?"
        for r in csv.DictReader(open(f)): out[(sh, r["name"])] = r
    return out


dev, epi = load("p2_1088883_sh*_deviation.csv"), load("p2_1088883_sh*_epitope.csv")
print(f"deviation rows {len(dev)}  epitope rows {len(epi)}")

rows, miss = [], {"ipsae": 0, "conf": 0, "plddt": 0, "join": 0}
for d in sorted(glob.glob(os.path.join(root, "sh*", "**", "predictions", "*"), recursive=True)):
    if not os.path.isdir(d): continue
    name = os.path.basename(d)
    m = re.search(r"/sh(\d+)/", d)
    if not m: continue
    k = (m.group(1), name)
    if k not in dev or k not in epi: miss["join"] += 1; continue
    ipf = glob.glob(os.path.join(d, "*_15_15.txt"))
    cff = glob.glob(os.path.join(d, "confidence_*.json"))
    plf = glob.glob(os.path.join(d, "plddt_*.npz"))
    cif = [x for x in glob.glob(os.path.join(d, "*_model_0.cif"))]
    if not ipf: miss["ipsae"] += 1; continue
    if not cff: miss["conf"] += 1; continue
    tb = parse_ipsae(ipf[0]); mx = tb.get("max")
    if not mx: miss["ipsae"] += 1; continue

    def g(dd, key):
        try: return float(dd.get(key, "nan"))
        except ValueError: return float("nan")

    conf = json.load(open(cff[0]))
    pl_b = pl_m = pl_d = float("nan")
    blen = int(dev[k]["binder_len"])
    if plf and cif:
        try:
            arr = np.load(plf[0])["plddt"].ravel()
            res = cif_chain_residues(cif[0])
            if len(res) == len(arr):
                first = res[0][0]
                idx = [i for i, (ch, _) in enumerate(res) if ch == first]
                if len(idx) == blen:
                    b = arr[idx]
                    pl_b, pl_m, pl_d = float(b.mean()), float(b[:10].mean()), float(b[10:].mean())
            if pl_b != pl_b: miss["plddt"] += 1
        except Exception:
            miss["plddt"] += 1
    else:
        miss["plddt"] += 1

    rows.append({
        "shard": k[0], "name": name, "binder_len": blen,
        "motif_rmsd": float(dev[k]["motif_rmsd_vs_af3"]),
        "receptor_rmsd": float(dev[k]["receptor_rmsd"]),
        "hotspots": int(epi[k]["n_hotspots_hit"]),
        "min_dist": round(float(epi[k]["min_hotspot_dist"] or 99), 2),
        "ipsae": g(mx, "ipSAE"), "ipsae_cons": min((g(v, "ipSAE") for v in tb.values()),
                                                   default=float("nan")),
        "pdockq": g(mx, "pDockQ"),
        "plddt_binder": round(pl_b, 4), "plddt_motif": round(pl_m, 4),
        "plddt_design": round(pl_d, 4),
        "plddt_complex": round(float(conf.get("complex_plddt", -1)), 4),
        "iplddt": round(float(conf.get("complex_iplddt", -1)), 4),
        "iptm": round(float(conf.get("iptm", -1)), 4),
    })
print(f"joined {len(rows)}   missing {miss}")
if not rows: sys.exit(1)

ok = [r for r in rows if r["plddt_binder"] == r["plddt_binder"]]
print(f"binder pLDDT resolved for {len(ok)}/{len(rows)}")
gate = [r for r in ok if r["motif_rmsd"] < 2.0 and r["hotspots"] >= 3
        and r["receptor_rmsd"] < 5.0]
print(f"gated {len(gate)}/{len(ok)} ({100*len(gate)/len(ok):.0f}%)")
if not gate: sys.exit(1)


def norm(v, invert=False):
    lo, hi = min(v), max(v)
    if hi == lo: return [1.0] * len(v)
    o = [(x - lo) / (hi - lo) for x in v]
    return [1 - x for x in o] if invert else o


a = norm([r["motif_rmsd"] for r in gate], invert=True)
b = norm([r["ipsae"] for r in gate])
c = norm([r["plddt_binder"] for r in gate])
for r, x, y, z in zip(gate, a, b, c):
    r["n_rmsd"], r["n_ipsae"], r["n_plddt"] = round(x, 3), round(y, 3), round(z, 3)
    r["score"] = round((x + y + z) / 3, 4)
gate.sort(key=lambda r: -r["score"])

cols = ["score", "n_rmsd", "n_ipsae", "n_plddt", "shard", "name", "binder_len",
        "motif_rmsd", "receptor_rmsd", "hotspots", "min_dist", "ipsae", "ipsae_cons",
        "plddt_binder", "plddt_motif", "plddt_design", "plddt_complex", "iplddt",
        "iptm", "pdockq"]
with open(out_csv, "w") as fh:
    fh.write(",".join(cols) + "\n")
    for r in gate: fh.write(",".join(str(r[c]) for c in cols) + "\n")

print("\ngated-set distributions")
for key in ("motif_rmsd", "ipsae", "ipsae_cons", "plddt_binder", "plddt_design",
            "plddt_complex", "iplddt", "iptm", "pdockq"):
    v = [r[key] for r in gate if r[key] == r[key] and r[key] >= 0]
    if v: print(f"  {key:14s} min {min(v):6.3f}  med {st.median(v):6.3f}  max {max(v):6.3f}")

print("\ndoes each metric separate on-target from off-target?")
for key in ("ipsae", "plddt_binder", "iptm", "plddt_complex"):
    for h in (4, 3, 2, 1, 0):
        v = [r[key] for r in rows if r["hotspots"] == h and r[key] == r[key] and r[key] >= 0]
        if v: print(f"  {key:14s} {h} hits  n={len(v):4d}  med {st.median(v):.3f}")
    print()

print(f"{'score':>6s} {'design':34s} {'motif':>6s} {'hit':>4s} {'ipSAE':>6s} "
      f"{'pLDDTb':>7s} {'pLDDTd':>7s}")
for r in gate[:20]:
    print(f"{r['score']:6.3f} {r['name'].replace('prod_r3_prod_r3_','')[:34]:34s} "
          f"{r['motif_rmsd']:6.2f} {r['hotspots']:4d} {r['ipsae']:6.3f} "
          f"{r['plddt_binder']:7.3f} {r['plddt_design']:7.3f}")
print(f"\nwrote {out_csv}")
