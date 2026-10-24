#!/usr/bin/env python
"""Rank the P2 campaign on RMSD, ipSAE and pLDDT.

Gate first, then rank. The epitope gate comes first because ipSAE provably cannot see it:
across 1800 folds earlier in this campaign, designs contacting 4 of 4 hotspots scored a
median ipSAE of 0.885 and designs contacting NONE scored 0.872. A confidence metric that
cannot separate on-target from off-target cannot be the primary filter.

Gate   : motif RMSD < 2.0 A, >= 3 hotspots contacted, receptor RMSD < 5.0 A
Rank   : mean of three min-max normalised terms, all oriented so higher is better
         - motif RMSD (inverted)
         - ipSAE, the `max` row at cutoffs 15/15
         - complex pLDDT
Reports the conservative ipSAE direction too, since the `max` row is the optimistic one.
"""
import csv, glob, json, os, re, statistics as st, sys

analysis, root, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]


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
print(f"deviation rows {len(dev)}, epitope rows {len(epi)}")

rows, no_ipsae, no_conf = [], 0, 0
for d in glob.glob(os.path.join(root, "sh*", "**", "predictions", "*"), recursive=True):
    if not os.path.isdir(d): continue
    name = os.path.basename(d)
    m = re.search(r"/sh(\d+)/", d)
    if not m: continue
    sh = m.group(1)
    k = (sh, name)
    if k not in dev or k not in epi: continue
    ip = glob.glob(os.path.join(d, "*_15_15.txt"))
    cf = glob.glob(os.path.join(d, "confidence_*.json"))
    if not ip: no_ipsae += 1; continue
    if not cf: no_conf += 1; continue
    tb = parse_ipsae(ip[0])
    mx = tb.get("max")
    if not mx: no_ipsae += 1; continue
    def g(dd, key):
        try: return float(dd.get(key, "nan"))
        except ValueError: return float("nan")
    conf = json.load(open(cf[0]))
    rows.append({
        "shard": sh, "name": name,
        "binder_len": int(dev[k]["binder_len"]),
        "motif_rmsd": float(dev[k]["motif_rmsd_vs_af3"]),
        "receptor_rmsd": float(dev[k]["receptor_rmsd"]),
        "hotspots": int(epi[k]["n_hotspots_hit"]),
        "min_dist": float(epi[k]["min_hotspot_dist"] or 99),
        "ipsae": g(mx, "ipSAE"),
        "ipsae_cons": min((g(v, "ipSAE") for v in tb.values()), default=float("nan")),
        "pdockq": g(mx, "pDockQ"),
        "plddt": float(conf.get("complex_plddt", -1)),
        "iptm": float(conf.get("iptm", -1)),
    })
print(f"joined {len(rows)}  (missing ipsae {no_ipsae}, missing confidence {no_conf})")
if not rows: sys.exit(1)

gate = [r for r in rows if r["motif_rmsd"] < 2.0 and r["hotspots"] >= 3
        and r["receptor_rmsd"] < 5.0]
print(f"passing the gate: {len(gate)}/{len(rows)} ({100*len(gate)/len(rows):.0f}%)")
if not gate: sys.exit(1)


def norm(vals, invert=False):
    lo, hi = min(vals), max(vals)
    if hi == lo: return [1.0] * len(vals)
    out = [(v - lo) / (hi - lo) for v in vals]
    return [1 - x for x in out] if invert else out


nm = norm([r["motif_rmsd"] for r in gate], invert=True)
ni = norm([r["ipsae"] for r in gate])
np_ = norm([r["plddt"] for r in gate])
for r, a, b, c in zip(gate, nm, ni, np_):
    r["score"] = round((a + b + c) / 3, 4)
gate.sort(key=lambda r: -r["score"])

cols = ["score", "shard", "name", "binder_len", "motif_rmsd", "receptor_rmsd",
        "hotspots", "min_dist", "ipsae", "ipsae_cons", "plddt", "iptm", "pdockq"]
with open(out_csv, "w") as fh:
    fh.write(",".join(cols) + "\n")
    for r in gate: fh.write(",".join(str(r[c]) for c in cols) + "\n")

print("\ndistributions over the gated set:")
for key in ("motif_rmsd", "ipsae", "ipsae_cons", "plddt", "iptm", "pdockq"):
    v = [r[key] for r in gate if r[key] == r[key]]
    print(f"  {key:13s} min {min(v):6.3f}  med {st.median(v):6.3f}  max {max(v):6.3f}")
print(f"\n{'score':>6s} {'design':30s} {'sh':>3s} {'motif':>6s} {'hits':>5s} "
      f"{'ipSAE':>6s} {'cons':>6s} {'pLDDT':>6s}")
for r in gate[:15]:
    print(f"{r['score']:6.3f} {r['name'].replace('prod_r3_prod_r3_','')[:30]:30s} "
          f"{r['shard']:>3s} {r['motif_rmsd']:6.2f} {r['hotspots']:5d} "
          f"{r['ipsae']:6.3f} {r['ipsae_cons']:6.3f} {r['plddt']:6.3f}")
print(f"\nwrote {out_csv}")
