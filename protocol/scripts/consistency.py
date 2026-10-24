#!/usr/bin/env python
"""Three-seed consistency for the top 150 GLP-1R designs.

Seeds 1, 2 and 3 are independent Boltz predictions of byte-identical inputs, so
disagreement between them is prediction variance, not a design difference.

The question the EGFR submission treated as decisive is not the average across
folds but whether each criterion holds in EVERY fold: in that campaign only 2 of
9 restrained passes survived re-folding, and round 5 switched to requiring both
folds to pass rather than their mean. This applies the same test here for the
first time.
"""
import csv, glob, json, os, re, statistics as st, sys
import numpy as np

A, SEEDS, OUT = sys.argv[1], sys.argv[2].split(","), sys.argv[3]
SRC_HOTS = {"p3": 4, "p10": 3, "p12": 3, "p13": 4}


def parse_ipsae(p):
    hdr, rows = None, {}
    for line in open(p):
        t = line.split()
        if not t: continue
        if t[0] == "Chn1": hdr = t; continue
        if hdr and len(t) >= len(hdr) - 1:
            d = dict(zip(hdr, t))
            if d.get("Type"): rows[d["Type"]] = d
    return rows


def load_tables(pref):
    def grab(pat, key="name"):
        out = {}
        for f in sorted(glob.glob(os.path.join(A, pat))):
            for r in csv.DictReader(open(f)):
                out[r[key]] = r
        return out
    dev = grab(f"{pref}_sh*_deviation.csv")
    inter = grab(f"{pref}_sh*_interactions.csv")
    epi, src = {}, {}
    for s in SRC_HOTS:
        for f in sorted(glob.glob(os.path.join(A, f"{pref}_sh*_epitope-{s}.csv"))):
            for r in csv.DictReader(open(f)):
                epi[r["name"]] = r; src[r["name"]] = s
    return dev, inter, epi, src


def from_structures(root):
    """ipSAE and binder pLDDT straight from each prediction directory."""
    out = {}
    for d in glob.glob(os.path.join(root, "sh*", "boltz_results_yaml",
                                    "predictions", "*")):
        if not os.path.isdir(d): continue
        nm = os.path.basename(d)
        rec = {}
        ip = glob.glob(os.path.join(d, "*_15_15.txt"))
        if ip:
            mx = parse_ipsae(ip[0]).get("max")
            if mx:
                try: rec["ipsae"] = float(mx.get("ipSAE", "nan"))
                except ValueError: pass
        cf = glob.glob(os.path.join(d, "confidence_*.json"))
        if cf:
            c = json.load(open(cf[0]))
            rec["iptm"] = float(c.get("iptm", -1))
            rec["plddt_complex"] = float(c.get("complex_plddt", -1))
            rec["iplddt"] = float(c.get("complex_iplddt", -1))
        pl = glob.glob(os.path.join(d, "plddt_*.npz"))
        if pl:
            try:
                arr = np.load(pl[0])["plddt"].ravel()
                if len(arr) >= 140: rec["plddt_binder"] = float(arr[:140].mean())
            except Exception: pass
        out[nm] = rec
    return out


SEED_SRC = {}
for s in SEEDS:
    pref, root = s.split("=")
    dev, inter, epi, src = load_tables(pref)
    extra = from_structures(root)
    SEED_SRC[pref.split("_")[0]] = (dev, inter, epi, src, extra)
    print(f"{pref}: dev {len(dev)} epi {len(epi)} inter {len(inter)} struct {len(extra)}")

names = set.intersection(*[set(v[0]) & set(v[2]) for v in SEED_SRC.values()])
print(f"\ndesigns measured in all {len(SEED_SRC)} seeds: {len(names)}")
labels = list(SEED_SRC)

rows = []
for nm in sorted(names):
    rec = {"name": nm}
    ok = True
    for lab in labels:
        dev, inter, epi, src, extra = SEED_SRC[lab]
        d, e, i, x = dev[nm], epi[nm], inter.get(nm, {}), extra.get(nm, {})
        nt = SRC_HOTS[src[nm]]
        try:
            rec[f"{lab}_motif"] = float(d["motif_rmsd_vs_af3"])
            rec[f"{lab}_recept"] = float(d["receptor_rmsd"])
            rec[f"{lab}_hits"] = int(e["n_hotspots_hit"])
            rec[f"{lab}_full"] = int(e["n_hotspots_hit"]) == nt
            rec[f"{lab}_sb"] = int(i.get("unique_salt_bridges", 0) or 0)
            rec[f"{lab}_hc"] = int(i.get("hydrophobic_contacts", 0) or 0)
            rec[f"{lab}_ipsae"] = x.get("ipsae")
            rec[f"{lab}_plddt"] = x.get("plddt_binder")
            rec[f"{lab}_iptm"] = x.get("iptm")
        except (ValueError, KeyError):
            ok = False
    rec["source"] = SEED_SRC[labels[0]][3][nm]
    rec["n_total"] = SRC_HOTS[rec["source"]]
    if ok: rows.append(rec)
print(f"fully joined: {len(rows)}\n")

def vals(r, k):
    return [r[f"{l}_{k}"] for l in labels if r.get(f"{l}_{k}") is not None]

print("=" * 78)
print("SPREAD ACROSS SEEDS, per design (median over designs)")
print("=" * 78)
print(f"{'metric':16s} {'median range':>14s} {'median SD':>11s} {'worst range':>12s}")
for k, f in (("ipsae", "{:.3f}"), ("plddt", "{:.3f}"), ("iptm", "{:.3f}"),
             ("motif", "{:.2f}"), ("recept", "{:.2f}"), ("sb", "{:.1f}"), ("hc", "{:.0f}")):
    rg, sd = [], []
    for r in rows:
        v = vals(r, k)
        if len(v) == len(labels):
            rg.append(max(v) - min(v))
            sd.append(st.pstdev(v))
    if rg:
        print(f"{k:16s} {f.format(st.median(rg)):>14s} {f.format(st.median(sd)):>11s} "
              f"{f.format(max(rg)):>12s}")

print("\n" + "=" * 78)
print("DOES EACH CRITERION HOLD IN EVERY SEED?")
print("=" * 78)
def count(pred):
    per = [sum(1 for l in labels if pred(r, l)) for r in rows]
    return sum(1 for p in per if p == len(labels)), sum(1 for p in per if p > 0), per

tests = [
    ("full epitope contact", lambda r, l: r[f"{l}_full"]),
    (">= 3 hotspots", lambda r, l: r[f"{l}_hits"] >= 3),
    ("motif RMSD < 2 A", lambda r, l: r[f"{l}_motif"] < 2.0),
    ("receptor RMSD < 5 A", lambda r, l: r[f"{l}_recept"] < 5.0),
    ("binder pLDDT >= 0.80", lambda r, l: (r[f"{l}_plddt"] or 0) >= 0.80),
    ("ipSAE >= 0.50", lambda r, l: (r[f"{l}_ipsae"] or 0) >= 0.50),
    ("salt bridges >= 5", lambda r, l: r[f"{l}_sb"] >= 5),
]
print(f"{'criterion':24s} {'all 3 seeds':>12s} {'any seed':>10s} {'fragile':>9s}")
for lbl, pred in tests:
    a3, a1, per = count(pred)
    print(f"{lbl:24s} {a3:5d} ({100*a3/len(rows):3.0f}%) {a1:5d} ({100*a1/len(rows):3.0f}%) "
          f"{a1-a3:6d}")

print("\ncombined gate (motif<2 AND >=3 hits AND receptor<5), per seed and in all three:")
g = lambda r, l: (r[f"{l}_motif"] < 2.0 and r[f"{l}_hits"] >= 3 and r[f"{l}_recept"] < 5.0)
for l in labels:
    n = sum(1 for r in rows if g(r, l))
    print(f"  {l:8s} {n:4d}/{len(rows)} ({100*n/len(rows):3.0f}%)")
allg = sum(1 for r in rows if all(g(r, l) for l in labels))
anyg = sum(1 for r in rows if any(g(r, l) for l in labels))
print(f"  ALL 3    {allg:4d}/{len(rows)} ({100*allg/len(rows):3.0f}%)")
print(f"  any      {anyg:4d}/{len(rows)}   -> {anyg-allg} pass sometimes and fail sometimes")

print("\nEGFR-style gate (pLDDT>=0.80 AND ipSAE>=0.50) in all three:")
e = lambda r, l: (r[f"{l}_plddt"] or 0) >= 0.80 and (r[f"{l}_ipsae"] or 0) >= 0.50
ne = sum(1 for r in rows if all(e(r, l) for l in labels))
print(f"  {ne}/{len(rows)} ({100*ne/len(rows):.0f}%)")
both = sum(1 for r in rows if all(g(r, l) and e(r, l) for l in labels))
print(f"both gates in all three seeds: {both}/{len(rows)} ({100*both/len(rows):.0f}%)")

def spearman(x, y):
    def rank(v):
        o = sorted(range(len(v)), key=lambda i: v[i]); rk = [0.0]*len(v)
        i = 0
        while i < len(o):
            j = i
            while j+1 < len(o) and v[o[j+1]] == v[o[i]]: j += 1
            for k2 in range(i, j+1): rk[o[k2]] = (i+j)/2 + 1
            i = j+1
        return rk
    rx, ry = rank(x), rank(y)
    n = len(x); mx, my = sum(rx)/n, sum(ry)/n
    num = sum((a-mx)*(b-my) for a, b in zip(rx, ry))
    den = (sum((a-mx)**2 for a in rx) * sum((b-my)**2 for b in ry)) ** .5
    return num/den if den else 0.0

print("\n" + "=" * 78)
print("DO THE SEEDS AGREE ON THE RANKING? (Spearman between seed pairs)")
print("=" * 78)
print(f"{'metric':16s} " + "".join(f"{labels[i]}~{labels[j]}".rjust(14)
      for i in range(len(labels)) for j in range(i+1, len(labels))))
for k in ("ipsae", "plddt", "iptm", "motif", "sb", "hc"):
    line = f"{k:16s} "
    for i in range(len(labels)):
        for j in range(i+1, len(labels)):
            pair = [(r[f'{labels[i]}_{k}'], r[f'{labels[j]}_{k}']) for r in rows
                    if r.get(f'{labels[i]}_{k}') is not None
                    and r.get(f'{labels[j]}_{k}') is not None]
            line += f"{spearman([p[0] for p in pair],[p[1] for p in pair]):14.3f}" if pair else " " * 14
    print(line)

cols = ["name", "source", "n_total"] + [f"{l}_{k}" for l in labels for k in
       ("motif", "recept", "hits", "full", "sb", "hc", "ipsae", "plddt", "iptm")]
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(cols + ["n_seeds_gate", "n_seeds_full"])
    for r in rows:
        w.writerow([r.get(c, "") for c in cols]
                   + [sum(1 for l in labels if g(r, l)),
                      sum(1 for l in labels if r[f"{l}_full"])])
print(f"\nwrote {OUT}")
