#!/usr/bin/env python
"""Aggregate the combined-hotspot smoke test and test the bridging question.

A27 and A42 sit 6.1 A apart and are effectively one site; A86 is 12-16 A away.
So a high "3 of 3" rate could mean either a binder that genuinely bridges both
faces, or one parked on the easy A27/A42 pair with A86 scraping in. The decisive
number is the PER-HOTSPOT hit rate, and in particular A86 alone.
"""
import csv, glob, os, re, statistics as st, sys

A = sys.argv[1]


def load(run, pat):
    out = {}
    for f in sorted(glob.glob(os.path.join(A, f"{run}_sh*_{pat}.csv"))):
        sh = re.search(r"_sh(\d+)_", f).group(1)
        for r in csv.DictReader(open(f)):
            out[(sh, r["name"])] = r
    return out


dev, epi = load("p10_1089071", "deviation"), load("p10_1089071", "epitope")
keys = sorted(set(dev) & set(epi))
print(f"deviation {len(dev)}  epitope {len(epi)}  joined {len(keys)}")
hot_cols = [c for c in next(iter(epi.values())) if c.startswith("d_des_A")]
print(f"hotspot columns: {hot_cols}\n")

rows = []
for k in keys:
    d, e = dev[k], epi[k]
    try:
        r = {"shard": k[0], "name": k[1],
             "motif": float(d["motif_rmsd_vs_af3"]),
             "recept": float(d["receptor_rmsd"]),
             "hits": int(e["n_hotspots_hit"]),
             "mind": float(e["min_hotspot_dist"] or 99)}
        for c in hot_cols:
            r[c] = float(e[c]) if e[c] not in ("", None) else 99.0
        rows.append(r)
    except (ValueError, KeyError):
        pass
n = len(rows)
NH = len(hot_cols)

print("=" * 60)
print(f"P10 combined hotspots A27 + A42 + A86   ({n} refolds)")
u2 = [r for r in rows if r["motif"] < 2.0]
full = [r for r in rows if r["hits"] == NH]
gate = [r for r in rows if r["motif"] < 2.0 and r["hits"] >= 3 and r["recept"] < 5.0]
print(f"  motif under 2A        : {len(u2):4d} ({100*len(u2)/n:3.0f}%)")
print(f"  FULL {NH}/{NH} hotspots    : {len(full):4d} ({100*len(full)/n:3.0f}%)")
print(f"  PASS gate (>=3 hits)  : {len(gate):4d} ({100*len(gate)/n:3.0f}%)")
print(f"  complete misses (0)   : {sum(1 for r in rows if r['hits']==0):4d}")
print(f"  motif RMSD            : med {st.median([r['motif'] for r in rows]):.2f} A")
print(f"  hits histogram        : " + "  ".join(
    f"{h}:{sum(1 for r in rows if r['hits']==h)}" for h in range(NH, -1, -1)))

print(f"\nPER-HOTSPOT engagement (contact = within 5.0 A of the designed region)")
print(f"  {'hotspot':10s} {'hit rate':>10s} {'median dist':>12s} {'min':>7s}")
for c in sorted(hot_cols):
    ds = [r[c] for r in rows]
    hit = sum(1 for x in ds if x < 5.0)
    print(f"  {c.replace('d_des_',''):10s} {hit:4d}/{n} ({100*hit/n:3.0f}%) "
          f"{st.median(ds):12.2f} {min(ds):7.2f}")

# the bridging test: among designs holding the easy pair, how many also reach A86?
pair = [r for r in rows if r.get("d_des_A27", 99) < 5 and r.get("d_des_A42", 99) < 5]
both = [r for r in pair if r.get("d_des_A86", 99) < 5]
a86only = [r for r in rows if r.get("d_des_A86", 99) < 5
           and not (r.get("d_des_A27", 99) < 5 and r.get("d_des_A42", 99) < 5)]
print(f"\nBRIDGING TEST")
print(f"  hold the A27+A42 pair            : {len(pair):4d}/{n} ({100*len(pair)/n:3.0f}%)")
print(f"  of those, ALSO reach A86         : {len(both):4d}/{len(pair)} "
      f"({100*len(both)/max(len(pair),1):3.0f}%)")
print(f"  reach A86 without the pair       : {len(a86only):4d}")
print(f"  -> {'genuine bridging' if len(both) > 0.5*len(pair) else 'A86 is the weak point'}")

print(f"\nHEAD TO HEAD, full engagement of every hotspot in the set")
for run, lbl, nh in (("p2_1088883", "P2  A24/25/27/42", 4),
                     ("p3_1088996", "P3  A86-89", 4),
                     ("p10_1089071", "P10 A27/42/86", 3)):
    d2, e2 = load(run, "deviation"), load(run, "epitope")
    kk = sorted(set(d2) & set(e2))
    if not kk: continue
    hh = [int(e2[k]["n_hotspots_hit"]) for k in kk]
    mm = [float(d2[k]["motif_rmsd_vs_af3"]) for k in kk]
    fu = sum(1 for h in hh if h == nh)
    g3 = sum(1 for h in hh if h >= 3)
    z = sum(1 for h in hh if h == 0)
    print(f"  {lbl:18s} n={len(kk):4d}  full {nh}/{nh}: {100*fu/len(kk):3.0f}%   "
          f">=3: {100*g3/len(kk):3.0f}%   misses: {100*z/len(kk):3.0f}%   "
          f"motif<2A: {100*sum(1 for x in mm if x<2)/len(kk):3.0f}%")
print("\nNOTE: P10 has 3 hotspots, P2 and P3 have 4, so 'full' is an easier bar")
print("for P10. The per-hotspot and bridging numbers above are the fair comparison.")
