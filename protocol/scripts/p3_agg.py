#!/usr/bin/env python
"""Aggregate the P3 (Trp87 / 86-89 patch) campaign and compare it to P2 head to head.

Same gate as P2 so the two are comparable:
  motif RMSD vs AF3 < 2.0 A, >= 3 hotspots contacted, receptor RMSD < 5.0 A
Both campaigns have exactly 4 hotspots, so ">= 3" means the same thing in each.
P2 hotspots were A24/A25/A27/A42; P3 are A86/A87/A88/A89.
"""
import csv, glob, json, os, re, statistics as st, sys

analysis = sys.argv[1]


def load(run, pat):
    out = {}
    for f in sorted(glob.glob(os.path.join(analysis, f"{run}_sh*_{pat}.csv"))):
        sh = re.search(r"_sh(\d+)_", f).group(1)
        for r in csv.DictReader(open(f)):
            out[(sh, r["name"])] = r
    return out


def summarize(run, label):
    dev, epi = load(run, "deviation"), load(run, "epitope")
    keys = sorted(set(dev) & set(epi))
    rows = []
    for k in keys:
        try:
            rows.append({
                "shard": k[0], "name": k[1],
                "motif": float(dev[k]["motif_rmsd_vs_af3"]),
                "recept": float(dev[k]["receptor_rmsd"]),
                "hits": int(epi[k]["n_hotspots_hit"]),
                "mind": float(epi[k]["min_hotspot_dist"] or 99),
                "des": str(epi[k].get("designed_engages", "")).strip().lower() == "true",
            })
        except (ValueError, KeyError):
            pass
    n = len(rows)
    if not n:
        print(f"{label}: no rows"); return None
    u2 = [r for r in rows if r["motif"] < 2.0]
    h4 = [r for r in rows if r["hits"] == 4]
    h3 = [r for r in rows if r["hits"] >= 3]
    gate = [r for r in rows if r["motif"] < 2.0 and r["hits"] >= 3 and r["recept"] < 5.0]
    elite = [r for r in gate if r["hits"] == 4 and r["motif"] < 1.0]
    print(f"\n{label}   (shards {len(set(r['shard'] for r in rows))}, "
          f"dev {len(dev)} / epi {len(epi)} joined -> {n})")
    print(f"  refolds measured   : {n}")
    print(f"  motif under 2A     : {len(u2):5d} ({100*len(u2)/n:3.0f}%)")
    print(f"  4/4 hotspots       : {len(h4):5d} ({100*len(h4)/n:3.0f}%)")
    print(f"  >=3/4 hotspots     : {len(h3):5d} ({100*len(h3)/n:3.0f}%)")
    print(f"  designed engages   : {sum(r['des'] for r in rows):5d} "
          f"({100*sum(r['des'] for r in rows)/n:3.0f}%)")
    print(f"  PASS full gate     : {len(gate):5d}/{n} ({100*len(gate)/n:3.0f}%)")
    print(f"  ELITE              : {len(elite):5d}/{n} ({100*len(elite)/n:3.0f}%)")
    print(f"  motif RMSD         : med {st.median([r['motif'] for r in rows]):.2f} A  "
          f"min {min(r['motif'] for r in rows):.2f}")
    print(f"  min hotspot dist   : med {st.median([r['mind'] for r in rows]):.2f} A  "
          f"min {min(r['mind'] for r in rows):.2f}")
    print(f"  hits histogram     : " + "  ".join(
        f"{h}:{sum(1 for r in rows if r['hits']==h)}" for h in (4, 3, 2, 1, 0)))
    return {"label": label, "n": n, "gate": len(gate), "elite": len(elite),
            "h4": len(h4), "u2": len(u2), "rows": rows}


a = summarize("p3_1088996", "P3  Trp87 / 86-89 patch")
b = summarize("p2_1088883", "P2  24/25/27/42 patch")
if a and b:
    print("\n" + "=" * 62)
    print(f"{'':22s} {'P3 (86-89)':>14s} {'P2 (24-42)':>14s}")
    for key, lab in (("n", "refolds"), ("u2", "motif < 2A"),
                     ("h4", "4/4 hotspots"), ("gate", "PASS gate"),
                     ("elite", "ELITE")):
        pa = 100 * a[key] / a["n"]; pb = 100 * b[key] / b["n"]
        if key == "n":
            print(f"{lab:22s} {a[key]:14d} {b[key]:14d}")
        else:
            print(f"{lab:22s} {a[key]:6d} ({pa:3.0f}%) {b[key]:6d} ({pb:3.0f}%)")
    d = 100 * a["gate"] / a["n"] - 100 * b["gate"] / b["n"]
    print(f"\ngate survival difference: {d:+.0f} percentage points for P3")
    # a crude binomial check on whether that gap is real
    import math
    p = (a["gate"] + b["gate"]) / (a["n"] + b["n"])
    se = math.sqrt(p * (1 - p) * (1 / a["n"] + 1 / b["n"])) * 100
    print(f"pooled SE of the difference: {se:.1f} pp  -> "
          f"{'significant' if abs(d) > 2*se else 'NOT distinguishable'} at ~2 SE")
    out = os.path.join(analysis, "p3_vs_p2_summary.json")
    json.dump({k: {kk: v[kk] for kk in ("label","n","gate","elite","h4","u2")}
               for k, v in (("p3", a), ("p2", b))}, open(out, "w"), indent=2)
    print(f"\nwrote {out}")
