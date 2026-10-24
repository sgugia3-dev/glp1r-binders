#!/usr/bin/env python
"""Aggregate the charge-complementary redesign round against its own inputs.

The 200 inputs were all FULLY engaged designs drawn from p3/p10/p12/p13_best, so
"retention" here is the fraction that still contact every hotspot in their source
campaign's set after being redesigned. That is a stricter and more informative
question than a fresh campaign's hit rate.

Also merges the per-source epitope tables into one rank-compatible file per shard:
the four source campaigns use three different hotspot sets, so h5_epitope wrote
them separately, but the ranker needs a single {RUN}_sh{N}_epitope.csv.
"""
import csv, glob, json, os, re, statistics as st, sys
from collections import Counter, defaultdict

A, RUN = sys.argv[1], sys.argv[2]
SRC_HOTS = {"p3": 4, "p10": 3, "p12": 3, "p13": 4}
COMMON = ["name", "binder_len", "motif_rmsd", "receptor_rmsd",
          "n_hotspots_hit", "min_hotspot_dist", "designed_engages"]


def load(pat):
    out = {}
    for f in sorted(glob.glob(os.path.join(A, pat))):
        sh = re.search(r"_sh(\d+)[_.]", f).group(1)
        for r in csv.DictReader(open(f)):
            out[(sh, r["name"])] = r
    return out


dev = load(f"{RUN}_sh*_deviation.csv")
inter = load(f"{RUN}_sh*_interactions.csv")
epi, src_of = {}, {}
for pfx in SRC_HOTS:
    for f in sorted(glob.glob(os.path.join(A, f"{RUN}_sh*_epitope-{pfx}.csv"))):
        sh = re.search(r"_sh(\d+)_", f).group(1)
        for r in csv.DictReader(open(f)):
            epi[(sh, r["name"])] = r
            src_of[(sh, r["name"])] = pfx
print(f"deviation {len(dev)}  epitope {len(epi)}  interactions {len(inter)}")

# merged, rank-compatible epitope table per shard
by_shard = defaultdict(list)
for (sh, nm), r in epi.items():
    row = {k: r.get(k, "") for k in COMMON}
    row["source"] = src_of[(sh, nm)]
    row["n_hotspots_total"] = SRC_HOTS[src_of[(sh, nm)]]
    by_shard[sh].append(row)
for sh, rows in by_shard.items():
    out = os.path.join(A, f"{RUN}_sh{sh}_epitope.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COMMON + ["source", "n_hotspots_total"])
        w.writeheader(); w.writerows(sorted(rows, key=lambda r: r["name"]))
print(f"wrote {len(by_shard)} merged epitope tables for the ranker")

keys = sorted(set(dev) & set(epi))
rows = []
for k in keys:
    d, e = dev[k], epi[k]
    i = inter.get(k, {})
    try:
        nh = SRC_HOTS[src_of[k]]
        rows.append({
            "shard": k[0], "name": k[1], "source": src_of[k], "n_total": nh,
            "motif": float(d["motif_rmsd_vs_af3"]), "recept": float(d["receptor_rmsd"]),
            "hits": int(e["n_hotspots_hit"]), "full": int(e["n_hotspots_hit"]) == nh,
            "sb": int(i.get("unique_salt_bridges", 0) or 0),
            "sb_des": int(i.get("sb_designed", 0) or 0),
            "hb": int(i.get("hbonds", 0) or 0),
            "hc": int(i.get("hydrophobic_contacts", 0) or 0),
            "q_des": int(i.get("designed_charge", 0) or 0),
            "anchor": i.get("anchor_seq", ""),
        })
    except (ValueError, KeyError):
        pass
n = len(rows)
print(f"joined {n}\n")

print("=" * 86)
print(f"P14 charge-complementary redesign   ({n} refolds from 200 input backbones)")
print("=" * 86)
u2 = [r for r in rows if r["motif"] < 2.0]
fu = [r for r in rows if r["full"]]
g3 = [r for r in rows if r["hits"] >= 3]
gate = [r for r in rows if r["motif"] < 2.0 and r["hits"] >= 3 and r["recept"] < 5.0]
print(f"  motif under 2A          : {len(u2):4d} ({100*len(u2)/n:3.0f}%)")
print(f"  FULL epitope retained   : {len(fu):4d} ({100*len(fu)/n:3.0f}%)")
print(f"  >=3 hotspots            : {len(g3):4d} ({100*len(g3)/n:3.0f}%)")
print(f"  PASS gate               : {len(gate):4d} ({100*len(gate)/n:3.0f}%)")
print(f"  complete misses         : {sum(1 for r in rows if r['hits']==0):4d}")

print("\nRETENTION BY SOURCE CAMPAIGN (inputs were 100% fully engaged)")
print(f"  {'source':8s} {'n':>5s} {'full':>12s} {'>=3':>10s} {'motif<2A':>10s}")
for pfx in ("p3", "p10", "p12", "p13"):
    sub = [r for r in rows if r["source"] == pfx]
    if not sub: continue
    f_ = sum(1 for r in sub if r["full"])
    print(f"  {pfx:8s} {len(sub):5d} {f_:5d} ({100*f_/len(sub):3.0f}%) "
          f"{sum(1 for r in sub if r['hits']>=3):5d} ({100*sum(1 for r in sub if r['hits']>=3)/len(sub):3.0f}%) "
          f"{100*sum(1 for r in sub if r['motif']<2)/len(sub):8.0f}%")

print("\nINTERACTIONS (medians; the H-bond count is distance-only, no angle term)")
for key, lbl in (("sb", "unique salt bridges"), ("sb_des", "  of those designed-region"),
                 ("hb", "H-bond atom pairs"), ("hc", "hydrophobic C-C contacts")):
    v = [r[key] for r in rows]
    print(f"  {lbl:28s} median {st.median(v):5.0f}   p10 {sorted(v)[len(v)//10]:4.0f}  "
          f"p90 {sorted(v)[9*len(v)//10]:4.0f}  max {max(v):4.0f}")

print("\nCHARGE of the refolded designs")
q = sorted(r["q_des"] for r in rows)
print(f"  q(11-140)  median {st.median(q):+.0f}   p10 {q[len(q)//10]:+.0f}  "
      f"p90 {q[9*len(q)//10]:+.0f}   min {q[0]:+.0f}  max {q[-1]:+.0f}")
print(f"  pre-redesign baseline was -5.7 median")
hi = [r for r in rows if r["q_des"] > 16]
print(f"  designs with q(11-140) > +16 (the high tail): {len(hi)} ({100*len(hi)/n:.0f}%)")

print("\nREDESIGNED ANCHORS")
anc = Counter(r["anchor"] for r in rows if r["anchor"])
print(f"  distinct anchors: {len(anc)} across {sum(anc.values())} designs")
ok1 = sum(1 for r in rows if r["anchor"][:1] == "H")
ok3 = sum(1 for r in rows if r["anchor"][2:3] in ("D", "E"))
print(f"  position 1 == H : {ok1}/{n}")
print(f"  position 3 D/E  : {ok3}/{n}")
for s, c in anc.most_common(6):
    print(f"    {s}  x{c}")

print("\nDOES CHARGE TRACK INTERACTIONS OR EPITOPE RETENTION?")
for lo, hi_ in ((None, 4), (4, 9), (9, 16), (16, None)):
    sub = [r for r in rows
           if (lo is None or r["q_des"] >= lo) and (hi_ is None or r["q_des"] < hi_)]
    if len(sub) < 10: continue
    lbl = f"q {lo if lo is not None else '-inf'}..{hi_ if hi_ is not None else 'inf'}"
    print(f"  {lbl:14s} n={len(sub):4d}  salt bridges {st.median([r['sb'] for r in sub]):4.0f}  "
          f"full epitope {100*sum(1 for r in sub if r['full'])/len(sub):3.0f}%  "
          f"motif<2A {100*sum(1 for r in sub if r['motif']<2)/len(sub):3.0f}%")

json.dump({"n": n, "full": len(fu), "gate": len(gate),
           "sb_median": st.median([r["sb"] for r in rows]),
           "q_median": st.median(q)},
          open(os.path.join(A, f"{RUN}_summary.json"), "w"), indent=2)
print(f"\nwrote {RUN}_summary.json")
