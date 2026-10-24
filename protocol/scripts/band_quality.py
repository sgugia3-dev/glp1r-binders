#!/usr/bin/env python
"""Full quality profile of the q(11-140) +4..+9 band, against the other bands.

Joins every available measurement for these designs: confidence (ipSAE, pLDDT,
ipTM, pDockQ), geometry (motif and receptor RMSD, hotspot contacts), interaction
counts (salt bridges, H-bonds, hydrophobic contacts) and sequence properties
(pI, GRAVY, class fractions, aliphatic and instability indices).

Counts are reported as backbones as well as folds, since MPNN made 3 sequences per
input backbone. The ranked table only holds GATED designs, so confidence metrics
are available for the gated subset while geometry and interactions cover all.
"""
import csv, glob, os, re, statistics as st, sys
from collections import Counter, defaultdict
from Bio.SeqUtils.ProtParam import ProteinAnalysis

A, SCRATCH, RUN, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
LO, HI = 4, 9
CLASSES = {"hydrophobic": set("AVLIMFWC"), "polar": set("STNQYH"),
           "positive": set("KR"), "negative": set("DE"), "special": set("GP")}


def load(pat, key="name"):
    out = {}
    for f in sorted(glob.glob(os.path.join(A, pat))):
        m = re.search(r"_sh(\d+)[_.]", f)
        sh = m.group(1) if m else "?"
        for r in csv.DictReader(open(f)):
            out[(sh, r[key])] = r
    return out


dev, inter, epi = load(f"{RUN}_sh*_deviation.csv"), \
                  load(f"{RUN}_sh*_interactions.csv"), load(f"{RUN}_sh*_epitope.csv")
ranked = {}
for r in csv.DictReader(open(os.path.join(A, "p14_RANKED.csv"))):
    ranked[(r["shard"], r["name"])] = r
seqs = {}
for y in glob.glob(os.path.join(SCRATCH, "runs/glp1r", RUN, "sh*", "yaml", "*.yaml")):
    sh = re.search(r"/sh(\d+)/", y).group(1)
    t = open(y).read()
    m = re.search(r"id:\s*A\s*\n\s*sequence:\s*([A-Z]+)", t)
    if m: seqs[(sh, os.path.basename(y)[:-5])] = m.group(1)
print(f"deviation {len(dev)}  interactions {len(inter)}  epitope {len(epi)}  "
      f"ranked {len(ranked)}  sequences {len(seqs)}")

rows = []
for k in sorted(set(dev) & set(inter) & set(epi)):
    d, i, e = dev[k], inter[k], epi[k]
    try:
        qd = int(i["designed_charge"])
        nt = int(e.get("n_hotspots_total") or 4)
        r = {"shard": k[0], "name": k[1], "src": e.get("source", "?"),
             "bb": (k[0], re.sub(r"__s\d+$", "", k[1])),
             "q_des": qd, "motif": float(d["motif_rmsd_vs_af3"]),
             "recept": float(d["receptor_rmsd"]), "hits": int(e["n_hotspots_hit"]),
             "n_total": nt, "full": int(e["n_hotspots_hit"]) == nt,
             "mind": float(e["min_hotspot_dist"] or 99),
             "sb": int(i.get("unique_salt_bridges", 0) or 0),
             "sb_des": int(i.get("sb_designed", 0) or 0),
             "hb": int(i.get("hbonds", 0) or 0),
             "hc": int(i.get("hydrophobic_contacts", 0) or 0),
             "anchor": i.get("anchor_seq", ""), "seq": seqs.get(k, "")}
        rk = ranked.get(k)
        for f_ in ("ipsae", "ipsae_cons", "plddt_binder", "plddt_design",
                   "plddt_complex", "iplddt", "iptm", "pdockq", "score"):
            r[f_] = float(rk[f_]) if rk and rk.get(f_) else None
        r["gated"] = rk is not None
        rows.append(r)
    except (ValueError, KeyError):
        pass

band = [r for r in rows if LO <= r["q_des"] < HI]
mid = [r for r in rows if HI <= r["q_des"] < 16]
high = [r for r in rows if r["q_des"] >= 16]
nb = lambda s: len({r["bb"] for r in s})
print(f"\nband q(11-140) {LO}..{HI}: {len(band)} folds / {nb(band)} backbones")


def dist(sel, key, fmt="{:6.3f}"):
    v = [r[key] for r in sel if r.get(key) is not None]
    if not v: return "        n/a"
    s = sorted(v)
    return (fmt + " " + fmt + " " + fmt).format(s[len(s)//10], st.median(s), s[9*len(s)//10])


print("\n" + "=" * 92)
print(f"QUALITY METRICS, q(11-140) {LO}..{HI}   (p10 / median / p90)")
print("=" * 92)
print(f"{'metric':22s} {'p10':>8s} {'median':>8s} {'p90':>8s}   {'mid band':>26s}")
for key, fmt in (("motif", "{:6.2f}"), ("recept", "{:6.2f}"), ("mind", "{:6.2f}"),
                 ("ipsae", "{:6.3f}"), ("ipsae_cons", "{:6.3f}"),
                 ("plddt_binder", "{:6.3f}"), ("plddt_design", "{:6.3f}"),
                 ("iplddt", "{:6.3f}"), ("iptm", "{:6.3f}"), ("pdockq", "{:6.3f}"),
                 ("sb", "{:6.0f}"), ("sb_des", "{:6.0f}"), ("hb", "{:6.0f}"),
                 ("hc", "{:6.0f}")):
    print(f"{key:22s} {dist(band, key, fmt)}   {dist(mid, key, fmt):>26s}")

print("\nGEOMETRY AND EPITOPE")
g = [r for r in band if r["motif"] < 2.0 and r["hits"] >= 3 and r["recept"] < 5.0]
f_ = [r for r in g if r["full"]]
print(f"  motif < 2 A        {sum(1 for r in band if r['motif']<2):4d}/{len(band)} "
      f"({100*sum(1 for r in band if r['motif']<2)/len(band):3.0f}%)")
print(f"  full epitope       {sum(1 for r in band if r['full']):4d}/{len(band)} "
      f"({100*sum(1 for r in band if r['full'])/len(band):3.0f}%)")
print(f"  pass gate          {len(g):4d}/{len(band)} ({100*len(g)/len(band):3.0f}%)   "
      f"backbones {nb(g)}")
print(f"  gate + full        {len(f_):4d}/{len(band)} ({100*len(f_)/len(band):3.0f}%)   "
      f"backbones {nb(f_)}")
print(f"  hits histogram     {dict(Counter(r['hits'] for r in band))}")

print("\nBY SOURCE CAMPAIGN")
print(f"  {'src':5s} {'folds':>6s} {'bb':>4s} {'gate+full':>10s} {'ipSAE':>7s} "
      f"{'pLDDTb':>7s} {'saltbr':>7s} {'motif':>7s}")
for s_ in ("p3", "p10", "p12", "p13"):
    sub = [r for r in band if r["src"] == s_]
    if not sub: continue
    gf = [r for r in sub if r["motif"] < 2 and r["hits"] >= 3
          and r["recept"] < 5 and r["full"]]
    ip = [r["ipsae"] for r in sub if r["ipsae"] is not None]
    pl = [r["plddt_binder"] for r in sub if r["plddt_binder"] is not None]
    print(f"  {s_:5s} {len(sub):6d} {nb(sub):4d} {len(gf):5d} ({100*len(gf)/len(sub):3.0f}%) "
          f"{(st.median(ip) if ip else 0):7.3f} {(st.median(pl) if pl else 0):7.3f} "
          f"{st.median([r['sb'] for r in sub]):7.0f} "
          f"{st.median([r['motif'] for r in sub]):7.2f}")

print("\nSEQUENCE PROPERTIES of the band (full 140-aa chain)")
props = []
for r in band:
    if not r["seq"]: continue
    p = ProteinAnalysis(r["seq"]); c = Counter(r["seq"]); n = len(r["seq"])
    props.append({"pI": p.isoelectric_point(), "gravy": p.gravy(),
                  "q74": p.charge_at_pH(7.4), "instab": p.instability_index(),
                  "arom": p.aromaticity(),
                  **{f"f_{k}": sum(c[a] for a in v)/n for k, v in CLASSES.items()}})
if props:
    print(f"  {'property':16s} {'p10':>8s} {'median':>8s} {'p90':>8s}")
    for k in ("pI", "q74", "gravy", "instab", "arom", "f_hydrophobic", "f_polar",
              "f_positive", "f_negative", "f_special"):
        v = sorted(x[k] for x in props)
        print(f"  {k:16s} {v[len(v)//10]:8.3f} {st.median(v):8.3f} {v[9*len(v)//10]:8.3f}")

print("\nTOP 12 IN THE BAND by composite score (gated only)")
gs = sorted([r for r in band if r["score"] is not None], key=lambda r: -r["score"])
print(f"  {'score':>6s} {'design':40s} {'src':4s} {'q':>4s} {'motif':>6s} {'hit':>4s} "
      f"{'ipSAE':>6s} {'pLDDTb':>7s} {'sb':>3s} {'hc':>4s}")
for r in gs[:12]:
    print(f"  {r['score']:6.3f} {r['name'][:40]:40s} {r['src']:4s} {r['q_des']:+4d} "
          f"{r['motif']:6.2f} {r['hits']:d}/{r['n_total']:d}  {r['ipsae']:6.3f} "
          f"{r['plddt_binder']:7.3f} {r['sb']:3d} {r['hc']:4d}")

cols = ["shard", "name", "src", "bb_tag", "q_des", "motif", "recept", "hits",
        "n_total", "full", "mind", "sb", "sb_des", "hb", "hc", "gated", "score",
        "ipsae", "ipsae_cons", "plddt_binder", "plddt_design", "plddt_complex",
        "iplddt", "iptm", "pdockq", "anchor", "seq"]
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(cols)
    for r in sorted(band, key=lambda x: -(x["score"] or -9)):
        r2 = dict(r); r2["bb_tag"] = r["bb"][1]
        w.writerow([r2.get(c, "") for c in cols])
print(f"\nwrote {OUT}  ({len(band)} rows)")
