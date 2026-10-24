#!/usr/bin/env python
"""Collect ipSAE and join it to the existing deviation/epitope tables.

Project rule: ipSAE at cutoffs 15/15 is the ranking metric and ipTM must not be
substituted. The `max` row is the binder-scored (B->A) direction and is the more
optimistic of the two asymmetric values, so the target-scored direction is reported
alongside it rather than hidden.
"""
import csv, glob, json, os, re, statistics as st, sys

root, analysis, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]

def parse(p):
    """Return {type: row} from an ipsae *_15_15.txt table."""
    rows, hdr = {}, None
    for line in open(p):
        t = line.split()
        if not t: continue
        if t[0] == "Chn1": hdr = t; continue
        if hdr and len(t) >= len(hdr) - 1:
            d = dict(zip(hdr, t))
            ty = d.get("Type")
            if ty: rows[ty] = d
    return rows

# deviation + epitope, keyed (shard, name)
def load(pat, key):
    out = {}
    for f in sorted(glob.glob(os.path.join(analysis, pat))):
        m = re.search(r"_sh(\d+)_", f); sh = m.group(1) if m else "?"
        for r in csv.DictReader(open(f)): out[(sh, r["name"])] = r
    return out
dev = load("p1_*_deviation.csv", "name")
epi = load("p1_*_epitope.csv", "name")

rows, nofile = [], 0
for txt in glob.glob(os.path.join(root, "sh*", "**", "predictions", "*", "*_15_15.txt"),
                     recursive=True):
    name = os.path.basename(os.path.dirname(txt))
    m = re.search(r"/sh(\d+)/", txt); sh = m.group(1) if m else "?"
    tb = parse(txt)
    mx = tb.get("max")
    if not mx: nofile += 1; continue
    def g(d, k):
        try: return float(d.get(k, "nan"))
        except ValueError: return float("nan")
    d, e = dev.get((sh, name), {}), epi.get((sh, name), {})
    rows.append({"shard": sh, "name": name,
                 "ipsae_max": g(mx, "ipSAE"),
                 "ipsae_asym": min(g(v, "ipSAE") for v in tb.values()) if tb else float("nan"),
                 "pdockq": g(mx, "pDockQ"), "lis": g(mx, "LIS"),
                 "iptm_af": g(mx, "ipTM_af"),
                 "motif_rmsd": d.get("motif_rmsd_vs_af3", ""),
                 "hotspots_hit": e.get("n_hotspots_hit", ""),
                 "binder_len": d.get("binder_len", "")})
print(f"ipsae tables parsed: {len(rows)}  (no max row: {nofile})")
if not rows: sys.exit(1)
cols = list(rows[0].keys())
with open(out_csv, "w") as fh:
    fh.write(",".join(cols) + "\n")
    for r in rows: fh.write(",".join(str(r[c]) for c in cols) + "\n")

v = [r["ipsae_max"] for r in rows if r["ipsae_max"] == r["ipsae_max"]]
a = [r["ipsae_asym"] for r in rows if r["ipsae_asym"] == r["ipsae_asym"]]
print(f"\nipSAE (max / binder-scored): min {min(v):.3f} med {st.median(v):.3f} max {max(v):.3f}")
print(f"ipSAE (conservative dirn)  : min {min(a):.3f} med {st.median(a):.3f} max {max(a):.3f}")
for t in (0.3, 0.5, 0.7, 0.8):
    print(f"  >= {t}: {sum(1 for x in v if x >= t):5d}/{len(v)}")
# does ipSAE track epitope adherence? the recorded warning says it need not
both = [r for r in rows if r["hotspots_hit"] not in ("", None)]
if both:
    for k in (4, 3, 0):
        s = [r["ipsae_max"] for r in both if str(r["hotspots_hit"]) == str(k)]
        if len(s) >= 5:
            print(f"  hotspots={k}: n={len(s):4d}  ipSAE median {st.median(s):.3f}")
top = sorted(rows, key=lambda r: -r["ipsae_max"])[:10]
print(f"\n{'design':34s} {'sh':>3s} {'ipSAE':>7s} {'pDockQ':>7s} {'hits':>5s} {'motif':>7s}")
for r in top:
    print(f"{r['name'].replace('glp1r_ecd4_glp1r_','')[:34]:34s} {r['shard']:>3s} "
          f"{r['ipsae_max']:7.3f} {r['pdockq']:7.3f} {str(r['hotspots_hit']):>5s} {str(r['motif_rmsd'])[:6]:>7s}")
print(f"\nwrote {out_csv}")
