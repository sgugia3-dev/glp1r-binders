#!/usr/bin/env python
"""Aggregate the charge-complementary redesign round, discarding the high-charge tail.

Charge is taken from the REFOLDED binder sequence, not from the design-stage
record, so the number reflects what was actually predicted. Full-chain charge is
the designed-region charge reported by h7 plus the charge of the redesigned anchor.

Designs with full-chain net charge above +10 are removed before any statistics are
computed: the per-shard medians were on target at +5 to +7, but the tails reached
+22 to +33 where selection traded charge for contacts.

Epitope results arrive in four files per shard, one per source campaign, because the
200 input backbones carry three different hotspot sets. They are merged on name and
each design keeps its own set's hotspot count.
"""
import csv, glob, json, os, re, statistics as st, sys
from collections import Counter, defaultdict

A, RUN, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
MAX_Q = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0
SRC_HOTSPOTS = {"p3": 4, "p10": 3, "p12": 3, "p13": 4}


def q(seq):
    c = Counter(seq)
    return c["K"] + c["R"] - c["D"] - c["E"]


def load(pat):
    out = {}
    for f in sorted(glob.glob(os.path.join(A, pat))):
        sh = re.search(r"_sh(\d+)[_.]", f).group(1)
        for r in csv.DictReader(open(f)):
            out[(sh, r["name"])] = r
    return out


dev = load(f"{RUN}_sh*_deviation.csv")
inter = load(f"{RUN}_sh*_interactions.csv")
epi = {}
for pfx in SRC_HOTSPOTS:
    for f in sorted(glob.glob(os.path.join(A, f"{RUN}_sh*_{pfx}_epitope.csv"))):
        sh = re.search(r"_sh(\d+)_", f).group(1)
        for r in csv.DictReader(open(f)):
            epi[(sh, r["name"])] = dict(r, _src=pfx)
print(f"deviation {len(dev)}  interactions {len(inter)}  epitope {len(epi)}")

rows = []
for k, iv in inter.items():
    dv, ev = dev.get(k), epi.get(k)
    if not dv or not ev: continue
    anc = iv.get("anchor_seq", "")
    try:
        des_q = int(iv["designed_charge"])
    except (KeyError, ValueError):
        continue
    rows.append({
        "shard": k[0], "name": k[1], "src": ev["_src"],
        "n_hot_total": SRC_HOTSPOTS[ev["_src"]],
        "q_full": des_q + q(anc), "q_designed": des_q, "q_anchor": q(anc),
        "anchor_seq": anc,
        "salt_bridges": int(iv["unique_salt_bridges"]),
        "sb_designed": int(iv["sb_designed"]),
        "hbonds": int(iv["hbonds"]),
        "hydrophobic": int(iv["hydrophobic_contacts"]),
        "motif": float(dv["motif_rmsd_vs_af3"]),
        "recept": float(dv["receptor_rmsd"]),
        "hits": int(ev["n_hotspots_hit"]),
    })
n_all = len(rows)
print(f"joined {n_all} refolds\n")

kept = [r for r in rows if r["q_full"] <= MAX_Q]
drop = [r for r in rows if r["q_full"] > MAX_Q]
print("=" * 86)
print(f"CHARGE FILTER: removing full-chain net charge above +{MAX_Q:g}")
print("=" * 86)
print(f"  removed : {len(drop):4d}/{n_all} ({100*len(drop)/max(n_all,1):.0f}%)")
if drop:
    dq = sorted(r["q_full"] for r in drop)
    print(f"            their charge: median {st.median(dq):+.0f}  "
          f"min {dq[0]:+.0f}  max {dq[-1]:+.0f}")
print(f"  kept    : {len(kept):4d}/{n_all} ({100*len(kept)/max(n_all,1):.0f}%)")
if not kept:
    print("nothing survived the filter"); sys.exit(1)
kq = sorted(r["q_full"] for r in kept)
print(f"            their charge: median {st.median(kq):+.0f}  "
      f"min {kq[0]:+.0f}  max {kq[-1]:+.0f}")
pos = sum(1 for r in kept if r["q_full"] > 0)
print(f"            net positive: {pos}/{len(kept)} ({100*pos/len(kept):.0f}%)")


def tab(sel, lbl):
    f = lambda key: st.median([r[key] for r in sel])
    print(f"{lbl:22s} {len(sel):5d} {f('q_full'):+7.0f} {f('q_designed'):+8.0f} "
          f"{f('q_anchor'):+8.0f} {f('salt_bridges'):8.0f} {f('hbonds'):7.0f} "
          f"{f('hydrophobic'):8.0f} {f('motif'):8.2f}")


print("\n" + "=" * 86)
print("KEPT DESIGNS, medians")
print("=" * 86)
print(f"{'group':22s} {'n':>5s} {'q_full':>7s} {'q_des':>8s} {'q_anch':>8s} "
      f"{'saltbr':>8s} {'hbond':>7s} {'hydroph':>8s} {'motifA':>8s}")
tab(kept, "ALL")
for pfx in ("p3", "p10", "p12", "p13"):
    sel = [r for r in kept if r["src"] == pfx]
    if sel: tab(sel, f"  from {pfx}")

print("\n" + "=" * 86)
print("EPITOPE RETENTION after redesign, per source campaign")
print("=" * 86)
print(f"{'source':10s} {'nhot':>5s} {'n':>5s} {'full':>11s} {'>=3':>11s} "
      f"{'miss':>9s} {'motif<2A':>10s}")
BEFORE = {"p3": 23, "p10": 82, "p12": 75, "p13": 39}
for pfx in ("p3", "p10", "p12", "p13"):
    sel = [r for r in kept if r["src"] == pfx]
    if not sel: continue
    nh = SRC_HOTSPOTS[pfx]
    fu = sum(1 for r in sel if r["hits"] == nh)
    g3 = sum(1 for r in sel if r["hits"] >= 3)
    ms = sum(1 for r in sel if r["hits"] == 0)
    u2 = sum(1 for r in sel if r["motif"] < 2.0)
    print(f"{pfx:10s} {nh:5d} {len(sel):5d} {100*fu/len(sel):6.0f}% "
          f"(was {BEFORE[pfx]:2d}%) {100*g3/len(sel):10.0f}% "
          f"{100*ms/len(sel):8.0f}% {100*u2/len(sel):9.0f}%")

print("\n" + "=" * 86)
print("INTERACTIONS vs the pre-redesign baseline (P12/P13 best sets, 250 complexes)")
print("=" * 86)
print("  baseline binder designed face: charge +2, 39-42 interface residues")
print(f"  redesigned salt bridges : median {st.median([r['salt_bridges'] for r in kept]):.0f}  "
      f"max {max(r['salt_bridges'] for r in kept)}")
print(f"    of those in 11-140    : median {st.median([r['sb_designed'] for r in kept]):.0f}")
print(f"  H-bond pairs (<3.5A)    : median {st.median([r['hbonds'] for r in kept]):.0f}  "
      f"(distance only, no angle term)")
print(f"  hydrophobic C-C (<4.5A) : median {st.median([r['hydrophobic'] for r in kept]):.0f}")
print(f"  motif RMSD              : median {st.median([r['motif'] for r in kept]):.2f} A  "
      f"(fixed-anchor campaigns: 1.28-1.62 A)")

anc = Counter(r["anchor_seq"] for r in kept)
print(f"\n  distinct redesigned anchors : {len(anc)} over {len(kept)} designs")
print(f"  His at position 1           : "
      f"{sum(1 for r in kept if r['anchor_seq'][:1]=='H')}/{len(kept)}")
print(f"  Asp/Glu at position 3       : "
      f"{sum(1 for r in kept if r['anchor_seq'][2:3] in ('D','E'))}/{len(kept)}")
print("  most common anchors:")
for s2, n2 in anc.most_common(6):
    print(f"      {s2}  x{n2}   (native was HAEGTFTSDV)")

cols = ["shard", "name", "src", "n_hot_total", "hits", "q_full", "q_designed",
        "q_anchor", "anchor_seq", "salt_bridges", "sb_designed", "hbonds",
        "hydrophobic", "motif", "recept"]
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
    w.writeheader(); w.writerows(sorted(kept, key=lambda r: (r["src"], r["name"])))
print(f"\nwrote {OUT}  ({len(kept)} designs passing the charge filter)")
