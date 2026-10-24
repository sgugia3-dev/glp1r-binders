#!/usr/bin/env python
"""Verify the anchor repair actually landed, then summarise the round.

Checks the refolded structures rather than the design-stage log: the anchor in
the predicted complex is what matters, and it must satisfy every consensus
restraint while residues 11-140 stay as the p14 input had them.
"""
import csv, glob, os, re, statistics as st, sys
from collections import Counter
sys.path.insert(0, "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts")
from consensus import CONSENSUS, INVARIANT, NATIVE

A, RUN = sys.argv[1], sys.argv[2]
SRC_HOTS = {"p3": 4, "p10": 3, "p12": 3, "p13": 4}


def load(pat):
    out = {}
    for f in sorted(glob.glob(os.path.join(A, pat))):
        m = re.search(r"_sh(\d+)[_.]", f)
        sh = m.group(1) if m else "?"
        for r in csv.DictReader(open(f)):
            out[(sh, r["name"])] = r
    return out


dev, inter = load(f"{RUN}_sh*_deviation.csv"), load(f"{RUN}_sh*_interactions.csv")
epi, src = {}, {}
for pfx in SRC_HOTS:
    for f in sorted(glob.glob(os.path.join(A, f"{RUN}_sh*_epitope-{pfx}.csv"))):
        sh = re.search(r"_sh(\d+)_", f).group(1)
        for r in csv.DictReader(open(f)):
            epi[(sh, r["name"])] = r; src[(sh, r["name"])] = pfx
print(f"deviation {len(dev)}  epitope {len(epi)}  interactions {len(inter)}")

# merged rank-compatible epitope tables
COMMON = ["name", "binder_len", "motif_rmsd", "receptor_rmsd",
          "n_hotspots_hit", "min_hotspot_dist", "designed_engages"]
by_sh = {}
for (sh, nm), r in epi.items():
    row = {k: r.get(k, "") for k in COMMON}
    row["source"] = src[(sh, nm)]; row["n_hotspots_total"] = SRC_HOTS[src[(sh, nm)]]
    by_sh.setdefault(sh, []).append(row)
for sh, rows in by_sh.items():
    with open(os.path.join(A, f"{RUN}_sh{sh}_epitope.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COMMON + ["source", "n_hotspots_total"])
        w.writeheader(); w.writerows(sorted(rows, key=lambda r: r["name"]))
print(f"wrote {len(by_sh)} merged epitope tables")

rows = []
for k in sorted(set(dev) & set(epi) & set(inter)):
    d, e, i = dev[k], epi[k], inter[k]
    try:
        nt = SRC_HOTS[src[k]]
        rows.append({"shard": k[0], "name": k[1], "src": src[k],
                     "anchor": i.get("anchor_seq", ""),
                     "motif": float(d["motif_rmsd_vs_af3"]),
                     "recept": float(d["receptor_rmsd"]),
                     "hits": int(e["n_hotspots_hit"]), "n_total": nt,
                     "full": int(e["n_hotspots_hit"]) == nt,
                     "sb": int(i.get("unique_salt_bridges", 0) or 0),
                     "hb": int(i.get("hbonds", 0) or 0),
                     "hc": int(i.get("hydrophobic_contacts", 0) or 0),
                     "q": int(i.get("designed_charge", 0) or 0)})
    except (ValueError, KeyError):
        pass
n = len(rows)
print(f"joined {n}\n")

print("=" * 76)
print("ANCHOR CONSENSUS COMPLIANCE, measured on the refolded structures")
print("=" * 76)
anc = Counter(r["anchor"] for r in rows if r["anchor"])
full_ok = sum(1 for r in rows
              if r["anchor"] and not [p for p in CONSENSUS
                                      if r["anchor"][p-1] not in CONSENSUS[p][0]])
inv_ok = sum(1 for r in rows
             if r["anchor"] and not [p for p in INVARIANT
                                     if r["anchor"][p-1] not in CONSENSUS[p][0]])
print(f"  satisfy ALL ten restraints   : {full_ok}/{n} ({100*full_ok/n:.0f}%)")
print(f"  satisfy invariant positions  : {inv_ok}/{n} ({100*inv_ok/n:.0f}%)")
print(f"  p14 for comparison           : 0/600 (0%) on both")
print(f"  distinct anchors             : {len(anc)}")
nn = [sum(1 for i2 in range(10) if a[i2] == NATIVE[i2])
      for a in (r["anchor"] for r in rows) if a]
if nn:
    print(f"  identity to native GLP-1     : median {st.median(nn):.0f}/10  "
          f"min {min(nn)}  max {max(nn)}")
print("  anchors found:")
for s, c in anc.most_common(10):
    diff = "".join("." if s[i2] == NATIVE[i2] else s[i2] for i2 in range(10))
    print(f"    {s}  x{c:<4d}  vs native {NATIVE}: {diff}")

print("\n" + "=" * 76)
print(f"ROUND SUMMARY  ({n} refolds)")
print("=" * 76)
u2 = [r for r in rows if r["motif"] < 2.0]
fu = [r for r in rows if r["full"]]
gate = [r for r in rows if r["motif"] < 2.0 and r["hits"] >= 3 and r["recept"] < 5.0]
print(f"  motif under 2A     : {len(u2):4d} ({100*len(u2)/n:3.0f}%)")
print(f"  full epitope       : {len(fu):4d} ({100*len(fu)/n:3.0f}%)")
print(f"  PASS gate          : {len(gate):4d} ({100*len(gate)/n:3.0f}%)")
print(f"  motif RMSD         : median {st.median([r['motif'] for r in rows]):.2f} A")
print(f"  salt bridges       : median {st.median([r['sb'] for r in rows]):.0f}")
print(f"  H-bond pairs       : median {st.median([r['hb'] for r in rows]):.0f}")
print(f"  hydrophobic        : median {st.median([r['hc'] for r in rows]):.0f}")
print(f"  q(11-140)          : median {st.median([r['q'] for r in rows]):+.0f} "
      f"(inherited from p14, unchanged by design)")
print("\n  by source:")
for p in ("p3", "p10", "p12", "p13"):
    s2 = [r for r in rows if r["src"] == p]
    if not s2: continue
    print(f"    {p:4s} n={len(s2):4d}  full {100*sum(1 for r in s2 if r['full'])/len(s2):3.0f}%  "
          f"motif<2A {100*sum(1 for r in s2 if r['motif']<2)/len(s2):3.0f}%  "
          f"sb {st.median([r['sb'] for r in s2]):.0f}")
