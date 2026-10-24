#!/usr/bin/env python
"""How many candidates survive in the low-charge group?

Two definitions are in play, so both are reported:
  q(full 140)  the controlled quantity, target +5.5 +-1.5 -> the window is +4..+7
  q(11-140)    the designed region, which the earlier bin table used (+4..+9)

Counts are given as BACKBONES as well as folds: ProteinMPNN produced 3 sequences
per input backbone, so a fold count overstates distinct candidates ~3x.
"""
import csv, glob, os, re, sys
from collections import Counter, defaultdict

A, RUN = sys.argv[1], sys.argv[2]
SRC_HOTS = {"p3": 4, "p10": 3, "p12": 3, "p13": 4}


def load(pat):
    out = {}
    for f in sorted(glob.glob(os.path.join(A, pat))):
        sh = re.search(r"_sh(\d+)[_.]", f).group(1)
        for r in csv.DictReader(open(f)):
            out[(sh, r["name"])] = r
    return out


dev, inter, epi = load(f"{RUN}_sh*_deviation.csv"), \
                  load(f"{RUN}_sh*_interactions.csv"), \
                  load(f"{RUN}_sh*_epitope.csv")
q = lambda s: Counter(s)["K"] + Counter(s)["R"] - Counter(s)["D"] - Counter(s)["E"]

rows = []
for k in sorted(set(dev) & set(inter) & set(epi)):
    d, i, e = dev[k], inter[k], epi[k]
    anc = i.get("anchor_seq", "")
    try:
        qd = int(i["designed_charge"])
        nh = int(e.get("n_hotspots_total") or SRC_HOTS.get(e.get("source", ""), 4))
        rows.append({
            "bb": (k[0], re.sub(r"__s\d+$", "", k[1])),
            "src": e.get("source", "?"),
            "q_des": qd, "q_full": qd + q(anc),
            "motif": float(d["motif_rmsd_vs_af3"]), "recept": float(d["receptor_rmsd"]),
            "hits": int(e["n_hotspots_hit"]), "n_total": nh,
            "sb": int(i.get("unique_salt_bridges", 0) or 0),
        })
    except (ValueError, KeyError):
        pass
n = len(rows)
nb = lambda s: len({r["bb"] for r in s})
print(f"{n} folds from {nb(rows)} backbones\n")


def report(label, sel):
    gate = [r for r in sel if r["motif"] < 2.0 and r["hits"] >= 3 and r["recept"] < 5.0]
    full = [r for r in gate if r["hits"] == r["n_total"]]
    print(f"{label}")
    print(f"  folds            {len(sel):4d}          backbones {nb(sel):4d}")
    print(f"  pass the gate    {len(gate):4d} "
          f"({100*len(gate)/max(len(sel),1):3.0f}%)   backbones {nb(gate):4d}")
    print(f"  gate + full epi  {len(full):4d} "
          f"({100*len(full)/max(len(sel),1):3.0f}%)   backbones {nb(full):4d}")
    if sel:
        sbv = sorted(r["sb"] for r in sel)
        print(f"  salt bridges     median {sbv[len(sbv)//2]}")
        print(f"  source mix       {dict(Counter(r['src'] for r in sel))}")
    print()


print("=" * 66)
print("LOW-CHARGE GROUP, on the controlled quantity q(full 140)")
print("=" * 66)
report("inside the target window  +4 <= q(full) <= +7",
       [r for r in rows if 4 <= r["q_full"] <= 7])
report("broader low band          +2 <= q(full) <= +9",
       [r for r in rows if 2 <= r["q_full"] <= 9])

print("=" * 66)
print("For reference: the bins used in the earlier table, on q(11-140)")
print("=" * 66)
for lo, hi, lbl in ((4, 9, "q(11-140) +4..+9   (low)"),
                    (9, 16, "q(11-140) +9..+16  (mid)"),
                    (16, 999, "q(11-140) >+16     (high tail)")):
    report(lbl, [r for r in rows if lo <= r["q_des"] < hi])

# how many of the packaged 50 are low-charge
pkg = os.path.join(os.path.dirname(A), "p14_best", "metrics.csv")
if os.path.exists(pkg):
    pr = list(csv.DictReader(open(pkg)))
    qf = [q(r["seq"]) for r in pr]
    print("=" * 66)
    print(f"Of the {len(pr)} packaged designs in p14_best:")
    for lo, hi, lbl in ((4, 7, "inside the target window +4..+7"),
                        (2, 9, "broader low band +2..+9"),
                        (10, 999, "above +9")):
        c = sum(1 for x in qf if lo <= x <= hi)
        print(f"  {lbl:34s} {c:3d}/{len(pr)}")
