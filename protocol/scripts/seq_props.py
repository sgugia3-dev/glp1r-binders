#!/usr/bin/env python
"""Physicochemical analysis of every GLP-1R campaign design and of the receptor.

Sequences come from the Boltz input YAMLs (chain A = binder, chain B = receptor),
so every refold in every campaign is covered, not just the packaged best sets.

Residue classes partition all 20 amino acids exactly once, stated explicitly so
the fractions are unambiguous:
  hydrophobic  A V L I M F W C   (apolar side chains, aromatics F/W and C included)
  polar        S T N Q Y H       (polar uncharged at pH 7; Y and H are borderline)
  positive     K R
  negative     D E
  special      G P               (conformationally special, poor fit to either class)
GRAVY (Kyte-Doolittle) is reported alongside as the continuous hydropathy measure.

The binder is 140 aa: a FIXED 10-residue anchor (GLP-1 7-16, HAEGTFTSDV, identical
in every design) plus 130 designed residues. Properties are reported for the
designed region as well as the whole chain, because the constant anchor dilutes
real differences between campaigns.
"""
import csv, glob, json, os, re, statistics as st, sys
from collections import Counter
from Bio.SeqUtils.ProtParam import ProteinAnalysis

SCRATCH, G, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
CAMPAIGNS = [("p2_1088883", "P2  A24/25/27/42"), ("p3_1088996", "P3  A86-89"),
             ("p10_1089071", "P10 A27/42/86 smoke"), ("p12_1089085", "P12 A27/42/86"),
             ("p13_1089119", "P13 A27/42/86/95")]
CLASSES = {"hydrophobic": set("AVLIMFWC"), "polar": set("STNQYH"),
           "positive": set("KR"), "negative": set("DE"), "special": set("GP")}
assert set().union(*CLASSES.values()) == set("ACDEFGHIKLMNPQRSTVWY")
assert sum(len(v) for v in CLASSES.values()) == 20
AROMATIC = set("FWY")
PH_CURVE = [4.0, 5.0, 6.0, 7.0, 7.4, 8.0, 9.0]


def props(seq):
    p = ProteinAnalysis(seq)
    n = len(seq)
    c = Counter(seq)
    d = {"len": n, "mw_kda": p.molecular_weight() / 1000.0,
         "pI": p.isoelectric_point(), "gravy": p.gravy(),
         "aromaticity": p.aromaticity(), "instability": p.instability_index()}
    for ph in PH_CURVE:
        d[f"q_pH{ph:g}"] = p.charge_at_pH(ph)
    for k, s in CLASSES.items():
        d[f"f_{k}"] = sum(c[a] for a in s) / n
    d["f_aromatic"] = sum(c[a] for a in AROMATIC) / n
    d["f_charged"] = d["f_positive"] + d["f_negative"]          # FCR
    d["ncpr"] = (sum(c[a] for a in "KR") - sum(c[a] for a in "DE")) / n
    # Ikai aliphatic index: thermostability proxy
    pc = {a: 100.0 * c[a] / n for a in "AVIL"}
    d["aliphatic_idx"] = pc["A"] + 2.9 * pc["V"] + 3.9 * (pc["I"] + pc["L"])
    h, t, s2 = p.secondary_structure_fraction()
    d["ss_helix"], d["ss_turn"], d["ss_sheet"] = h, t, s2
    d["n_cys"] = c["C"]
    return d


FIELDS = ["len", "mw_kda", "pI"] + [f"q_pH{p:g}" for p in PH_CURVE] + [
    "ncpr", "f_charged", "f_hydrophobic", "f_polar", "f_positive", "f_negative",
    "f_special", "f_aromatic", "gravy", "aliphatic_idx", "aromaticity",
    "instability", "ss_helix", "ss_turn", "ss_sheet", "n_cys"]

# ---- collect every design sequence from the YAMLs -------------------------
rows, receptor = [], None
for run, lbl in CAMPAIGNS:
    ys = sorted(glob.glob(os.path.join(SCRATCH, "runs/glp1r", run, "sh*", "yaml", "*.yaml")))
    seen = set()
    for y in ys:
        sh = re.search(r"/sh(\d+)/", y)
        txt = open(y).read()
        seqs = re.findall(r"id:\s*([AB])\s*\n\s*sequence:\s*([A-Z]+)", txt)
        dd = dict(seqs)
        b = dd.get("A")
        if receptor is None and dd.get("B"):
            receptor = dd["B"]
        if not b: continue
        name = os.path.basename(y)[:-5]
        k = (run, sh.group(1) if sh else "?", name)
        if k in seen: continue
        seen.add(k)
        rows.append({"campaign": run, "label": lbl, "shard": k[1], "name": name,
                     "seq": b})
    print(f"{lbl:22s} {len(seen):5d} sequences from {len(ys)} yamls")
print(f"receptor length {len(receptor)}\n")

# unique designed sequences per campaign (MPNN makes several per backbone)
for run, lbl in CAMPAIGNS:
    sub = [r for r in rows if r["campaign"] == run]
    u = len({r["seq"] for r in sub})
    print(f"  {lbl:22s} {len(sub):5d} folds, {u:5d} unique binder sequences")

# ---- per-design properties ------------------------------------------------
for r in rows:
    r["full"] = props(r["seq"])
    r["des"] = props(r["seq"][10:])

with open(os.path.join(OUT, "campaign_seq_properties.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["campaign", "label", "shard", "name", "region"] + FIELDS)
    for r in rows:
        for reg, key in (("full_140", "full"), ("designed_130", "des")):
            w.writerow([r["campaign"], r["label"], r["shard"], r["name"], reg]
                       + [round(r[key][f], 4) for f in FIELDS])
print(f"\nwrote campaign_seq_properties.csv ({2*len(rows)} rows)")


def summary(sel, key):
    out = {}
    for f in FIELDS:
        v = [x[key][f] for x in sel]
        out[f] = (st.median(v), st.mean(v), st.pstdev(v) if len(v) > 1 else 0.0,
                  min(v), max(v))
    return out


# ---- receptor and its subregions -----------------------------------------
HOT = {"P2 epitope 24/25/27/42": [24, 25, 27, 42],
       "P3 epitope 86-89": [86, 87, 88, 89],
       "P12 epitope 27/42/86": [27, 42, 86],
       "P13 epitope 27/42/86/95": [27, 42, 86, 95]}
REGIONS = [("receptor FULL 1-409", receptor),
           ("receptor ECD 1-117", receptor[:117]),
           ("receptor TMD 118-409", receptor[117:])]
for nm, hs in HOT.items():
    idx = sorted({i for h in hs for i in range(max(0, h - 6), min(len(receptor), h + 5))})
    REGIONS.append((f"{nm} +-5", "".join(receptor[i] for i in idx)))

print("\n" + "=" * 104)
print("RECEPTOR AND EPITOPE REGIONS")
print("=" * 104)
hdr = f"{'region':28s} {'len':>4s} {'pI':>6s} {'q7.4':>7s} {'NCPR':>7s} " \
      f"{'FCR':>6s} {'hydr':>6s} {'polar':>6s} {'pos':>6s} {'neg':>6s} {'GRAVY':>7s}"
print(hdr)
recrows = []
for nm, s in REGIONS:
    d = props(s)
    recrows.append((nm, s, d))
    print(f"{nm:28s} {d['len']:4d} {d['pI']:6.2f} {d['q_pH7.4']:7.2f} "
          f"{d['ncpr']:7.3f} {d['f_charged']:6.3f} {d['f_hydrophobic']:6.3f} "
          f"{d['f_polar']:6.3f} {d['f_positive']:6.3f} {d['f_negative']:6.3f} "
          f"{d['gravy']:7.3f}")
with open(os.path.join(OUT, "receptor_properties.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["region", "sequence"] + FIELDS)
    for nm, s, d in recrows:
        w.writerow([nm, s] + [round(d[f], 4) for f in FIELDS])

# ---- per-campaign tables -------------------------------------------------
for reg_lbl, key in (("DESIGNED REGION (130 aa, the part that varies)", "des"),
                     ("FULL BINDER CHAIN (140 aa, includes the constant anchor)", "full")):
    print("\n" + "=" * 104)
    print(reg_lbl)
    print("=" * 104)
    print(f"{'campaign':22s} {'n':>5s} {'pI':>13s} {'q@7.4':>14s} {'NCPR':>8s} "
          f"{'FCR':>7s} {'hydrophob':>10s} {'polar':>8s} {'GRAVY':>8s}")
    for run, lbl in CAMPAIGNS:
        sel = [r for r in rows if r["campaign"] == run]
        if not sel: continue
        s = summary(sel, key)
        print(f"{lbl:22s} {len(sel):5d} "
              f"{s['pI'][0]:5.2f}+-{s['pI'][2]:<4.2f}  "
              f"{s['q_pH7.4'][0]:6.2f}+-{s['q_pH7.4'][2]:<4.2f} "
              f"{s['ncpr'][0]:8.3f} {s['f_charged'][0]:7.3f} "
              f"{s['f_hydrophobic'][0]:10.3f} {s['f_polar'][0]:8.3f} {s['gravy'][0]:8.3f}")

print("\n" + "=" * 104)
print("FULL METRIC TABLE, designed region, median (min - max)")
print("=" * 104)
print(f"{'metric':16s} " + "".join(f"{l.split()[0]:>17s}" for _, l in CAMPAIGNS))
for f in FIELDS:
    line = f"{f:16s} "
    for run, lbl in CAMPAIGNS:
        sel = [r for r in rows if r["campaign"] == run]
        if not sel: line += f"{'-':>17s}"; continue
        v = [x["des"][f] for x in sel]
        line += f"{st.median(v):>8.3f}{'':1s}({min(v):.2f}-{max(v):.2f})".rjust(17)
    print(line)

# ---- composition, pooled -------------------------------------------------
print("\n" + "=" * 104)
print("AMINO-ACID COMPOSITION of the designed region, % (receptor for reference)")
print("=" * 104)
rec_c = Counter(receptor); rec_n = len(receptor)
print(f"{'aa':>3s} {'class':12s}" + "".join(f"{l.split()[0]:>9s}" for _, l in CAMPAIGNS)
      + f"{'RECEPTOR':>10s}")
cls_of = {a: k for k, v in CLASSES.items() for a in v}
for a in "ACDEFGHIKLMNPQRSTVWY":
    line = f"{a:>3s} {cls_of[a]:12s}"
    for run, lbl in CAMPAIGNS:
        sel = [r for r in rows if r["campaign"] == run]
        if not sel: line += f"{'-':>9s}"; continue
        tot = sum(len(r["seq"]) - 10 for r in sel)
        cnt = sum(Counter(r["seq"][10:])[a] for r in sel)
        line += f"{100.0*cnt/tot:9.2f}"
    line += f"{100.0*rec_c[a]/rec_n:10.2f}"
    print(line)
json.dump({"receptor": {nm: {f: d[f] for f in FIELDS} for nm, s, d in recrows}},
          open(os.path.join(OUT, "receptor_properties.json"), "w"), indent=2)
print(f"\nwrote receptor_properties.csv / .json")
