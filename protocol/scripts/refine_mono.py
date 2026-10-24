#!/usr/bin/env python
"""Re-measure the unbound-vs-bound difference over the compact domain only.

A median whole-chain RMSD of 10.3 A would normally mean the fold changes on
binding, but this construct is a long stalk plus a compact domain (established
earlier in the campaign, which is why segment Rg from residue 30 replaced
whole-chain Rg). Superposing across a hinge inflates RMSD even when the domain
itself is unchanged, so the domain has to be measured on its own before any
conformational-change claim is made.

Reports RMSD over several windows and checks whether the distribution is bimodal
(some designs keeping their fold, others not) rather than broadly shifted.
"""
import csv, glob, json, os, statistics as st, sys
import numpy as np

MONO_DIR, CPLX, MAN, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]


def ca(path):
    names, out = [], {}
    for l in open(path):
        if l.startswith("_atom_site."):
            names.append(l.strip().split(".")[1]); continue
        if not (l.startswith("ATOM") or l.startswith("HETATM")): continue
        t = l.split()
        if len(t) != len(names): continue
        d = dict(zip(names, t))
        if (d.get("label_atom_id") or "").strip() != "CA": continue
        ch = d.get("auth_asym_id") or d.get("label_asym_id")
        try: rn = int(d.get("auth_seq_id") or d.get("label_seq_id"))
        except (TypeError, ValueError): continue
        out.setdefault(ch, {})[rn] = (float(d["Cartn_x"]), float(d["Cartn_y"]),
                                      float(d["Cartn_z"]))
    return out


def rmsd(P, Q):
    P, Q = np.asarray(P, float), np.asarray(Q, float)
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    if np.linalg.det(V @ Wt) < 0: V[:, -1] *= -1
    return float(np.sqrt((((Pc @ (V @ Wt)) - Qc) ** 2).sum() / len(P)))


WINDOWS = [("all_140", 1, 140), ("designed_11_140", 11, 140),
           ("domain_30_140", 30, 140), ("domain_40_140", 40, 140),
           ("core_50_120", 50, 120)]
rows = []
for m in json.load(open(MAN))["monomers"]:
    mc = glob.glob(os.path.join(MONO_DIR, m["mono"], "*_model_0.cif"))
    if not mc: continue
    mm = ca(mc[0]); mono = mm[max(mm, key=lambda c: len(mm[c]))]
    for mem in m["members"]:
        cc = glob.glob(os.path.join(CPLX, "sh" + mem["shard"], "boltz_results_yaml",
                                    "predictions", mem["name"], "*_model_0.cif"))
        if not cc: continue
        cm = ca(cc[0]); bound = cm[min(cm, key=lambda c: len(cm[c]))]
        r = {"mono": m["mono"], "name": mem["name"],
             "source": mem["name"].split("_")[0]}
        for lbl, lo, hi in WINDOWS:
            idx = [i for i in range(lo, hi + 1) if i in mono and i in bound]
            r[lbl] = round(rmsd([mono[i] for i in idx], [bound[i] for i in idx]), 3) \
                if len(idx) >= 40 else ""
        rows.append(r)

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(f"compared {len(rows)} designs\n")
print(f"{'window':18s} {'p10':>7s} {'median':>8s} {'p90':>7s}   {'<2A':>8s} {'<3A':>8s}")
for lbl, _, _ in WINDOWS:
    v = sorted(r[lbl] for r in rows if isinstance(r[lbl], float))
    if not v: continue
    u2 = sum(1 for x in v if x < 2.0); u3 = sum(1 for x in v if x < 3.0)
    print(f"{lbl:18s} {v[len(v)//10]:7.2f} {st.median(v):8.2f} {v[9*len(v)//10]:7.2f}   "
          f"{u2:4d} ({100*u2/len(v):3.0f}%) {u3:4d} ({100*u3/len(v):3.0f}%)")

key = "domain_30_140"
v = sorted(r[key] for r in rows if isinstance(r[key], float))
print(f"\nis {key} bimodal? histogram of RMSD (A)")
bins = [(0, 1), (1, 2), (2, 3), (3, 5), (5, 8), (8, 12), (12, 20), (20, 999)]
for lo, hi in bins:
    c = sum(1 for x in v if lo <= x < hi)
    bar = "#" * int(60 * c / max(len(v), 1))
    print(f"  {lo:3d}-{hi if hi<999 else '+':>3} {c:4d}  {bar}")
keep = [r for r in rows if isinstance(r[key], float) and r[key] < 3.0]
print(f"\ndomain within 3 A of its bound conformation: {len(keep)}/{len(rows)} "
      f"({100*len(keep)/len(rows):.0f}%)")
print("by source:")
for s in ("p3", "p10", "p12", "p13"):
    sub = [r for r in rows if r["source"] == s and isinstance(r[key], float)]
    if not sub: continue
    k = sum(1 for r in sub if r[key] < 3.0)
    print(f"  {s:4s} {k:3d}/{len(sub):3d} ({100*k/len(sub):3.0f}%)   "
          f"median {st.median([r[key] for r in sub]):6.2f} A")
print(f"\nwrote {OUT}")
