#!/usr/bin/env python
"""Consensus candidates + sequence diversity for the P2 ranking.

Two problems with the flat mean of min-max normalised terms:
  1. The three metrics have very different dynamic ranges on the gated set
     (motif RMSD spans 7x; ipSAE spans 0.13 absolute). Min-max rescaling hands
     them equal weight, so ipSAE noise counts as much as real RMSD signal.
  2. ipSAE barely tracks the epitope (4 hits 0.882 vs 0 hits 0.871), so that
     equal weight is being spent on a metric with almost no on-target signal.
So alongside the ranking, take the weighting-free CONSENSUS set: designs in the
top quartile of all three metrics independently, which no reweighting can change.
Then cluster by sequence, because N top rows is not N candidates if they are
all the same backbone.
"""
import csv, glob, os, sys, statistics as st

ranked, root, out = sys.argv[1], sys.argv[2], sys.argv[3]
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}

rows = list(csv.DictReader(open(ranked)))
for r in rows:
    for k in ("score", "motif_rmsd", "ipsae", "plddt_binder", "receptor_rmsd",
              "ipsae_cons", "plddt_design", "iptm", "pdockq", "min_dist"):
        r[k] = float(r[k])
    r["hotspots"] = int(r["hotspots"])
print(f"gated designs: {len(rows)}")


def q(vals, p):
    s = sorted(vals)
    return s[min(len(s) - 1, int(p * len(s)))]


# top quartile on each metric, independently
t_rmsd = q([r["motif_rmsd"] for r in rows], 0.25)          # lower is better
t_ips = q([r["ipsae"] for r in rows], 0.75)
t_pl = q([r["plddt_binder"] for r in rows], 0.75)
print(f"top-quartile thresholds: motif_rmsd <= {t_rmsd:.3f}  "
      f"ipsae >= {t_ips:.3f}  plddt_binder >= {t_pl:.3f}")

cons = [r for r in rows if r["motif_rmsd"] <= t_rmsd and r["ipsae"] >= t_ips
        and r["plddt_binder"] >= t_pl]
cons4 = [r for r in cons if r["hotspots"] == 4]
print(f"consensus (top quartile on all three): {len(cons)}")
print(f"  of those with 4/4 hotspots:          {len(cons4)}")
exp = len(rows) * 0.25 ** 3
print(f"  expected if the three were independent: {exp:.1f}  "
      f"-> enrichment {len(cons)/exp:.1f}x")


def seq_of(d):
    cif = glob.glob(os.path.join(d, "*_model_0.cif"))
    if not cif: return None
    names, seen, out_s, first = [], set(), [], None
    for line in open(cif[0]):
        if line.startswith("_atom_site."):
            names.append(line.strip().split(".")[1]); continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")): continue
        t = line.split()
        if len(t) != len(names): continue
        dd = dict(zip(names, t))
        if dd.get("label_atom_id") != "CA": continue
        ch = dd.get("auth_asym_id") or dd.get("label_asym_id")
        if first is None: first = ch
        if ch != first: break
        rn = dd.get("auth_seq_id") or dd.get("label_seq_id")
        if (ch, rn) in seen: continue
        seen.add((ch, rn))
        out_s.append(AA3.get(dd.get("label_comp_id", "").upper(), "X"))
    return "".join(out_s)


# sequences for the consensus set plus the top 40 by score
want = {(r["shard"], r["name"]): r for r in cons}
for r in rows[:40]: want.setdefault((r["shard"], r["name"]), r)
print(f"\nreading sequences for {len(want)} designs")
got = 0
for (sh, nm), r in want.items():
    d = os.path.join(root, f"sh{sh}", "boltz_results_yaml", "predictions", nm)
    s = seq_of(d) if os.path.isdir(d) else None
    r["seq"] = s or ""
    if s: got += 1
print(f"  resolved {got}/{len(want)}")


def ident(a, b):
    n = min(len(a), len(b))
    if not n: return 0.0
    return sum(x == y for x, y in zip(a[:n], b[:n])) / n


# greedy cluster the consensus set at 70% identity over the designed region
pool = sorted([r for r in cons if r["seq"]], key=lambda r: -r["score"])
clusters = []
for r in pool:
    des = r["seq"][10:]
    for c in clusters:
        if ident(des, c[0]["seq"][10:]) >= 0.70:
            c.append(r); break
    else:
        clusters.append([r])
print(f"\nconsensus set clusters at 70% identity over the designed region: {len(clusters)}")
print(f"  cluster sizes: {[len(c) for c in clusters][:25]}")

reps = [c[0] for c in clusters]
cols = ["rank", "cluster_size", "score", "shard", "name", "binder_len", "motif_rmsd",
        "receptor_rmsd", "hotspots", "min_dist", "ipsae", "ipsae_cons",
        "plddt_binder", "plddt_design", "iptm", "pdockq", "seq"]
with open(out, "w") as fh:
    w = csv.writer(fh); w.writerow(cols)
    for i, (c, r) in enumerate(zip(clusters, reps), 1):
        w.writerow([i, len(c)] + [r[k] for k in cols[2:]])

print(f"\n{'#':>3s} {'siz':>3s} {'score':>6s} {'design':30s} {'sh':>3s} {'motif':>6s} "
      f"{'hit':>4s} {'ipSAE':>6s} {'pLDDTb':>7s} {'L':>4s}")
for i, (c, r) in enumerate(zip(clusters, reps), 1):
    print(f"{i:3d} {len(c):3d} {r['score']:6.3f} "
          f"{r['name'].replace('prod_r3_prod_r3_',''):30s} {r['shard']:>3s} "
          f"{r['motif_rmsd']:6.2f} {r['hotspots']:4d} {r['ipsae']:6.3f} "
          f"{r['plddt_binder']:7.3f} {len(r['seq']):4d}")

print("\ntop cluster representatives, designed region (motif = first 10):")
for i, r in enumerate(reps[:12], 1):
    print(f"  {i:2d} {r['name'].replace('prod_r3_prod_r3_','')} sh{r['shard']}")
    print(f"     {r['seq'][:10]} | {r['seq'][10:]}")
print(f"\nwrote {out}")
