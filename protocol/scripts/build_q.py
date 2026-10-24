#!/usr/bin/env python
"""Bundle the three replicate structures for a given list of designs."""
import csv, glob, json, os, sys
from collections import defaultdict

NAMES, A, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
RUNS = {"s1": "p15_1099477", "s2": "p18s2_1106925", "s3": "p18s3_1106926"}
BASE = "/lustre/scratch/gsgueglia/challenge2026/runs/glp1r"
HOT = {"p3": [86, 87, 88, 89], "p10": [27, 42, 86],
       "p12": [27, 42, 86], "p13": [27, 42, 86, 95]}
os.makedirs(OUT, exist_ok=True)


def read_cif(p):
    names, rows = [], []
    for l in open(p):
        if l.startswith("_atom_site."): names.append(l.strip().split(".")[1]); continue
        if not l.startswith(("ATOM", "HETATM")): continue
        t = l.split()
        if len(t) == len(names): rows.append(dict(zip(names, t)))
    return rows


def to_pdb(rows):
    byc = defaultdict(set)
    for r in rows:
        c = r.get("auth_asym_id") or r.get("label_asym_id")
        byc[c].add(int(r.get("auth_seq_id") or r.get("label_seq_id")))
    order = sorted(byc, key=lambda c: -len(byc[c]))
    newid = {order[0]: "B", order[1]: "A"} if len(order) >= 2 else {order[0]: "A"}
    n, out = 0, []
    for r in rows:
        el = (r.get("type_symbol") or "").strip()
        if el == "H": continue
        c = r.get("auth_asym_id") or r.get("label_asym_id")
        an = (r.get("label_atom_id") or "").strip().strip('"')
        nm = an.ljust(4) if (len(an) >= 4 or len(el) > 1) else " " + an.ljust(3)
        try:
            rn = int(r.get("auth_seq_id") or r.get("label_seq_id"))
            x, y, z = float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])
            b = float(r.get("B_iso_or_equiv", 0) or 0)
        except Exception:
            continue
        n += 1
        out.append("ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
                   % (n, nm, (r.get("label_comp_id") or "UNK").upper()[:3],
                      newid[c], rn, x, y, z, 1.0, b, el[:2].rjust(2)))
    out.append("END\n")
    return "".join(out), {newid[c]: len(v) for c, v in byc.items()}


cons = {r["name"]: r for r in csv.DictReader(open(os.path.join(A, "seed_consistency.csv")))}
seqs = {}
for m in sorted(glob.glob(os.path.join(BASE, "p15_1099477", "sh*", "manifest.json"))):
    for y in json.load(open(m))["yamls"]:
        seqs[y["name"]] = y["binder_seq"]

order = [l.split()[1] for l in open(NAMES) if l.strip() and l.split()[0].isdigit()]
bundles = {k: {} for k in RUNS}
designs, bad = [], []
for i, name in enumerate(order, 1):
    src = name.split("_")[0]
    tag = f"q{i:02d}"
    rec = {"tag": tag, "name": name, "parent": "_".join(name.split("_")[:2]),
           "source": src, "rank": i, "hot_list": HOT.get(src, []),
           "seq": seqs.get(name, ""), "anchor": seqs.get(name, "")[:10]}
    c = cons.get(name)
    for k, pref in (("s1", "p15"), ("s2", "p18s2"), ("s3", "p18s3")):
        cif = glob.glob(os.path.join(BASE, RUNS[k], "sh*", "boltz_results_yaml",
                                     "predictions", name, "*_model_0.cif"))
        if not cif:
            bad.append((tag, k, "no structure")); continue
        txt, counts = to_pdb(read_cif(cif[0]))
        if counts.get("A") != 140 or counts.get("B") != 409:
            bad.append((tag, k, counts)); continue
        bundles[k][tag] = txt
        if c:
            rec[k] = {m: (float(c[f"{pref}_{m}"]) if m not in ("hits", "sb", "hc")
                          else int(float(c[f"{pref}_{m}"])))
                      for m in ("motif", "recept", "hits", "sb", "hc", "ipsae",
                                "plddt", "iptm")}
    if c:
        rec["n_seeds_gate"] = int(c["n_seeds_gate"])
        rec["n_seeds_full"] = int(c["n_seeds_full"])
        rec["n_total"] = int(c["n_total"])
        v = [rec[k]["ipsae"] for k in RUNS if k in rec]
        mo = [rec[k]["motif"] for k in RUNS if k in rec]
        pl = [rec[k]["plddt"] for k in RUNS if k in rec]
        rec["mean_ipsae"] = round(sum(v) / len(v), 4)
        rec["range_ipsae"] = round(max(v) - min(v), 4)
        rec["mean_motif"] = round(sum(mo) / len(mo), 3)
        rec["range_motif"] = round(max(mo) - min(mo), 3)
        rec["mean_plddt"] = round(sum(pl) / len(pl), 4)
    for k in RUNS: rec[f"b_{k}"] = f"q_{k}.json"
    designs.append(rec)

for k, b in bundles.items():
    p = os.path.join(OUT, f"q_{k}.json")
    json.dump(b, open(p, "w"))
    print(f"  q_{k}.json: {len(b)} structures {os.path.getsize(p)/1e6:.2f} MB")
import statistics as st
ok = [d for d in designs if "mean_ipsae" in d]
meta = {"n_designs": len(designs), "replicates": list(RUNS.values()),
        "n_gate3": sum(1 for d in ok if d["n_seeds_gate"] == 3),
        "n_full3": sum(1 for d in ok if d["n_seeds_full"] == 3),
        "median_mean_ipsae": round(st.median([d["mean_ipsae"] for d in ok]), 3),
        "median_mean_plddt": round(st.median([d["mean_plddt"] for d in ok]), 3),
        "median_mean_motif": round(st.median([d["mean_motif"] for d in ok]), 2),
        "median_range_motif": round(st.median([d["range_motif"] for d in ok]), 2),
        "median_range_ipsae": round(st.median([d["range_ipsae"] for d in ok]), 3)}
json.dump({"meta": meta, "designs": designs}, open(os.path.join(OUT, "data.json"), "w"), indent=1)
print(f"designs {len(designs)}  with metrics {len(ok)}  problems {bad or 'none'}")
print(" ", {k: v for k, v in meta.items() if k != "replicates"})
