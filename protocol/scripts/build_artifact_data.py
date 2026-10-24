#!/usr/bin/env python
"""Assemble one JSON of joined metrics plus the structure files for the viewer.

Joins the four measurement tables that describe the p15 set: the packaged
ranking, the unbound fold comparison, the per-window RMSD, and the
self-association run. Structures are copied for the designs the viewer can show.
"""
import csv, glob, json, os, re, shutil, sys
from collections import Counter

G, CPLX, MONO, OUT, NVIEW = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
A = os.path.join(G, "analysis")
os.makedirs(os.path.join(OUT, "struct"), exist_ok=True)

def rows(p):
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []

pkg = rows(os.path.join(G, "p15_best", "metrics.csv"))
win = {r["name"]: r for r in rows(os.path.join(A, "p16_monomer_windows.csv"))}
mon = {r["name"]: r for r in rows(os.path.join(A, "p16_monomer_vs_bound.csv"))}
selfa = {}
for r in rows(os.path.join(A, "p17_self_association.csv")):
    for nm in r.get("members", "").split(";"):
        if nm: selfa[nm] = r
print(f"packaged {len(pkg)}  windows {len(win)}  monomer {len(mon)}  self {len(selfa)}")

HOTS = {"p3": [86, 87, 88, 89], "p10": [27, 42, 86],
        "p12": [27, 42, 86], "p13": [27, 42, 86, 95]}
f = lambda v, d=None: (float(v) if v not in ("", None) else d)
i_ = lambda v, d=0: (int(float(v)) if v not in ("", None) else d)

out = []
for r in pkg:
    nm = r["name"]
    src = nm.split("_")[0]
    w, m, s = win.get(nm, {}), mon.get(nm, {}), selfa.get(nm, {})
    seq = r["seq"]
    c = Counter(seq)
    out.append({
        "tag": r["tag"], "name": nm, "shard": r["shard"], "source": src,
        "rank": i_(r["rank"]), "score": f(r["score"]),
        "hotspots": i_(r["hotspots"]), "hot_list": HOTS.get(src, []),
        "motif_rmsd": f(r["motif_rmsd"]), "receptor_rmsd": f(r["receptor_rmsd"]),
        "min_dist": f(r.get("min_dist")),
        "ipsae": f(r["ipsae"]), "ipsae_cons": f(r.get("ipsae_cons")),
        "plddt_binder": f(r["plddt_binder"]), "plddt_design": f(r.get("plddt_design")),
        "iplddt": f(r.get("iplddt")), "iptm": f(r.get("iptm")), "pdockq": f(r.get("pdockq")),
        "anchor": seq[:10], "seq": seq,
        "q_full": c["K"] + c["R"] - c["D"] - c["E"],
        "q_designed": Counter(seq[10:])["K"] + Counter(seq[10:])["R"]
                      - Counter(seq[10:])["D"] - Counter(seq[10:])["E"],
        "mono_plddt": f(m.get("mono_plddt")),
        "mono_plddt_designed": f(m.get("mono_plddt_designed")),
        "rmsd_domain": f(w.get("domain_30_140")),
        "rmsd_all": f(w.get("all_140")),
        "self_iptm": f(s.get("iptm_self")),
        "self_contacts": i_(s.get("contacts_AB")),
        "mono_id": (m.get("mono") or s.get("mono") or ""),
    })
out.sort(key=lambda r: -(r["score"] or 0))

copied = 0
for r in out:
    if copied >= NVIEW: break
    cd = os.path.join(CPLX, "sh" + r["shard"], "boltz_results_yaml",
                      "predictions", r["name"])
    cc = glob.glob(os.path.join(cd, "*_model_0.cif"))
    if not cc: continue
    base = f"{r['tag']}"
    shutil.copy(cc[0], os.path.join(OUT, "struct", base + "_complex.cif"))
    r["file_complex"] = f"struct/{base}_complex.cif"
    if r["mono_id"]:
        mc = glob.glob(os.path.join(MONO, r["mono_id"], "*_model_0.cif"))
        if mc:
            shutil.copy(mc[0], os.path.join(OUT, "struct", base + "_alone.cif"))
            r["file_alone"] = f"struct/{base}_alone.cif"
    copied += 1

agg = {
    "n_packaged": len(out), "n_viewable": copied,
    "campaign": "p15_1099477",
    "cohort": {
        "refolds": 210, "gated": 187, "full_epitope": 195,
        "consensus_compliant": 210,
        "mono_plddt_median": 0.880, "mono_ge_080": 200,
        "domain_rmsd_median": 0.42, "domain_within_3A": 203,
        "self_iptm_median": 0.262, "receptor_iptm_median": 0.923,
        "self_ge_080": 1,
    },
}
json.dump({"meta": agg, "designs": out},
          open(os.path.join(OUT, "data.json"), "w"), indent=1)
tot = sum(os.path.getsize(os.path.join(OUT, "struct", x))
          for x in os.listdir(os.path.join(OUT, "struct")))
print(f"viewable designs: {copied}   structure files: "
      f"{len(os.listdir(os.path.join(OUT,'struct')))}  ({tot/1e6:.1f} MB)")
print(f"wrote {OUT}/data.json")
