#!/usr/bin/env python
"""Compare each binder folded alone with the same chain inside its complex.

Reports the monomer's own confidence and how far it moves from the bound
conformation. A design that keeps its fold alone is a better candidate than one
that only becomes ordered against the receptor.

Superposition is on binder CA atoms only, all 140 and then the designed region
11-140 separately, since the 10-residue anchor is a short terminal segment whose
placement is driven by the receptor pocket in the complex.
"""
import csv, glob, json, os, sys
import numpy as np

MONO_DIR, CPLX_ROOT, MAN, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}


def ca_of(path, want_chain=None):
    names, out, seq = [], {}, {}
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
        seq.setdefault(ch, {})[rn] = ONE.get((d.get("label_comp_id") or "").upper(), "X")
    return out, seq


def kabsch_rmsd(P, Q):
    P, Q = np.asarray(P, float), np.asarray(Q, float)
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    if np.linalg.det(V @ Wt) < 0: V[:, -1] *= -1
    R = V @ Wt
    return float(np.sqrt((((Pc @ R) - Qc) ** 2).sum() / len(P)))


man = json.load(open(MAN))["monomers"]
rows = []
for m in man:
    md = os.path.join(MONO_DIR, m["mono"])
    mcif = glob.glob(os.path.join(md, "*_model_0.cif"))
    mconf = glob.glob(os.path.join(md, "confidence_*.json"))
    if not mcif:
        continue
    mca, _ = ca_of(mcif[0])
    mch = max(mca, key=lambda c: len(mca[c]))
    mono = mca[mch]
    plddt = None
    if mconf:
        plddt = json.load(open(mconf[0])).get("complex_plddt")
    pl_arr = glob.glob(os.path.join(md, "plddt_*.npz"))
    pl_des = pl_anchor = None
    if pl_arr:
        try:
            a = np.load(pl_arr[0])["plddt"].ravel()
            if len(a) == m["length"]:
                pl_anchor, pl_des = float(a[:10].mean()), float(a[10:].mean())
        except Exception:
            pass
    for mem in m["members"]:
        cd = os.path.join(CPLX_ROOT, "sh" + mem["shard"], "boltz_results_yaml",
                          "predictions", mem["name"])
        ccif = glob.glob(os.path.join(cd, "*_model_0.cif"))
        if not ccif: continue
        cca, _ = ca_of(ccif[0])
        bch = min(cca, key=lambda c: len(cca[c]))     # binder is the short chain
        bound = cca[bch]
        common = sorted(set(mono) & set(bound))
        if len(common) < 50: continue
        full = kabsch_rmsd([mono[i] for i in common], [bound[i] for i in common])
        des = [i for i in common if i > 10]
        dr = kabsch_rmsd([mono[i] for i in des], [bound[i] for i in des]) if len(des) > 50 else None
        rows.append({"mono": m["mono"], "shard": mem["shard"], "name": mem["name"],
                     "source": mem["name"].split("_")[0],
                     "n_sharing_seq": m["n_complex_designs"],
                     "mono_plddt": round(plddt, 4) if plddt is not None else "",
                     "mono_plddt_anchor": round(pl_anchor, 4) if pl_anchor else "",
                     "mono_plddt_designed": round(pl_des, 4) if pl_des else "",
                     "rmsd_mono_vs_bound_full": round(full, 3),
                     "rmsd_mono_vs_bound_designed": round(dr, 3) if dr else "",
                     "n_res_compared": len(common)})

if not rows:
    print("no comparisons made"); sys.exit(1)
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

import statistics as st
num = lambda k: [r[k] for r in rows if isinstance(r[k], float)]
print(f"compared {len(rows)} complex designs against {len({r['mono'] for r in rows})} monomers")
for k, lbl in (("mono_plddt", "monomer pLDDT (whole chain)"),
               ("mono_plddt_anchor", "  anchor 1-10"),
               ("mono_plddt_designed", "  designed 11-140"),
               ("rmsd_mono_vs_bound_full", "RMSD monomer vs bound, all 140"),
               ("rmsd_mono_vs_bound_designed", "  designed region only")):
    v = sorted(num(k))
    if v:
        print(f"  {lbl:34s} p10 {v[len(v)//10]:6.3f}  median {st.median(v):6.3f}  "
              f"p90 {v[9*len(v)//10]:6.3f}")
hi = [r for r in rows if isinstance(r["mono_plddt"], float) and r["mono_plddt"] >= 0.80]
lo = [r for r in rows if isinstance(r["mono_plddt"], float) and r["mono_plddt"] < 0.70]
print(f"  monomer pLDDT >= 0.80 : {len(hi)}/{len(rows)}")
print(f"  monomer pLDDT <  0.70 : {len(lo)}/{len(rows)}")
st2 = [r for r in rows if isinstance(r["rmsd_mono_vs_bound_designed"], float)
       and r["rmsd_mono_vs_bound_designed"] < 2.0]
print(f"  designed region within 2 A of bound: {len(st2)}/{len(rows)}")
print(f"wrote {OUT}")
