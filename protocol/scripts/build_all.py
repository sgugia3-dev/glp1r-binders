#!/usr/bin/env python
"""Collect all three predictions for every packaged design, as PDB.

PDB rather than mmCIF because 3Dmol's mmCIF reader silently truncated these
files in the browser. The B-factor column carries per-residue pLDDT on a 0-100
scale, which is what the viewer's pLDDT colouring reads.

The single-chain and two-chain runs numbered their inputs independently, so the
two-chain structure is matched by binder SEQUENCE, never by index.
"""
import csv, glob, json, os, sys
from collections import Counter, defaultdict

G, CPLX, MONOD, MONOMAN, DIMD, DIMMAN, OUT = sys.argv[1:8]
A = os.path.join(G, "analysis")
os.makedirs(os.path.join(OUT, "pdb"), exist_ok=True)


def read_cif(p):
    names, rows = [], []
    for l in open(p):
        if l.startswith("_atom_site."): names.append(l.strip().split(".")[1]); continue
        if not l.startswith(("ATOM", "HETATM")): continue
        t = l.split()
        if len(t) == len(names): rows.append(dict(zip(names, t)))
    return rows


def to_pdb(rows, out, biggest_is):
    byc = defaultdict(set)
    for r in rows:
        c = r.get("auth_asym_id") or r.get("label_asym_id")
        byc[c].add(int(r.get("auth_seq_id") or r.get("label_seq_id")))
    order = sorted(byc, key=lambda c: -len(byc[c]))
    if biggest_is == "receptor" and len(order) >= 2:
        newid = {order[0]: "B", order[1]: "A"}
    else:
        newid = {c: chr(ord("A") + i) for i, c in enumerate(sorted(byc))}
    n, lines = 0, []
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
        lines.append("ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
                     % (n, nm, (r.get("label_comp_id") or "UNK").upper()[:3],
                        newid[c], rn, x, y, z, 1.0, b, el[:2].rjust(2)))
    lines.append("END\n")
    open(out, "w").writelines(lines)
    return {newid[c]: len(v) for c, v in byc.items()}


pkg = list(csv.DictReader(open(os.path.join(G, "p15_best", "metrics.csv"))))
mono = {m["seq"]: m for m in json.load(open(MONOMAN))["monomers"]}
dim = {}
for e in json.load(open(DIMMAN))["monomers"]:
    dim.setdefault(e["seq"], e["mono"])
win = {r["name"]: r for r in csv.DictReader(open(os.path.join(A, "p16_monomer_windows.csv")))}
mvb = {r["name"]: r for r in csv.DictReader(open(os.path.join(A, "p16_monomer_vs_bound.csv")))}
sf = {}
for r in csv.DictReader(open(os.path.join(A, "p17_self_association.csv"))):
    for nm in r.get("members", "").split(";"):
        if nm: sf[nm] = r
HOTS = {"p3": [86, 87, 88, 89], "p10": [27, 42, 86],
        "p12": [27, 42, 86], "p13": [27, 42, 86, 95]}
f = lambda v, d=None: (float(v) if v not in ("", None) else d)
i_ = lambda v, d=0: (int(float(v)) if v not in ("", None) else d)

out, bad = [], []
for r in pkg:
    nm, tag = r["name"], r["tag"]
    src = nm.split("_")[0]
    seq = r["seq"]
    c = Counter(seq); cd = Counter(seq[10:])
    w, m, s = win.get(nm, {}), mvb.get(nm, {}), sf.get(nm, {})
    rec = {"tag": tag, "name": nm, "shard": r["shard"], "source": src,
           "rank": i_(r["rank"]), "score": f(r["score"]),
           "hotspots": i_(r["hotspots"]), "hot_list": HOTS.get(src, []),
           "motif_rmsd": f(r["motif_rmsd"]), "receptor_rmsd": f(r["receptor_rmsd"]),
           "min_dist": f(r.get("min_dist")), "ipsae": f(r["ipsae"]),
           "ipsae_cons": f(r.get("ipsae_cons")), "plddt_binder": f(r["plddt_binder"]),
           "plddt_design": f(r.get("plddt_design")), "iplddt": f(r.get("iplddt")),
           "iptm": f(r.get("iptm")), "pdockq": f(r.get("pdockq")),
           "anchor": seq[:10], "seq": seq,
           "q_full": c["K"] + c["R"] - c["D"] - c["E"],
           "q_designed": cd["K"] + cd["R"] - cd["D"] - cd["E"],
           "mono_plddt": f(m.get("mono_plddt")),
           "mono_plddt_designed": f(m.get("mono_plddt_designed")),
           "rmsd_domain": f(w.get("domain_30_140")), "rmsd_all": f(w.get("all_140")),
           "self_iptm": f(s.get("iptm_self")), "self_contacts": i_(s.get("contacts_AB"))}

    cc = glob.glob(os.path.join(CPLX, "sh" + r["shard"], "boltz_results_yaml",
                                "predictions", nm, "*_model_0.cif"))
    if cc:
        cnt = to_pdb(read_cif(cc[0]), os.path.join(OUT, "pdb", tag + "_complex.pdb"),
                     "receptor")
        if cnt.get("A") == 140 and cnt.get("B") == 409:
            rec["file_complex"] = f"pdb/{tag}_complex.pdb"
        else: bad.append((tag, "complex", cnt))
    mo = mono.get(seq)
    if mo:
        mc = glob.glob(os.path.join(MONOD, mo["mono"], "*_model_0.cif"))
        if mc:
            cnt = to_pdb(read_cif(mc[0]), os.path.join(OUT, "pdb", tag + "_alone.pdb"), "self")
            if cnt.get("A") == 140: rec["file_alone"] = f"pdb/{tag}_alone.pdb"
            else: bad.append((tag, "alone", cnt))
    di = dim.get(seq)
    if di:
        dc = glob.glob(os.path.join(DIMD, di, "*_model_0.cif"))
        if dc:
            cnt = to_pdb(read_cif(dc[0]), os.path.join(OUT, "pdb", tag + "_pair.pdb"), "self")
            if cnt.get("A") == 140 and cnt.get("B") == 140:
                rec["file_pair"] = f"pdb/{tag}_pair.pdb"
            else: bad.append((tag, "pair", cnt))
    out.append(rec)

out.sort(key=lambda r: -(r["score"] or 0))
nc = sum(1 for r in out if r.get("file_complex"))
na = sum(1 for r in out if r.get("file_alone"))
np_ = sum(1 for r in out if r.get("file_pair"))
meta = {"n_packaged": len(out), "n_viewable": nc, "campaign": "p15_1099477",
        "cohort": {"refolds": 210, "gated": 187, "full_epitope": 195,
                   "consensus_compliant": 210, "mono_plddt_median": 0.880,
                   "mono_ge_080": 200, "domain_rmsd_median": 0.42,
                   "domain_within_3A": 203, "self_iptm_median": 0.262,
                   "receptor_iptm_median": 0.923, "self_ge_080": 1}}
json.dump({"meta": meta, "designs": out}, open(os.path.join(OUT, "data.json"), "w"), indent=1)
tot = sum(os.path.getsize(os.path.join(OUT, "pdb", x))
          for x in os.listdir(os.path.join(OUT, "pdb")))
print(f"designs {len(out)}   complex {nc}   alone {na}   pair {np_}")
print(f"files {len(os.listdir(os.path.join(OUT,'pdb')))}   total {tot/1e6:.1f} MB")
print("chain-count problems:", bad or "none")
