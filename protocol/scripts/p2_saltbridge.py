#!/usr/bin/env python
"""Salt bridges between receptor Arg167 and binder acidic side chains.

A salt bridge is counted as a guanidinium nitrogen of the target Arg (NE, NH1, NH2) within
CUT of a carboxylate oxygen of a binder Asp (OD1/OD2) or Glu (OE1/OE2). Side-chain atoms
only: a backbone contact is not a salt bridge.

Reported separately for acids in the FIXED MOTIF (binder 1-10, which carries Glu3 and Asp9
natively) and acids in the DESIGNED region, because only the designed ones are a result.
Residue identity of the target position is asserted, not assumed.
"""
import argparse, csv, glob, json, math, os, re, sys

ARG_N = {"NE", "NH1", "NH2"}
ACID_O = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}


def read(path):
    hdr, rows, inl = [], [], False
    for line in open(path):
        if line.startswith("_atom_site."):
            hdr.append(line.strip().split(".")[1]); inl = True; continue
        if inl and line.startswith(("ATOM", "HETATM")): rows.append(line.split())
        elif inl and line.startswith("#"): break
    i = {h: k for k, h in enumerate(hdr)}
    need = ("auth_asym_id", "auth_seq_id", "label_comp_id", "label_atom_id",
            "Cartn_x", "Cartn_y", "Cartn_z")
    if any(n not in i for n in need): return None
    out = {}
    for r in rows:
        if len(r) <= max(i[n] for n in need): continue
        try: n = int(r[i["auth_seq_id"]])
        except ValueError: continue
        out.setdefault(r[i["auth_asym_id"]], {}).setdefault(n, []).append(
            (r[i["label_comp_id"]], r[i["label_atom_id"]],
             (float(r[i["Cartn_x"]]), float(r[i["Cartn_y"]]), float(r[i["Cartn_z"]]))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-root", required=True, help="dir holding sh*/ shard dirs")
    ap.add_argument("--analysis", required=True)
    ap.add_argument("--target-resid", type=int, default=167)
    ap.add_argument("--cut", type=float, default=4.0)
    ap.add_argument("--out-csv", required=True)
    a = ap.parse_args()

    # candidate set from the aggregate, keyed (shard, name)
    cand = set()
    summ = os.path.join(a.analysis, "p1_summary.json")
    if os.path.exists(summ):
        for c in json.load(open(summ))["candidates"]:
            cand.add((str(c["shard"]), c["name"]))
    print(f"candidate set from p1_summary.json: {len(cand)}")

    rows, wrong_id, n = [], 0, 0
    for shdir in sorted(glob.glob(os.path.join(a.run_root, "sh*"))):
        sh = re.sub(r".*sh", "", os.path.basename(shdir))
        man_p = os.path.join(shdir, "manifest.json")
        if not os.path.exists(man_p): continue
        man = {y["name"]: y for y in json.load(open(man_p))["yamls"]}
        for cif in glob.glob(os.path.join(shdir, "**", "predictions", "*", "*_model_0.cif"),
                             recursive=True):
            name = os.path.basename(os.path.dirname(cif))
            rec = man.get(name)
            if not rec: continue
            ch = read(cif)
            if not ch or len(ch) < 2: continue
            tgt = max(ch, key=lambda c: len(ch[c]))
            bnd = [c for c in ch if c != tgt][0]
            at = ch[tgt].get(a.target_resid)
            if not at: continue
            if at[0][0] != "ARG":
                wrong_id += 1; continue
            gn = [p for (cn, an, p) in at if an in ARG_N]
            if not gn: continue
            motif = set(rec["motif"])
            best_m, best_d, nb_m, nb_d = 99.0, 99.0, 0, 0
            for pos, atoms in ch[bnd].items():
                cn = atoms[0][0]
                if cn not in ACID_O: continue
                os_ = [p for (c2, an, p) in atoms if an in ACID_O[cn]]
                if not os_: continue
                dmin = min(math.dist(p, q) for p in os_ for q in gn)
                if pos in motif:
                    best_m = min(best_m, dmin)
                    if dmin <= a.cut: nb_m += 1
                else:
                    best_d = min(best_d, dmin)
                    if dmin <= a.cut: nb_d += 1
            n += 1
            rows.append({"shard": sh, "name": name, "binder_len": rec["binder_len"],
                         "is_candidate": (sh, name) in cand,
                         "n_bridges_designed": nb_d, "n_bridges_motif": nb_m,
                         "min_dist_designed": round(best_d, 2),
                         "min_dist_motif": round(best_m, 2),
                         "any_bridge": (nb_d + nb_m) > 0})

    if not rows: print("FATAL: nothing measured"); return 3
    cols = list(rows[0].keys())
    with open(a.out_csv, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows: fh.write(",".join(str(r[c]) for c in cols) + "\n")

    import statistics as st
    print(f"\nscored {n} refolds (target ARG{a.target_resid}, cut {a.cut} A), "
          f"wrong residue identity: {wrong_id}")
    for label, sel in (("ALL scored", rows),
                       ("TOP candidates only", [r for r in rows if r["is_candidate"]])):
        if not sel: continue
        bd = sum(1 for r in sel if r["n_bridges_designed"] > 0)
        bm = sum(1 for r in sel if r["n_bridges_motif"] > 0)
        ba = sum(1 for r in sel if r["any_bridge"])
        dd = [r["min_dist_designed"] for r in sel if r["min_dist_designed"] < 99]
        print(f"\n--- {label}: n={len(sel)} ---")
        print(f"  salt bridge from a DESIGNED acid : {bd} ({100*bd/len(sel):.1f}%)")
        print(f"  salt bridge from a MOTIF acid    : {bm} ({100*bm/len(sel):.1f}%)")
        print(f"  any salt bridge                  : {ba} ({100*ba/len(sel):.1f}%)")
        if dd: print(f"  nearest designed acid: min {min(dd):.2f} med {st.median(dd):.2f} A")
    top = [r for r in rows if r["is_candidate"] and r["n_bridges_designed"] > 0]
    top.sort(key=lambda r: r["min_dist_designed"])
    if top:
        print(f"\ntop candidates WITH a designed salt bridge to ARG{a.target_resid}: {len(top)}")
        print("%-40s %3s %4s %7s %6s" % ("design", "sh", "len", "minDist", "nBr"))
        for r in top[:15]:
            print("%-40s %3s %4d %7.2f %6d" % (r["name"].replace("glp1r_ecd4_glp1r_", "")[:40],
                  r["shard"], r["binder_len"], r["min_dist_designed"], r["n_bridges_designed"]))
    print(f"\nwrote {a.out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
