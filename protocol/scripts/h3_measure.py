#!/usr/bin/env python
"""Stage 3b -- deviation of each refold from the AlphaFold3 predicted structure.

Method. Superpose the refold onto the AF3 model using the RECEPTOR CA atoms only
(Kabsch), then, with that single transform, measure where the ten motif residues ended
up relative to where AF3 put the peptide's first ten residues. Superposing on the
receptor and scoring the motif keeps the two independent: an identical fold gives ~0 A
and a binder that dragged the motif elsewhere gives a large number.

Chain assignment is resolved by LENGTH and sequence identity, never by label: Boltz
renumbers chains from 1 and the binder's label depends on YAML order.
"""
import argparse, glob, json, math, os, sys

A3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
      "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
      "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}


def read_ca(path):
    """{chain: {resnum: (x,y,z)}} plus {chain: seq-string} using CA atoms only."""
    hdr, rows, inl = [], [], False
    for line in open(path):
        if line.startswith("_atom_site."):
            hdr.append(line.strip().split(".")[1]); inl = True; continue
        if inl and line.startswith(("ATOM", "HETATM")): rows.append(line.split())
        elif inl and line.startswith("#"): break
    i = {h: k for k, h in enumerate(hdr)}
    need = ("auth_asym_id", "auth_seq_id", "label_comp_id", "label_atom_id",
            "Cartn_x", "Cartn_y", "Cartn_z")
    if any(n not in i for n in need): return {}, {}
    ca, seq = {}, {}
    for r in rows:
        if len(r) <= max(i[n] for n in need): continue
        if r[i["label_atom_id"]] != "CA": continue
        try: n = int(r[i["auth_seq_id"]])
        except ValueError: continue
        c = r[i["auth_asym_id"]]
        ca.setdefault(c, {})[n] = (float(r[i["Cartn_x"]]), float(r[i["Cartn_y"]]),
                                   float(r[i["Cartn_z"]]))
        seq.setdefault(c, {})[n] = A3.get(r[i["label_comp_id"]], "X")
    return ca, {c: "".join(v[k] for k in sorted(v)) for c, v in seq.items()}


def kabsch(P, Q):
    """Return (R, tP, tQ) aligning P onto Q. Pure python, no numpy dependency games."""
    n = len(P)
    cp = [sum(p[k] for p in P) / n for k in range(3)]
    cq = [sum(q[k] for q in Q) / n for k in range(3)]
    P = [[p[k] - cp[k] for k in range(3)] for p in P]
    Q = [[q[k] - cq[k] for k in range(3)] for q in Q]
    H = [[sum(P[m][r] * Q[m][c] for m in range(n)) for c in range(3)] for r in range(3)]
    import numpy as np
    U, S, Vt = np.linalg.svd(np.array(H))
    d = 1.0 if np.linalg.det(Vt.T @ U.T) > 0 else -1.0
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, cp, cq


def apply(R, cp, cq, p):
    import numpy as np
    v = np.array([p[k] - cp[k] for k in range(3)])
    w = R @ v
    return tuple(float(w[k]) + cq[k] for k in range(3))


def rmsd(A, B):
    return math.sqrt(sum(sum((a[k] - b[k]) ** 2 for k in range(3))
                         for a, b in zip(A, B)) / len(A))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred-dir", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--af3-model", required=True)
    ap.add_argument("--out-csv", required=True)
    a = ap.parse_args()

    man = json.load(open(a.manifest))
    rec_len = man["receptor_len"]
    by = {y["name"]: y for y in man["yamls"]}

    ref_ca, ref_seq = read_ca(a.af3_model)
    # AF3 model: chain A = receptor (409), chain B = peptide (30)
    rA = max(ref_ca, key=lambda c: len(ref_ca[c]))
    rB = [c for c in ref_ca if c != rA]
    if not rB: print("FATAL: AF3 model has one chain"); return 2
    rB = rB[0]
    print(f"AF3 model: receptor chain {rA} ({len(ref_ca[rA])} CA), peptide chain {rB} "
          f"({len(ref_ca[rB])} CA)")

    rows, missing = [], 0
    for cif in sorted(glob.glob(os.path.join(a.pred_dir, "*", "*_model_0.cif"))):
        name = os.path.basename(os.path.dirname(cif))
        rec = by.get(name)
        if not rec: continue
        ca, seq = read_ca(cif)
        if not ca: missing += 1; continue
        # identify chains by LENGTH, then confirm the binder by sequence identity
        tgt = max(ca, key=lambda c: len(ca[c]))
        bnd = [c for c in ca if c != tgt]
        if not bnd: missing += 1; continue
        bnd = bnd[0]
        ok_len = len(ca[tgt]) == rec_len and len(ca[bnd]) == rec["binder_len"]
        ok_seq = seq.get(bnd, "") == rec["binder_seq"]
        # receptor-only superposition
        common = sorted(set(ca[tgt]) & set(ref_ca[rA]))
        if len(common) < 50: missing += 1; continue
        R, cp, cq = kabsch([ca[tgt][k] for k in common], [ref_ca[rA][k] for k in common])
        tgt_rms = rmsd([apply(R, cp, cq, ca[tgt][k]) for k in common],
                       [ref_ca[rA][k] for k in common])
        # motif: binder positions rec["motif"] vs AF3 peptide residues 1..10
        mp = sorted(rec["motif"])
        ref_pep = sorted(ref_ca[rB])[: len(mp)]
        have = [k for k in mp if k in ca[bnd]]
        if len(have) != len(mp): missing += 1; continue
        motif_rms = rmsd([apply(R, cp, cq, ca[bnd][k]) for k in mp],
                         [ref_ca[rB][k] for k in ref_pep])
        rows.append({"name": name, "tag": rec["tag"], "binder_len": rec["binder_len"],
                     "motif_lo": mp[0], "motif_hi": mp[-1],
                     "receptor_rmsd": round(tgt_rms, 3),
                     "motif_rmsd_vs_af3": round(motif_rms, 3),
                     "chains_len_ok": ok_len, "binder_seq_ok": ok_seq,
                     "binder_chain": bnd, "target_chain": tgt})

    if not rows:
        print(f"FATAL: nothing measurable (skipped {missing})"); return 3
    rows.sort(key=lambda r: r["motif_rmsd_vs_af3"])
    cols = list(rows[0].keys())
    with open(a.out_csv, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows: fh.write(",".join(str(r[c]) for c in cols) + "\n")

    print(f"\nmeasured {len(rows)} refolds, skipped {missing}")
    bad = [r for r in rows if not (r["chains_len_ok"] and r["binder_seq_ok"])]
    print(f"chain/sequence identity verified: {len(rows)-len(bad)}/{len(rows)}")
    print(f"\n{'design':44s} {'len':>4s} {'recRMSD':>8s} {'motifRMSD':>10s}")
    for r in rows:
        print(f"{r['name'][20:64]:44s} {r['binder_len']:4d} "
              f"{r['receptor_rmsd']:8.2f} {r['motif_rmsd_vs_af3']:10.2f}")
    import statistics as st
    v = [r["motif_rmsd_vs_af3"] for r in rows]
    print(f"\nmotif RMSD vs AF3: min {min(v):.2f}  median {st.median(v):.2f}  max {max(v):.2f} A")
    print(f"wrote {a.out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
