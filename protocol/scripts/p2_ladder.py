#!/usr/bin/env python
"""Is the i+4 ladder MPNN recovered (ext 4/8/12/16) the receptor-facing stripe?

MPNN re-derived native GLP-1 residues at designed-region positions 4,5,8,9,12,16,
which is a helical face. If that face is the one packed against GLP1R, the recovery
is structural; if it points at solvent, it is just MPNN's A/L/F compositional bias.
Measures, per candidate, each ext position's minimum heavy-atom distance to the
receptor chain and its relative burial.
"""
import csv, glob, os, sys, statistics as st
import numpy as np

cand, outdir = sys.argv[1], sys.argv[2]
NATIVE_EXT = "SSYLEGQAAKEFIAWLVKGR"


def atoms(path):
    names, rec = [], []
    for line in open(path):
        if line.startswith("_atom_site."):
            names.append(line.strip().split(".")[1]); continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")): continue
        t = line.split()
        if len(t) != len(names): continue
        d = dict(zip(names, t))
        if d.get("type_symbol") == "H": continue
        try:
            xyz = (float(d["Cartn_x"]), float(d["Cartn_y"]), float(d["Cartn_z"]))
        except (KeyError, ValueError):
            continue
        ch = d.get("auth_asym_id") or d.get("label_asym_id")
        rn = d.get("auth_seq_id") or d.get("label_seq_id")
        rec.append((ch, rn, d.get("label_comp_id"), d.get("label_atom_id"), xyz))
    return rec


rows = list(csv.DictReader(open(cand)))
per_pos = {i: [] for i in range(1, 21)}
contact_counts = {i: 0 for i in range(1, 21)}
n_done = 0
for r in rows:
    tag = "rank%02d_sh%s_%s" % (int(r["rank"]), r["shard"],
                                r["name"].replace("prod_r3_prod_r3_", ""))
    f = os.path.join(outdir, tag + ".cif")
    if not os.path.exists(f): continue
    rec = atoms(f)
    chains = []
    for ch, rn, cm, an, xyz in rec:
        if not chains or chains[-1] != ch: chains.append(ch)
    binder_ch = chains[0]
    tgt = np.array([a[4] for a in rec if a[0] != binder_ch])
    if not len(tgt): continue
    # binder residues in order
    order, seen = [], set()
    for ch, rn, cm, an, xyz in rec:
        if ch != binder_ch: continue
        if (ch, rn) in seen: continue
        seen.add((ch, rn)); order.append(rn)
    for i in range(1, 21):
        rn = order[10 + i - 1]
        pts = np.array([a[4] for a in rec if a[0] == binder_ch and a[1] == rn])
        d = float(np.sqrt(((pts[:, None, :] - tgt[None, :, :]) ** 2).sum(-1)).min())
        per_pos[i].append(d)
        if d < 5.0: contact_counts[i] += 1
    n_done += 1

print(f"analysed {n_done} candidates\n")
print(f"{'ext':>4s} {'native':>6s} {'recov':>6s} {'med min-dist to receptor':>26s} "
      f"{'<5A in':>8s}")
rec_counts = {}
for i in range(1, 21):
    col = [r["seq"][10 + i - 1] for r in rows]
    rec_counts[i] = sum(c == NATIVE_EXT[i - 1] for c in col)
for i in range(1, 21):
    d = per_pos[i]
    if not d: continue
    star = "  <== i+4 ladder" if i in (4, 8, 12, 16) else ""
    print(f"{i:4d} {NATIVE_EXT[i-1]:>6s} {rec_counts[i]:3d}/{len(rows):<3d} "
          f"{st.median(d):21.2f} A {contact_counts[i]:5d}/{n_done}{star}")

lad = [i for i in (4, 8, 12, 16)]
other = [i for i in range(1, 21) if i not in lad]
ld = [x for i in lad for x in per_pos[i]]
od = [x for i in other for x in per_pos[i]]
print(f"\nladder positions 4/8/12/16 : median min-dist {st.median(ld):.2f} A, "
      f"{100*sum(x<5 for x in ld)/len(ld):.0f}% in contact")
print(f"all other ext positions    : median min-dist {st.median(od):.2f} A, "
      f"{100*sum(x<5 for x in od)/len(od):.0f}% in contact")
hi = [i for i in range(1, 21) if rec_counts[i] >= 11]
lo = [i for i in range(1, 21) if rec_counts[i] <= 4]
hd = [x for i in hi for x in per_pos[i]]
od2 = [x for i in lo for x in per_pos[i]]
print(f"\nhigh-recovery positions {hi}: median min-dist {st.median(hd):.2f} A")
print(f"low-recovery  positions {lo}: median min-dist {st.median(od2):.2f} A")
