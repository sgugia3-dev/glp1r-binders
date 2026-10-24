#!/usr/bin/env python
"""Assemble the P2 best-candidate deliverable directory, ordered by the ranking.

Takes p2_RANKED.csv (already gated: motif RMSD < 2 A, >= 3 hotspots, receptor
RMSD < 5 A; scored as the mean of min-max normalised motif RMSD / ipSAE /
binder pLDDT), then greedily clusters by sequence identity over the designed
region so that each entry in the directory is a DISTINCT backbone rather than
the same design refolded under a different seed.

Writes, for each kept design: the Boltz complex, a binder-only PDB for
downstream work, and a row in metrics.csv. Sequences go to one FASTA.
"""
import csv, glob, os, shutil, sys

ranked, root, out, topn, ident_cut = (sys.argv[1], sys.argv[2], sys.argv[3],
                                      int(sys.argv[4]), float(sys.argv[5]))
# 'tier' orders by hotspot count first, then score. Use it when full engagement
# is the scarce property and the composite score, which does not include hotspot
# count, would otherwise bury the fully-engaged designs.
SORT = sys.argv[6] if len(sys.argv) > 6 else "score"
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
MOTIF_LEN = 10


def read_cif(path):
    names, rows = [], []
    for line in open(path):
        if line.startswith("_atom_site."):
            names.append(line.strip().split(".")[1]); continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")): continue
        t = line.split()
        if len(t) != len(names): continue
        rows.append(dict(zip(names, t)))
    return rows


def binder_chain(rows):
    for r in rows:
        return r.get("auth_asym_id") or r.get("label_asym_id")
    return None


def seq_of(rows, ch):
    seen, out_s = set(), []
    for r in rows:
        c = r.get("auth_asym_id") or r.get("label_asym_id")
        if c != ch: continue
        if r.get("label_atom_id") != "CA": continue
        rn = r.get("auth_seq_id") or r.get("label_seq_id")
        if rn in seen: continue
        seen.add(rn)
        out_s.append(AA3.get((r.get("label_comp_id") or "").upper(), "X"))
    return "".join(out_s)


def write_binder_pdb(rows, ch, path):
    n = 0
    with open(path, "w") as fh:
        fh.write("REMARK   binder chain only, extracted from the Boltz complex\n")
        fh.write("REMARK   B-factor column holds per-residue pLDDT (0-100)\n")
        for r in rows:
            c = r.get("auth_asym_id") or r.get("label_asym_id")
            if c != ch: continue
            n += 1
            el = (r.get("type_symbol") or "").strip()
            an = (r.get("label_atom_id") or "").strip().strip('"')
            nm = an.ljust(4) if len(an) >= 4 or len(el) > 1 else " " + an.ljust(3)
            try:
                b = float(r.get("B_iso_or_equiv", 0) or 0)
            except ValueError:
                b = 0.0
            try:
                rs = int(r.get("auth_seq_id") or r.get("label_seq_id"))
            except (TypeError, ValueError):
                rs = n
            fh.write("ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
                     % (n, nm, (r.get("label_comp_id") or "UNK")[:3], "A", rs,
                        float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]),
                        1.0, b, el[:2].rjust(2)))
        fh.write("TER\nEND\n")
    return n


def ident(a, b):
    n = min(len(a), len(b))
    return sum(x == y for x, y in zip(a[:n], b[:n])) / n if n else 0.0


rows = list(csv.DictReader(open(ranked)))
if SORT == "tier":
    rows.sort(key=lambda r: (-int(r["hotspots"]), -float(r["score"])))
else:
    rows.sort(key=lambda r: -float(r["score"]))
print(f"ordering: {SORT}")
print(f"gated designs in ranking: {len(rows)}")

os.makedirs(os.path.join(out, "structures"), exist_ok=True)
os.makedirs(os.path.join(out, "binder_only"), exist_ok=True)

kept, clusters, scanned, missing = [], [], 0, 0
for r in rows:
    if len(kept) >= topn: break
    scanned += 1
    d = os.path.join(root, "sh" + r["shard"], "boltz_results_yaml", "predictions", r["name"])
    cif = glob.glob(os.path.join(d, "*_model_0.cif"))
    if not cif:
        missing += 1; continue
    at = read_cif(cif[0])
    ch = binder_chain(at)
    s = seq_of(at, ch)
    if not s: missing += 1; continue
    des = s[MOTIF_LEN:]
    dup = False
    for c in clusters:
        if ident(des, c["seq"][MOTIF_LEN:]) >= ident_cut:
            c["n"] += 1; dup = True; break
    if dup: continue
    rank = len(kept) + 1
    tag = "p2_%02d_sh%s_%s" % (rank, r["shard"], r["name"].replace("prod_r3_prod_r3_", ""))
    shutil.copy(cif[0], os.path.join(out, "structures", tag + "_complex.cif"))
    na = write_binder_pdb(at, ch, os.path.join(out, "binder_only", tag + "_binder.pdb"))
    r2 = dict(r); r2.update({"rank": rank, "tag": tag, "seq": s, "cluster_n": 1,
                             "binder_atoms": na})
    clusters.append({"seq": s, "n": 1, "r": r2})
    kept.append(r2)

for c in clusters:
    c["r"]["cluster_n"] = c["n"]
print(f"scanned {scanned} ranked rows -> kept {len(kept)} distinct backbones "
      f"(identity cutoff {ident_cut:.0%}); {missing} unreadable")

cols = ["rank", "tag", "score", "shard", "name", "cluster_n", "binder_len",
        "motif_rmsd", "receptor_rmsd", "hotspots", "min_dist", "ipsae",
        "ipsae_cons", "plddt_binder", "plddt_design", "plddt_complex", "iplddt",
        "iptm", "pdockq", "seq"]
with open(os.path.join(out, "metrics.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(cols)
    for r in kept: w.writerow([r.get(c, "") for c in cols])

with open(os.path.join(out, "candidates.fasta"), "w") as fh:
    for r in kept:
        fh.write(">%s score=%s motif_rmsd=%s hotspots=%s ipsae=%s plddt_binder=%s len=%d\n"
                 % (r["tag"], r["score"], r["motif_rmsd"], r["hotspots"], r["ipsae"],
                    r["plddt_binder"], len(r["seq"])))
        for i in range(0, len(r["seq"]), 60):
            fh.write(r["seq"][i:i + 60] + "\n")

# designed region only, which is what goes to synthesis if the anchor is grafted
with open(os.path.join(out, "designed_region.fasta"), "w") as fh:
    for r in kept:
        fh.write(">%s designed_region\n%s\n" % (r["tag"], r["seq"][MOTIF_LEN:]))

mot = set(r["seq"][:MOTIF_LEN] for r in kept)
n4 = sum(1 for r in kept if int(r["hotspots"]) == 4)
with open(os.path.join(out, "README.txt"), "w") as fh:
    fh.write(f"""P2 best candidates  (GLP-1R, campaign p2_1088883)
=================================================
{len(kept)} distinct backbones, ordered by the composite ranking.

CONTENTS
  metrics.csv            all metrics, one row per candidate, best first
  candidates.fasta       full 140-aa binder sequences
  designed_region.fasta  the 130 designed residues only (anchor stripped)
  structures/            Boltz complex, binder chain A + GLP1R chain B, mmCIF
  binder_only/           binder chain as PDB; B-factor column = pLDDT (0-100)

HOW THESE WERE CHOSEN
  Campaign: 1200 refolds. Protocol = rfd3 contig B1-10,130-130,/0,A1-409 with
  atom-level hotspots A24:OE1,NE2  A25:NH1,NH2  A27:CD1  A42:OG1, no potential;
  ProteinMPNN with the 10-residue anchor fixed; Boltz-2 with a cached receptor MSA.

  1. GATE (869/1200 passed, 72%)
       motif RMSD vs the AF3 peptide < 2.0 A
       >= 3 of 4 hotspots contacted, measured geometrically from the DESIGNED
         region only, with the anchor motif excluded
       receptor RMSD < 5.0 A
  2. SCORE  mean of three min-max normalised terms over the gated set:
       motif RMSD (inverted), ipSAE at 15/15 cutoffs, binder-chain mean pLDDT
  2b. ORDER  {"by hotspot count, then score" if SORT == "tier" else "by score"}
  3. DEDUPLICATE  greedy clustering at {ident_cut:.0%} identity over the designed
     region, keeping the highest-scoring member. cluster_n in metrics.csv is how
     many gated designs collapsed into that entry.

  {n4}/{len(kept)} contact all 4 hotspots. Anchor is identical in {len(kept)-len(mot)+1}/{len(kept)} ({sorted(mot)[0]}).

READ THE RANKING WITH THIS CAVEAT
  ipSAE does not detect epitope adherence on this target. Measured over all 1200
  folds: designs contacting 4/4 hotspots have a median ipSAE of 0.882, designs
  contacting NONE have 0.871 - a 0.011 separation. Binder pLDDT is the best of
  the confidence metrics here (0.880 vs 0.845) and is monotonic in hotspot count,
  but it is still weak. The geometric hotspot check in step 1, not the score in
  step 2, is what establishes that a design is on-target. Treat the score as an
  ordering within an already-verified set.

  Related: complex_plddt is reported but was NOT ranked on. The complex is 409
  receptor residues to 140 binder residues, so it is three-quarters a constant.

  ipsae_cons is the conservative ipSAE direction; ipsae is the optimistic `max`
  row. The gap between them is large for some designs - check it before trusting
  a single ipSAE value.
""")
print("\n%4s %-34s %6s %6s %4s %6s %7s %4s" % ("rank", "tag", "score", "motif",
                                               "hits", "ipSAE", "pLDDTb", "dup"))
for r in kept[:30]:
    print("%4d %-34s %6.3f %6.2f %4s %6.3f %7.3f %4d"
          % (r["rank"], r["tag"][:34], float(r["score"]), float(r["motif_rmsd"]),
             r["hotspots"], float(r["ipsae"]), float(r["plddt_binder"]), r["cluster_n"]))
print(f"\nwrote {out}")
