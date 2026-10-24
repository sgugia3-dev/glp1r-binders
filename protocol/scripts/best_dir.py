#!/usr/bin/env python
"""Assemble a best-candidate deliverable directory for one GLP-1R campaign.

Everything that differs between campaigns - name, hotspots, totals, gate rate,
ordering - is an argument or derived from the data, so the README in a given
directory always describes that directory.

Selection:
  1. GATE  already applied in the ranked CSV: motif RMSD < 2 A, >= 3 hotspots,
     receptor RMSD < 5 A.
  2. SCORE already applied: mean of min-max normalised motif RMSD (inverted),
     ipSAE at 15/15, and binder-chain mean pLDDT.
  3. ORDER  'score' or 'tier' (hotspot count first, then score). Tier matters
     when full engagement is scarce, because the score ignores hotspot count.
  4. DEDUPLICATE  greedy clustering over the designed region, so each entry is
     a distinct backbone rather than one design refolded under another seed.
"""
import argparse, csv, glob, os, re, shutil, sys

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}

ap = argparse.ArgumentParser()
ap.add_argument("--ranked", required=True)
ap.add_argument("--root", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--prefix", required=True, help="tag prefix, e.g. p3")
ap.add_argument("--campaign", required=True, help="run id, e.g. p3_1088996")
ap.add_argument("--hotspots", required=True, help="e.g. A86:...,A87:...")
ap.add_argument("--contig", default="")
ap.add_argument("--hotspot-count", type=int, default=None,
                help="number of hotspots to call full engagement. Needed when "
                     "--hotspots is prose or the round mixes sources with "
                     "different set sizes, where counting commas is wrong.")
ap.add_argument("--analysis", default="", help="dir holding the deviation/epitope CSVs")
ap.add_argument("--motif-len", type=int, default=10)
ap.add_argument("--topn", type=int, default=50)
ap.add_argument("--identity", type=float, default=0.70)
ap.add_argument("--order", choices=["score", "tier"], default="tier")
ap.add_argument("--ipsae-note", default="")
a = ap.parse_args()
ML = a.motif_len


def read_cif(path):
    names, rows = [], []
    for line in open(path):
        if line.startswith("_atom_site."):
            names.append(line.strip().split(".")[1]); continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")): continue
        t = line.split()
        if len(t) == len(names): rows.append(dict(zip(names, t)))
    return rows


def chain_of(r):
    return r.get("auth_asym_id") or r.get("label_asym_id")


def seq_of(rows, ch):
    seen, out = set(), []
    for r in rows:
        if chain_of(r) != ch or r.get("label_atom_id") != "CA": continue
        rn = r.get("auth_seq_id") or r.get("label_seq_id")
        if rn in seen: continue
        seen.add(rn); out.append(AA3.get((r.get("label_comp_id") or "").upper(), "X"))
    return "".join(out)


def write_binder_pdb(rows, ch, path):
    n = 0
    with open(path, "w") as fh:
        fh.write("REMARK   binder chain only, extracted from the Boltz complex\n")
        fh.write("REMARK   B-factor column holds per-residue pLDDT (0-100)\n")
        for r in rows:
            if chain_of(r) != ch: continue
            n += 1
            el = (r.get("type_symbol") or "").strip()
            an = (r.get("label_atom_id") or "").strip().strip('"')
            nm = an.ljust(4) if len(an) >= 4 or len(el) > 1 else " " + an.ljust(3)
            try: b = float(r.get("B_iso_or_equiv", 0) or 0)
            except ValueError: b = 0.0
            try: rs = int(r.get("auth_seq_id") or r.get("label_seq_id"))
            except (TypeError, ValueError): rs = n
            fh.write("ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
                     % (n, nm, (r.get("label_comp_id") or "UNK")[:3], "A", rs,
                        float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]),
                        1.0, b, el[:2].rjust(2)))
        fh.write("TER\nEND\n")
    return n


def ident(x, y):
    n = min(len(x), len(y))
    return sum(p == q for p, q in zip(x[:n], y[:n])) / n if n else 0.0


rows = list(csv.DictReader(open(a.ranked)))
n_gated = len(rows)
if a.order == "tier":
    rows.sort(key=lambda r: (-int(r["hotspots"]), -float(r["score"])))
else:
    rows.sort(key=lambda r: -float(r["score"]))
print(f"ordering: {a.order}   gated designs in ranking: {n_gated}")

# campaign-wide totals, straight from the measurement tables
n_total, hist = 0, {}
if a.analysis:
    for f in sorted(glob.glob(os.path.join(a.analysis, f"{a.campaign}_sh*_epitope.csv"))):
        for r in csv.DictReader(open(f)):
            n_total += 1
            h = int(r["n_hotspots_hit"]); hist[h] = hist.get(h, 0) + 1
print(f"campaign refolds: {n_total}   hits histogram: {dict(sorted(hist.items(), reverse=True))}")

os.makedirs(os.path.join(a.out, "structures"), exist_ok=True)
os.makedirs(os.path.join(a.out, "binder_only"), exist_ok=True)
for d in ("structures", "binder_only"):
    for old in glob.glob(os.path.join(a.out, d, "*")):
        os.remove(old)

kept, clusters, scanned, missing = [], [], 0, 0
for r in rows:
    if len(kept) >= a.topn: break
    scanned += 1
    d = os.path.join(a.root, "sh" + r["shard"], "boltz_results_yaml", "predictions", r["name"])
    cif = glob.glob(os.path.join(d, "*_model_0.cif"))
    if not cif: missing += 1; continue
    at = read_cif(cif[0])
    if not at: missing += 1; continue
    ch = chain_of(at[0])
    s = seq_of(at, ch)
    if not s: missing += 1; continue
    des = s[ML:]
    hit = next((c for c in clusters if ident(des, c["seq"][ML:]) >= a.identity), None)
    if hit: hit["n"] += 1; continue
    rank = len(kept) + 1
    # rfd3 names designs from the input JSON key, which doubles the run
    # prefix (prod_r3_prod_r3_m1_..., w2_r86_89_w2_r86_89_m1_...). Strip one copy.
    short = re.sub(r"^(.+?)_\1_", "", r["name"])
    tag = "%s_%02d_sh%s_%s" % (a.prefix, rank, r["shard"], short)
    shutil.copy(cif[0], os.path.join(a.out, "structures", tag + "_complex.cif"))
    na = write_binder_pdb(at, ch, os.path.join(a.out, "binder_only", tag + "_binder.pdb"))
    r2 = dict(r); r2.update({"rank": rank, "tag": tag, "seq": s, "cluster_n": 1,
                             "binder_atoms": na})
    clusters.append({"seq": s, "n": 1, "r": r2}); kept.append(r2)
for c in clusters: c["r"]["cluster_n"] = c["n"]
print(f"scanned {scanned} ranked rows -> kept {len(kept)} distinct backbones "
      f"(identity cutoff {a.identity:.0%}); {missing} unreadable")

cols = ["rank", "tag", "score", "shard", "name", "cluster_n", "binder_len",
        "motif_rmsd", "receptor_rmsd", "hotspots", "min_dist", "ipsae", "ipsae_cons",
        "plddt_binder", "plddt_design", "plddt_complex", "iplddt", "iptm", "pdockq", "seq"]
with open(os.path.join(a.out, "metrics.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(cols)
    for r in kept: w.writerow([r.get(c, "") for c in cols])
with open(os.path.join(a.out, "candidates.fasta"), "w") as fh:
    for r in kept:
        fh.write(">%s score=%s motif_rmsd=%s hotspots=%s ipsae=%s plddt_binder=%s len=%d\n"
                 % (r["tag"], r["score"], r["motif_rmsd"], r["hotspots"], r["ipsae"],
                    r["plddt_binder"], len(r["seq"])))
        for i in range(0, len(r["seq"]), 60): fh.write(r["seq"][i:i + 60] + "\n")
with open(os.path.join(a.out, "designed_region.fasta"), "w") as fh:
    for r in kept:
        fh.write(">%s designed_region\n%s\n" % (r["tag"], r["seq"][ML:]))

nhot = a.hotspot_count or len(a.hotspots.split(","))
# count full engagement against THIS campaign's hotspot count, not a literal 4:
# P10 has three hotspots, so "== 4" would report zero fully-engaged designs.
n4 = sum(1 for r in kept if int(r["hotspots"]) == nhot)
motifs = set(r["seq"][:ML] for r in kept)
gate_pct = (100.0 * n_gated / n_total) if n_total else float("nan")
ordering = ("hotspot count first, then composite score"
            if a.order == "tier" else "composite score")
hist_s = "  ".join(f"{h}:{hist[h]}" for h in sorted(hist, reverse=True)) or "n/a"
with open(os.path.join(a.out, "README.txt"), "w") as fh:
    fh.write(f"""{a.prefix.upper()} best candidates   (GLP-1R, campaign {a.campaign})
{'=' * 62}
{len(kept)} distinct backbones, ordered by {ordering}.

CONTENTS
  metrics.csv            all metrics, one row per candidate, best first
  candidates.fasta       full binder sequences
  designed_region.fasta  the designed residues only (first {ML}-aa anchor stripped)
  structures/            Boltz complex, binder chain A + GLP1R chain B, mmCIF
  binder_only/           binder chain as PDB; B-factor column = pLDDT (0-100)

PROTOCOL
  rfd3 contig {a.contig or '(see inputs)'}
  atom-level hotspots: {a.hotspots}
  then ProteinMPNN with the {ML}-residue anchor fixed,
  then Boltz-2 with a cached receptor MSA.

SELECTION
  1. GATE   {n_gated}/{n_total} refolds passed ({gate_pct:.0f}%)
       motif RMSD vs the AF3 peptide < 2.0 A
       >= 3 of {nhot} hotspots contacted, measured geometrically from the
         DESIGNED region only, with the anchor motif excluded
       receptor RMSD < 5.0 A
  2. SCORE  mean of three min-max normalised terms over the gated set:
       motif RMSD (inverted), ipSAE at 15/15 cutoffs, binder-chain mean pLDDT
  3. ORDER  {ordering}
  4. DEDUPLICATE  greedy clustering at {a.identity:.0%} identity over the designed
     region, keeping the best-ordered member. cluster_n in metrics.csv is how
     many gated designs collapsed into that entry.

  {n4}/{len(kept)} of these contact all {nhot} hotspots.
  Campaign-wide hotspot-count histogram: {hist_s}
  Anchor sequence in this set: {', '.join(sorted(motifs))}

READ THE RANKING WITH THIS CAVEAT
  ipSAE barely detects epitope adherence on this target, so it must not be used
  to decide whether a design is on-target - only to order designs that the
  geometric hotspot check in step 1 has already verified.
{('  Measured on this campaign: ' + a.ipsae_note) if a.ipsae_note else ''}
  Binder pLDDT is the best-behaved of the confidence metrics here, but it is
  still weak. Step 1, not step 2, is what establishes epitope engagement.

  complex_plddt is reported but was NOT ranked on: the complex is 409 receptor
  residues to ~140 binder residues, so it is three-quarters a constant.

  ipsae_cons is the conservative ipSAE direction; ipsae is the optimistic `max`
  row. The gap is large for some designs - check it before trusting one value.
""")
print("\n%4s %-40s %6s %6s %5s %6s %7s %4s"
      % ("rank", "tag", "score", "motif", "hits", "ipSAE", "pLDDTb", "dup"))
for r in kept[:30]:
    print("%4d %-40s %6.3f %6.2f %5s %6.3f %7.3f %4d"
          % (r["rank"], r["tag"][:40], float(r["score"]), float(r["motif_rmsd"]),
             r["hotspots"], float(r["ipsae"]), float(r["plddt_binder"]), r["cluster_n"]))
print(f"\nwrote {a.out}")
