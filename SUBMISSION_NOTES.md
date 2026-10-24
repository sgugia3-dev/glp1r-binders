# Submission notes

## Competition metadata is not recorded in this repository

**No deadline, track, cap or competition identifier can be determined from any file on disk.**
A search of the campaign tree for `deadline`, `track`, `proteinbase` and `adaptyv` returns one
hit, in `scripts/mk_p12.py`, and it is a Slurm wall-clock variable (`DEADLINE_S`), not a
submission date. The sibling EGFR repository records its own deadline and track in its notes;
this one has no equivalent and none is invented here. **Confirm the deadline, the track and the
design cap against the competition brief before submitting.**

---

## The file to upload

`SUBMISSION.csv` — **10 designs, 10 data rows plus a header, 40 columns.** Row 1 is the most
preferred; the file is already in rank order, so no reordering is needed.

### Column contract

The three required fields are present and first:

| column | value |
|---|---|
| `name` | `glp1r-01` … `glp1r-10` |
| `sequence` | single-chain amino-acid sequence, **140 aa for all ten**, all ten distinct |
| `molecule_class` | `protein` for all ten (single-chain, not nanobody/scFv/Fab) |

37 further columns carry provenance and metrics, so the file can be submitted as-is without
stripping:

| columns | what they hold |
|---|---|
| `rank` | 1-10, the submitted preference order; equals row order |
| `parent`, `source`, `design_id` | design lineage. `parent` is the backbone ancestor (`p13_04`, `p3_20`, …), `source` is its lineage (`p3`/`p10`/`p12`/`p13`), `design_id` is the full campaign name and the key into `metrics/seed_consistency.csv` and `metrics/p15_RANKED.csv` |
| `anchor` | the 10-residue agonist anchor, binder positions 1-10 |
| `n_seeds_gate`, `n_seeds_full` | how many of the three replicates hold the acceptance gate / the full hotspot set. 0-3 |
| `n_hotspots_total` | **how many hotspots that design's lineage has.** 3 or 4. Read every hotspot count against this — see below |
| `mean_ipsae`, `range_ipsae`, `mean_plddt`, `mean_motif_rmsd`, `range_motif_rmsd` | summaries over the three replicates. `range_*` is max − min |
| `seed1_*` … `seed3_*` (8 each) | the per-replicate values: `motif` (anchor RMSD, Å), `recept` (receptor RMSD, Å), `hits` (hotspots contacted), `sb` (interface salt bridges), `hc` (interface hydrophobic contacts), `ipsae`, `plddt` (binder chain), `iptm` |

`SUBMISSION.csv` is a superset of `metrics.csv` in the campaign's `glp1r_10_candidates/`
directory: the ten sequences match row for row, and `design_id` equals that file's `name`.

---

## Four things a reader needs before using the file

### 1. `n_hotspots_total` is 3 for five designs and 4 for five designs

The lineages targeted different receptor residues. Ranks 1-5 (p13 and p3 parents) carry a
**four**-hotspot set; ranks 6-10 (p10 and p12 parents) carry a **three**-hotspot set. A `hits`
value of 3 against `n_hotspots_total` = 3 is **full engagement**, and is not weaker than a 4
against 4. The epitope sets are tabulated in `README.md` section 1.

### 2. The three "seeds" are three predictions of identical input, not three designs

`seed1`, `seed2` and `seed3` are Boltz-2 runs of the same sequence against the same receptor,
differing only in random seed (campaign runs `p15_1099477`, `p18s2_1106925`, `p18s3_1106926`).
Disagreement between them is **prediction variance, not a property of the design**. They do not
vary the model, the MSA or the receptor construct, so they say nothing about whether the
predictor is right — only about how stable its answer is. `designs/<NN>/` holds all three
structures for each design.

### 3. Two designs do not hold everything in all three replicates, and are submitted anyway

| design | what happens |
|---|---|
| `glp1r-08` (rank 8, `p12_41`) | `n_seeds_gate` = **1**. Anchor RMSD 0.962 Å in replicate 1, 2.398 Å and 2.188 Å in replicates 2 and 3 — over the 2.0 Å gate in both. It keeps all 3 of its hotspots in every replicate, and its ipSAE (0.901-0.914) and binder pLDDT (0.915-0.918) are among the best in the set. |
| `glp1r-05` (rank 5, `p13_23`) | `n_seeds_full` = **2**. Drops from 4 hotspots to 3 in replicate 2; holds the gate in all three. |

Both are flagged in the file itself rather than silently dropped, because the evidence in
`analysis/replicate_consistency.md` is that **no single prediction identifies which designs will
behave this way** (replicate-1 ipSAE discriminates at AUC 0.487, binder pLDDT at 0.468), so
removing them on a criterion that cannot predict them would be a false economy. A reader who
wants the strictest set should take the 8 rows with `n_seeds_gate` = 3 **and** `n_seeds_full` = 3.

### 4. The ranking is a preference order, not a prediction

`rank` carries over the p15 composite score (anchor RMSD inverted, ipSAE, binder pLDDT, each
min-max normalised over the gated set, hotspot count ordered first). Two measured facts bound
how much weight it can carry:

- **ipSAE barely separates on-target from off-target on this receptor** — median 0.882 at 4/4
  hotspots against 0.871 at 0 hotspots over 1200 folds, and 0.894 against 0.870 over 480 folds
  (each round's `p*_best/README.txt` in the campaign tree). What establishes epitope engagement is the geometric
  per-hotspot check, measured from the designed region with the anchor excluded, not the score.
- **ipSAE reproduces across a reseed at Pearson r +0.309 to +0.450**, against +0.945 to +0.954
  for binder pLDDT.

The inputs are in `metrics/p15_RANKED.csv` (187 gated p15 refolds, with the `score` column and
its three normalised terms), so the order can be recomputed or re-weighted by any reader.

---

## Eligibility, as far as the files establish it

| check | status |
|---|---|
| de novo | **Yes, with one stated exception.** Backbones are RFdiffusion3; sequences are ProteinMPNN. Residues 11-140 — 130 of 140 — are designed. **Residues 1-10 are not de novo:** they are the GLP-1 agonist N-terminus (binder position 1 = GLP-1 His7), held fixed through ProteinMPNN and restrained to a 23-agonist consensus. All ten anchors differ from native `HAEGTFTSDV` (eight distinct anchors across the ten), but they are a natural-ligand motif, not an invention. This should be declared. |
| epitope | Extracellular domain, receptor residues 1-117. All nine hotspot residues used across the campaign (Q24, R25, L27, T42, L86, W87, L88, Q89, L95) lie inside it. In the predicted complexes a median 56% of receptor interface residues fall in 1-117 and the remainder in 118-409 (`metrics/interface_properties.csv`, 250 designs). |
| length | 140 aa, all ten. |
| molecule class | single-chain protein, all ten. |
| cysteines | **zero across all ten sequences** — no free-thiol or disulfide-scrambling liability. |
| metrics, structures, methodology | included: `metrics/`, `designs/`, `METHODS.md`. |
| design cap | **Not determinable from files on disk.** Ten designs is under every cap used in the sibling EGFR submission (20 or 40), but the cap for this challenge is not recorded here. |

---

## Known weaknesses, stated plainly

- **No experimental data of any kind.** Everything is Boltz-2 prediction or a sequence property
  computed from one.
- **ipSAE, the headline confidence number, is close to uninformative about epitope adherence on
  this target and reproduces poorly across a reseed.** See above and `README.md` section 5.
- **One design holds the gate in only 1 of 3 replicates** (`glp1r-08`) and one loses a hotspot in
  1 of 3 (`glp1r-05`).
- **All ten designs exceed the alanine ceiling of the only quantitative composition envelope
  available** (median 0.200 against 0.104), and 8 of 10 fall outside its hydrophobic band. That
  envelope comes from a different target and its own documentation says it should not be
  imported, so `analysis/composition.md` reports the comparison without treating it as a verdict.
- **Whole-chain unbound RMSD is large** (median 10.283 Å over 1-140) because the 10-residue
  agonist arm is unconstrained without the receptor. The folded body is stable — median 0.243 Å
  over residues 50-120. `analysis/controls.md` explains which window to read and why.
- **Self-association is low but was never selected against.** Median self ipTM 0.262 over 210
  binders; the ten run 0.123-0.467.
- **No length exploration.** Every design in the campaign is exactly 140 aa.
- **Competition deadline, track and cap are not recorded anywhere in this repository.**
