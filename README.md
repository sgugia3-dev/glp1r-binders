# De novo binders to the human GLP-1 receptor

Ten de novo designs against human GLP-1R, each validated by **three independent Boltz-2
predictions of byte-identical input**. Backbones come from RFdiffusion3 built around the
GLP-1 agonist N-terminus, sequences from ProteinMPNN with that anchor held fixed, and
validation from unrestrained Boltz-2 co-folding against the full 409-residue receptor chain.

Everything below is computed from files in this repository or in the campaign tree it was
assembled from (`/lustre/home/gsgueglia/challenge2026/glp1r2026/`). Where a figure could not
be sourced, the document says so rather than supplying one.

---

## 1. Objective

Bind the extracellular face of the human GLP-1 receptor with a single-chain protein that
carries the agonist peptide's own N-terminal anchor.

The receptor chain used throughout is **409 residues** (`metrics/receptor_properties.csv`,
region `receptor FULL 1-409`). Its **extracellular domain, residues 1-117**, is the designable
face; the remaining 118-409 is the transmembrane domain. The three regions are not alike —
the ECD is acidic (pI 4.88, charge −5.18 at pH 7) while the TMD is basic (pI 8.96, +7.39),
so the full chain nets out at pI 7.93 and +2.45. Everything the designs are meant to touch
sits in the acidic ECD.

**Binder architecture: 140 residues in two parts.**

| residues | role |
|---|---|
| 1-10 | **the anchor** — the GLP-1 agonist N-terminus. Binder position 1 = GLP-1 His7, so the anchor spans GLP-1 7-16; the native sequence is `HAEGTFTSDV` (`protocol/scripts/consensus.py`). Held fixed during ProteinMPNN. |
| 11-140 | **the designed region** — 130 residues, free to ProteinMPNN, and the only part the epitope gate is allowed to score. |

Every epitope measurement in this campaign is made **from the designed region alone, with the
anchor excluded**, so a design cannot earn an epitope hit for contacts its fixed agonist arm
makes anyway. That rule is stated in every round's own `README.txt`.

### The epitope sets differ between lineages — read any hotspot count against its own set

This matters more here than anything else in the results table. Different rounds targeted
different receptor residues, and the number of hotspots in the gate is **not** the same number
for every design. A design reported as hitting **3 of 3 is not weaker than one hitting 4 of 4**.

| lineage | receptor hotspot residues | residue identities | hotspots in that lineage's gate |
|---|---|---|---|
| p2 | A24, A25, A27, A42 | Q24, R25, L27, T42 | its README is inconsistent — see section 4 |
| p3 | A86, A87, A88, A89 | L86, W87, L88, Q89 | **4** |
| p10 | A27, A42, A86 | L27, T42, L86 | **3** |
| p12 | A27, A42, A86 | L27, T42, L86 | **3** |
| p13 | A27, A42, A86, A95 | L27, T42, L86, L95 | **4** |
| p14, p15 | inherited per parent: p3 → A86-A89; p10/p12 → A27/A42/A86; p13 → A27/A42/A86/A95 | | **3 or 4, per parent** |

Hotspot residue identities read off `metrics/receptor_properties.csv`; hotspot sets off each
round's `p*_best/README.txt` in the campaign tree. All nine hotspot residues lie inside the ECD (1-117).

Of the ten submitted designs, **five carry a four-hotspot set** (ranks 1-5, p13 and p3 parents)
and **five carry a three-hotspot set** (ranks 6-10, p10 and p12 parents) — the column
`n_hotspots_total` in `SUBMISSION.csv` records which. This trap had bitten the campaign
tree itself: `p15_best/README.txt` previously read "34/50 of these contact all 4 hotspots",
because the packaging script was given a fixed count of 4 for a directory that mixes
lineages. Those 34 were exactly the 14 p3 and 20 p13 entries, i.e. **every p15 design that
has four hotspots at all**; the other 16 are p10/p12 designs with a three-hotspot set and
all 16 contact all three. That file now reports 50/50 against each design's own set.
**Nothing in that round failed to reach its own full set.**

---

## 2. What is in this repository

```
SUBMISSION.csv          the ten designs, ranked (name, sequence, molecule_class + 37 metric columns)
METHODS.md              the protocol, stage by stage
designs/                one directory per design, three PDB structures each (seed1/2/3),
                        chain A = binder 1-140, chain B = GLP1R 1-409, B-factor = pLDDT x100
designs/sequences.fasta the ten sequences
metrics/                the nine evidence tables, indexed in metrics/README.txt
protocol/scripts/       the pipeline as run
protocol/envs/          conda environment exports and tool commits
protocol/inputs/        receptor structures, RFdiffusion3 contig definitions, MSA note
analysis/               this document's supporting analyses:
                          replicate_consistency.md  the three-replicate result
                          controls.md               no-receptor and self-pairing controls
                          composition.md            sequence composition, and its limits
```

---

## 3. Results

All ten designs are **140 aa**, **single-chain**, and **cysteine-free** (0 cysteines across all
ten sequences). All ten anchors satisfy every position restraint of the agonist consensus
(section 4); none is identical to native GLP-1 `HAEGTFTSDV`, and eight distinct anchors appear
across the ten.

Figures from `SUBMISSION.csv` / `metrics/seed_consistency.csv`, across the ten designs:

| | across the ten | median |
|---|---|---|
| mean ipSAE (of 3 replicates) | 0.8832 - 0.9180 | **0.9018** |
| ipSAE spread across replicates | 0.0025 - 0.0205 | 0.0124 |
| mean binder pLDDT | 0.8820 - 0.9172 | **0.9062** |
| mean anchor RMSD vs the AF3 peptide | 0.649 - 1.849 Å | **0.772 Å** |
| anchor RMSD spread across replicates | 0.245 - 1.436 Å | 0.458 Å |
| replicates holding the acceptance gate | **9 of 10 hold all three** | |
| replicates holding the full hotspot set | **9 of 10 hold all three** | |

Pooled over all **30 predictions** (ten designs × three replicates):

| | range | median |
|---|---|---|
| ipSAE | 0.881 - 0.924 | 0.902 |
| binder pLDDT | 0.880 - 0.918 | 0.906 |
| ipTM | 0.899 - 0.950 | 0.938 |
| anchor RMSD | 0.410 - 2.398 Å | 0.837 Å |
| receptor RMSD | 2.901 - 4.758 Å | 3.138 Å |
| interface salt bridges | 5 - 10 | 8 |
| interface hydrophobic contacts | 64 - 110 | 83.5 |

**The acceptance gate** is: anchor RMSD < 2.0 Å **and** ≥ 3 hotspots contacted **and** receptor
RMSD < 5.0 Å, all measured on an unrestrained co-fold. The two designs that do not hold it, or
do not hold the full hotspot set, in every replicate are named rather than hidden:

- **`glp1r-08` (rank 8, parent `p12_41`)** holds the gate in **1 of 3** replicates. Its anchor
  RMSD is 0.962 Å in replicate 1 but 2.398 Å and 2.188 Å in replicates 2 and 3 — over the 2.0 Å
  line in both. It keeps all three of its hotspots in every replicate.
- **`glp1r-05` (rank 5, parent `p13_23`)** holds the gate in all three but drops from 4 hotspots
  to 3 in replicate 2.

Both are submitted. The case for keeping them, and the argument that a single prediction
could not have told them apart from the rest, is in `analysis/replicate_consistency.md`.

### Ranking basis

`rank` is carried over from the p15 composite score — the mean of three min-max normalised
terms over the gated set (anchor RMSD inverted, ipSAE at 15/15 cutoffs, binder-chain mean
pLDDT), with hotspot count ordered first. It is **a preference order, not a prediction**, for a
reason measured on this target and stated in section 5: ipSAE barely separates on-target from
off-target here.

---

## 4. How the campaign narrowed

Seven design rounds, then two replicate rounds. Refold and gate counts are from each round's own
`p*_best/README.txt` in the campaign tree; the gate is the same three-part rule throughout.

| round | tag | epitope | gate rule | refolds | gated | rate |
|---|---|---|---|---|---|---|
| p2 | `p2_1088883` | A24/A25/A27/A42 | see note | 1200 | 869 | 72% |
| p3 | `p3_1088996` | A86-A89 | ≥ 3 of 4 | 480 | 309 | 64% |
| p10 | `p10_1089071` | A27/A42/A86 | ≥ 3 of 3 | 90 | 64 | 71% |
| p12 | `p12_1089085` | A27/A42/A86 | ≥ 3 of 3 | 1197 | 779 | 65% |
| p13 | `p13_1089119` | A27/A42/A86/A95 | ≥ 3 of 4 | 1200 | 822 | 68% |
| p14 | `p14_1092297` | per parent | ≥ 3 of 4 | 600 | 432 | 72% |
| p15 | `p15_1099477` | per parent | ≥ 3 of 4 | 210 | **187** | **89%** |
| replicate 2 | `p18s2_1106925` | per parent | ≥ 3 of 4 | 150 | 138 | 92% |
| replicate 3 | `p18s3_1106926` | per parent | ≥ 3 of 4 | 150 | 135 | 90% |

(Replicate gate counts recomputed from `metrics/seed_consistency.csv` on the same rule; the
column `n_seeds_gate` reproduces exactly on all 150 rows.)

**p2 died.** The A24/A25/A27/A42 epitope was the first one tried and it is the only lineage
with no descendant anywhere in the final set: the 150 designs carried into the replicate rounds
are p13 (55), p12 (42), p10 (38) and p3 (15), and the ten submitted designs are p13 (3),
p10 (3), p3 (2) and p12 (2). p2 is also the only round whose `README.txt` is internally
inconsistent — the gate line reads "≥ 3 of 3 hotspots" while its own hotspot histogram runs to
a 4 bin (4:774 3:202 2:86 1:61 0:77). Nothing in this submission depends on that round, so it
is recorded rather than resolved.

**p14 inverted the charge, and broke the anchor.** p14 rebuilt the p3/p10/p12/p13 bodies with
the designed-region charge flipped and the anchor left free to ProteinMPNN. The charge
inversion is unambiguous: the formal net charge of residues 11-140 across each round's 50 best
designs moves from a median of −8 (p2), −5.5 (p3), −5 (p10, p13) and −4 (p12) to **+10 (p14)**,
and every one of p14's 50 best carries a positive body against 0-5 of 50 in each earlier round.
The interaction gain was real but small — designed-region salt bridges have median 9 across the
142 rows of `metrics/p14_band_q4to9.csv` against `sb_designed` medians of 7.0 (p12, p13),
7.5 (p10) and 8.0 (p3) in the per-round baseline tables, which are not a controlled comparison.

The cost was not small. Running `protocol/scripts/consensus.py` over the round's own interaction tables:
of **600 p14 anchors (233 distinct), 0 satisfy the agonist consensus, and 0 satisfy even its
five invariant positions.** Position 4 (invariant G) is wrong in 600/600 and position 5
(invariant T) in 600/600; position 9 (invariant D) in 568/600 and position 8 (invariant S) in
548/600. The most common p14 anchor is `HAEEEFFDEV` — `HAEGTFTSDV` with the agonist's
`GT..SD` core overwritten.

**p15 repaired the anchor and kept the body.** p15 re-imposed the consensus on positions 1-10
while keeping the charge-inverted p14 body: all ten submitted sequences have residues 11-140
**byte-identical** to a p14 row (every one appears in `metrics/p14_band_q4to9.csv`), and
50 of 50 p15 best anchors satisfy every consensus restraint against 0 of 50 for p14. p15 has
the **highest gate rate of any round, 187/210 = 89%**, against 64-72% for every other round.
Re-anchoring under the consensus did not cost epitope engagement; it bought it.

**The funnel, end to end.** 4977 refolds across the seven design rounds → 187 gated p15 refolds
(`metrics/p15_RANKED.csv`) → the **top 150 by composite score** carried into two further
replicates (verified: the 150 rows of `seed_consistency.csv` are exactly the top 150 of the 187
by `score`, cutting at 0.4897) → ten submitted.

### The consensus restraints

`protocol/scripts/consensus.py` holds the position restraints, derived from a 23-agonist alignment
(frame: binder position 1 = GLP-1 His7).

| binder pos | GLP-1 pos | allowed | note from the script |
|---|---|---|---|
| 1 | 7 | H, Y, F | variable, 82.6% H |
| 2 | 8 | G, A, S | variable; Aib leads at 34.8% but is non-canonical |
| 3 | 9 | E, D | acid class 65.2%; Q excluded to keep the acid requirement |
| 4 | 10 | **G** | INVARIANT 100% |
| 5 | 11 | **T** | INVARIANT 100% |
| 6 | 12 | **F** | INVARIANT 100% |
| 7 | 13 | T, I | near-invariant, 95.7% T |
| 8 | 14 | **S** | INVARIANT 100% |
| 9 | 15 | **D** | INVARIANT 100% |
| 10 | 16 | Y, V, L, K | variable; hydrophobic class 56.5% |

Native `HAEGTFTSDV` satisfies every restraint, and so do all ten submitted anchors.

---

## 5. Honest limitations

These change how the ranking should be read.

- **There is no experimental data here at all.** Every number in this repository is a
  structure-prediction output or a sequence property computed from one. No binding assay, no
  expression measurement, no biophysical characterisation. The EGFR submission in the sibling
  repository could at least calibrate against an experimental benchmark it had assembled; this
  one cannot.

- **ipSAE barely detects epitope adherence on this target, and must not be read as evidence that
  a design is on-target.** Measured on the campaign's own folds and recorded in the round
  READMEs: median ipSAE **0.882 at 4/4 hotspots versus 0.871 at 0 hotspots** over 1200 p2 folds
  — a separation of 0.011 — and **0.894 versus 0.870** over 480 p3 folds, a separation of 0.024.
  The geometric per-hotspot check in step 1 of each round's gate, not the score in step 2, is
  what establishes epitope engagement. The ipSAE column is a tie-breaker among designs the
  geometry has already cleared.

- **ipSAE is also the least reproducible metric here.** Between the three replicates of the same
  150 inputs, ipSAE correlates at Pearson r **+0.309 to +0.450**, while binder pLDDT correlates
  at **+0.945 to +0.954** and ipTM at +0.795 to +0.859. Ranking on a metric that reproduces at
  r ≈ 0.4 across a reseed is ranking partly on noise. Details in
  `analysis/replicate_consistency.md`.

- **A single prediction could not have identified the unstable designs.** Across the 150, the
  replicate-1 ipSAE separates designs that hold the gate in all three replicates from those that
  do not at **AUC 0.487** — indistinguishable from chance — and binder pLDDT at 0.468. Only the
  two RMSD terms carry any signal at all (anchor RMSD AUC 0.662, receptor RMSD 0.643, both in
  the "lower is better" direction), and neither is strong enough to act on. This is the main
  argument for having run the replicates.

- **The replicate figures are read off a set that was selected on replicate 1, so 132/150 is an
  upper bound.** All 150 pass the gate in replicate 1 **by construction** — they were chosen
  from the p15 gated set. The honest statement is that 18 of 150 designs that a first prediction
  passed were not confirmed by two more; it is not an unbiased reproducibility rate.

- **About 44% of the predicted receptor interface is outside the ECD.** Across the 250 designs in
  `metrics/interface_properties.csv`, the median fraction of receptor interface residues lying
  within the ECD (1-117) is **0.562** (range 0.487-0.636); 246 of 250 are at or above half. The
  binders were designed for the ECD face and the hotspots are all in it, but in the predicted
  complexes a substantial minority of contacts fall in 118-409. That table covers only the
  pre-inversion rounds (p2/p3/p10/p12/p13) — it was not recomputed for p14, p15 or the ten.

- **Sequence composition sits far outside the only quantitative envelope available, which is
  from a different target.** All ten designs exceed the alanine ceiling of the EGFR
  composition envelope (median alanine **0.200** against a ceiling of 0.104, 1.92×), and 8 of 10
  fall outside its hydrophobic band. That envelope's own documentation states it should not be
  imported to another target, and this document respects that — the comparison and the caveat
  are in `analysis/composition.md`, which declines to call it a verdict.

- **Whole-chain RMSD is the wrong measurement for this architecture, and it looks alarming if you
  use it.** Folded with no receptor, the 210 binders in `metrics/p16_monomer_vs_bound.csv` move
  by a median of **10.283 Å** over residues 1-140, which would read as a failed design set. Over
  residues 50-120 the same monomers move by a median of **0.243 Å**, with 98.6% under 1.0 Å. The
  difference is the 10-residue agonist arm, which has nothing to hold it when the receptor is
  absent. Which window is honest, and why, is set out in `analysis/controls.md`.

- **Self-association is low but not zero.** Folded against a copy of itself, the 210 binders give
  a median self ipTM of 0.262, but 24 of 210 reach ≥ 0.5 and one reaches 0.810. The ten
  submitted designs run 0.123-0.467, all below 0.5 — but none of them was *selected* on this
  axis, and nothing here establishes that a Boltz self-fold predicts real oligomerisation.
  See `analysis/controls.md`.

- **Campaign-wide amino-acid composition cannot be computed from what is shipped.**
  `metrics/campaign_seq_properties.csv` covers 4167 designs but carries no sequence and no
  per-amino-acid column, so alanine content across the whole campaign is not recoverable. Every
  alanine figure in this repository comes from the 350 sequences in the seven rounds' best-50
  tables, the 142 in `p14_band_q4to9.csv`, or the ten submitted. Note also that the
  `f_hydrophobic` column in that table uses the alphabet **AVILMFWC** — verified against
  `receptor_properties.csv` — which is not the AVILMFWY used for the envelope comparison; the two
  are not interchangeable.

- **Three predictions of one input is a thin replicate design.** The replicates differ only in
  Boltz seed. They do not vary the MSA, the model, the receptor construct or the method, so they
  measure sampling variance within one predictor and nothing about whether that predictor is
  right. A cross-method or cross-species arm would have tested something this one cannot.

- **The campaign never left 140 residues.** Every design in every round is exactly 140 aa. No
  length exploration was run, so nothing here says whether a longer or shorter binder would do
  better on this target.
