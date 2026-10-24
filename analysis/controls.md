# Two negative controls: no receptor, and no partner but itself

Two questions a co-folding score cannot answer. **Does the binder hold its fold when the
receptor is not there?** and **does it prefer the receptor to another copy of itself?** Both were
run as separate Boltz-2 jobs.

Sources: `metrics/p16_monomer_vs_bound.csv` and `metrics/p16_monomer_windows.csv` (210 rows
each, binder folded alone, `protocol/scripts/p16_mono.sbatch`); `metrics/p17_self_association.csv`
(210 rows, binder folded against a copy of itself, `protocol/scripts/p17_dimer.sbatch`).

**What the 210 are.** Both tables cover the complete p15 refold set — the 187 that passed the
p15 gate plus the 23 that did not, verified by name against `metrics/p15_RANKED.csv`. All 210
sequences are distinct (`n_sharing_seq` = 1 on every row), but they collapse to **71 distinct
name prefixes** — three sequences per prefix in 68 cases, two in three — so these are 210
sequences over 71 upstream design entities, not 210 independent backbones. All ten submitted
designs are in both tables.

---

## 1. The no-receptor control, and why whole-chain RMSD is the wrong ruler here

### The number you get if you measure the whole chain

Fold each binder with no receptor present and superpose it on its bound pose:

| window | n | min | median | p90 | max | < 1.0 Å | < 2.0 Å |
|---|---|---|---|---|---|---|---|
| **all, 1-140** | 210 | 0.718 | **10.283 Å** | 17.437 | 21.739 | 10 (4.8%) | 46 (21.9%) |
| designed, 11-140 | 210 | 0.306 | **5.811 Å** | 11.738 | 13.216 | 40 (19.0%) | 73 (34.8%) |
| 30-140 | 210 | 0.144 | 0.424 Å | 1.681 | 4.928 | 174 (82.9%) | 194 (92.4%) |
| 40-140 | 210 | 0.131 | **0.294 Å** | 0.718 | 4.726 | 200 (95.2%) | 207 (98.6%) |
| core, 50-120 | 210 | 0.111 | **0.243 Å** | 0.451 | 5.163 | 207 (98.6%) | 207 (98.6%) |

Read only the first row and this looks like a failed design set: a median 10.3 Å whole-chain
RMSD would normally mean the fold does not survive without its partner. Read the last row and
the same 210 structures are essentially rigid — median 0.243 Å, 98.6% under 1.0 Å.

### Why the whole-chain number is misleading for this architecture

The binder is 140 residues in two parts that are not alike. Residues **1-10 are the GLP-1
agonist N-terminus**, an extended peptide arm whose entire job is to lie in the receptor's
orthosteric groove. It is not packed against the binder's own body. Remove the receptor and
there is nothing to hold it, so it flops — and because RMSD is computed after a global
superposition, ten badly placed residues at one end of a 140-residue chain lever the whole
alignment and inflate the number for the other 130.

**The structures say this directly.** In the monomer folds the anchor is predicted with low
confidence and the designed region with high confidence, on the same chains:

| | min | median | max |
|---|---|---|---|
| monomer pLDDT, whole chain | 0.703 | 0.880 | 0.942 |
| monomer pLDDT, **anchor 1-10** | 0.410 | **0.584** | 0.812 |
| monomer pLDDT, **designed 11-140** | 0.722 | **0.901** | 0.956 |

A median gap of **0.318 pLDDT** between the two halves of the same chain, with the anchor below
0.70 in 182 of 210 monomers and the designed region at or above 0.85 in 189 of 210. The model
is telling you which part it cannot place.

**And the whole-chain number is not reporting a body that moved.** Of the 145 monomers with
whole-chain RMSD ≥ 3.0 Å, the core 50-120 window has median **0.250 Å** and 133 of the 145 are
under 0.5 Å. The large whole-chain values are the arm, essentially always.

### Which window is honest

**Residues 40-140, or 50-120.** Both are stable across every lineage and both exclude the
agonist arm and the short linking stretch behind it. 50-120 is the stricter core and 40-140 is
the more complete body; they agree (median 0.294 Å vs 0.243 Å).

Residues 11-140 — the "designed region" — is **not** an honest window despite its name, at a
median 5.811 Å. Excluding the anchor from the *residue set* does not undo its effect on the
*superposition*, and the region immediately behind the anchor moves with it: the window boundary
has to retreat to about residue 30 before the number settles (30-140: median 0.424 Å, 82.9%
under 1.0 Å; 40-140: median 0.294 Å, 95.2% under 1.0 Å).

Per lineage the contrast is the same everywhere, and the whole-chain figure mostly reflects how
far the arm happened to swing:

| lineage | n | median whole-chain | median core 50-120 |
|---|---|---|---|
| p10 | 45 | 5.273 Å | 0.231 Å |
| p13 | 65 | 8.474 Å | 0.254 Å |
| p12 | 58 | 10.185 Å | 0.233 Å |
| p3 | 42 | 15.592 Å | 0.270 Å |

### What actually fails

**Three of 210 monomers fail the core window** (core 50-120 ≥ 1.0 Å, at 5.045, 5.133 and
5.163 Å). All three share the prefix `p3_10_sh7_m1_0_model_1__s0__s2` — they are **three
sequences of one design**, not three designs. That design also has the lowest monomer pLDDT in
the set (0.703-0.727) and an anchor pLDDT of 0.437-0.460. It is not among the ten submitted.

On the looser 30-140 window, 36 of 210 exceed 1.0 Å and 16 exceed 2.0 Å — so the stability
claim is genuinely window-dependent and is stated here with the window attached.

### The ten submitted designs

| rank | parent | 1-140 | 11-140 | 30-140 | 40-140 | **50-120** | monomer pLDDT | anchor pLDDT |
|---|---|---|---|---|---|---|---|---|
| 1 | p13_04 | 15.044 | 7.783 | 0.353 | 0.357 | **0.345** | 0.854 | 0.541 |
| 2 | p3_20 | 1.289 | 0.306 | 0.224 | 0.222 | **0.235** | 0.878 | 0.650 |
| 3 | p13_06 | 15.693 | 8.796 | 0.269 | 0.271 | **0.284** | 0.870 | 0.603 |
| 4 | p3_19 | 17.583 | 12.472 | 2.128 | 0.477 | **0.342** | 0.837 | 0.531 |
| 5 | p13_23 | 15.329 | 9.826 | 0.640 | 0.271 | **0.256** | 0.865 | 0.517 |
| 6 | p10_27 | 3.563 | 1.503 | 0.571 | 0.455 | **0.329** | 0.903 | 0.621 |
| 7 | p10_05 | 1.638 | 0.920 | 0.457 | 0.441 | **0.242** | 0.933 | 0.793 |
| 8 | p12_41 | 13.598 | 9.958 | 0.628 | 0.220 | **0.214** | 0.860 | 0.456 |
| 9 | p10_22 | 15.361 | 10.003 | 0.681 | 0.454 | **0.392** | 0.861 | 0.536 |
| 10 | p12_29 | 4.324 | 1.919 | 0.209 | 0.196 | **0.189** | 0.860 | 0.518 |
| | median | 14.321 | 8.290 | 0.514 | 0.314 | **0.270** | 0.863 | 0.538 |

All ten hold the core to under 0.4 Å. Rank 4 (`p3_19`) is the one design whose 30-140 window
exceeds 2.0 Å while its 40-140 window is 0.477 Å — a clean illustration of the boundary effect:
residues 30-39 move with the arm, residues 40 onward do not.

### What this control does *not* establish

It is a second structure prediction, not a measurement. It shows that Boltz-2 places the
designed body the same way with and without the receptor, which is consistent with a
self-supporting fold and inconsistent with a body that only exists because the receptor is
holding it. It is not evidence of thermodynamic stability, of folding in solution, or of
expression. And it says nothing about the anchor, which this control cannot assess by
construction — the arm's conformation is defined by the thing that was removed.

---

## 2. The self-pairing control

Each binder folded against a copy of itself. If a designed interface is generically sticky
rather than receptor-specific, this is where it shows.

| | min | p25 | median | p75 | p90 | max |
|---|---|---|---|---|---|---|
| self ipTM | 0.096 | 0.188 | **0.262** | 0.386 | 0.507 | **0.810** |
| pair ipTM A-B | 0.096 | 0.186 | 0.257 | 0.382 | 0.504 | 0.810 |
| dimer pLDDT | 0.639 | 0.783 | 0.810 | 0.837 | 0.861 | 0.919 |
| interface pLDDT | 0.545 | 0.668 | 0.709 | 0.743 | 0.790 | 0.885 |
| pTM | 0.411 | 0.503 | 0.533 | 0.603 | 0.662 | 0.874 |
| contacts A-B | 147 | 538 | 732 | 888 | 1114 | 1693 |
| interface residues on A | 13 | 29 | 35 | 44 | 51 | 73 |

### The comparison that matters

| | n | median ipTM |
|---|---|---|
| **binder + receptor** (`seed_consistency.csv`, 450 predictions) | 450 | **0.926** (min 0.882, max 0.961) |
| **binder + itself** (`p17_self_association.csv`) | 210 | **0.262** (min 0.096, max 0.810) |

The receptor complex scores 0.926 at its *worst* prediction of 450. The self-dimer scores 0.262
at its median and reaches 0.882 — the receptor complex's floor — in none of the 210. These
binders are predicted to prefer the receptor to themselves by a wide margin.

### But it is not zero, and the tail is real

| self ipTM threshold | designs |
|---|---|
| ≥ 0.3 | 92 / 210 (43.8%) |
| ≥ 0.4 | 50 / 210 (23.8%) |
| ≥ 0.5 | 24 / 210 (11.4%) |
| ≥ 0.6 | 10 / 210 (4.8%) |
| ≥ 0.7 | 4 / 210 (1.9%) |
| ≥ 0.8 | 1 / 210 (0.5%) |

The strongest is `p10_15_sh1_m1_0_model_0__s0__s1_model_0__s0` at **0.810**, with 877 A-B
contacts and an interface pLDDT of 0.730 — a confidently predicted homodimer. Two of the top
four self-associators are the same p10_15 prefix at different sequences, so the behaviour
tracks the body rather than being seed noise. By lineage the medians are p10 0.334, p12 0.259,
p13 0.247, p3 0.243; p10 also contributes 5 of the 10 designs above 0.6.

**The ten submitted designs all sit below 0.5** (range 0.123 - 0.467, median 0.200), with
rank 5 (`p13_23`, 0.467) and rank 7 (`p10_05`, 0.428) highest. That is a favourable outcome but
**not a selected one**: self-association was never an acceptance criterion in any round, so the
ten are below 0.5 by the same process that left 24 of 210 above it. A reader should treat this
as a check that passed, not as a filter that was applied.

### Self-association does not track anything else that was measured

| correlation | Pearson r |
|---|---|
| monomer pLDDT vs self ipTM | **+0.004** |
| core 50-120 monomer RMSD vs self ipTM | **−0.036** |

Both are zero. A design that folds cleanly on its own is neither more nor less likely to be
predicted to dimerise. This is the useful part of the control: self-association is an
independent axis that none of the campaign's other metrics was reporting on, which is why it
needed its own job.

### What this control does *not* establish

A Boltz-2 self-fold is not an oligomerisation assay. The model will produce an interface for
almost any pair of chains, and the absolute ipTM scale for a homodimer is not calibrated against
anything here. What can be said is relative and internal: against the same predictor, the
receptor complex scores far higher than the self-pair for every one of these designs, and the
ten submitted designs are in the lower half of the self-pairing distribution. Whether any of
them actually oligomerises is unknown and untested.
