# Sequence composition, and the limits of the only envelope available

This campaign never had a composition filter. Nothing in any round's gate looked at amino-acid
content, so what follows is a description of what the pipeline produced, not an audit against a
standard it was built to meet.

There is one quantitative envelope available to compare against, and it is from a **different
target**. The sibling repository `submission/egfr-ph-selective-binders/` records it in
`analysis/composition_envelope.md`. That document is explicit that its bounds should not be
carried across targets. This document therefore reports the comparison and declines to call it
a verdict. Section 3 explains why.

Sources: sequences from `SUBMISSION.csv` and the seven rounds' `p*_best/metrics.csv` in the campaign tree (not shipped here);
`metrics/p14_band_q4to9.csv`; `metrics/campaign_seq_properties.csv`;
`metrics/receptor_properties.csv`; `metrics/interface_properties.csv`.

---

## 1. What the ten submitted designs are made of

All ten are 140 aa. Alanine fraction and the AVILMFWY hydrophobic fraction are counted directly
from the sequences; formal net charge is (K+R) − (D+E).

| rank | parent | %Ala | AVILMFWY | formal net charge | charged D+E+K+R | %Gly | %Cys |
|---|---|---|---|---|---|---|---|
| 1 | p13_04 | 0.200 | 0.493 | +6 | 0.314 | 0.079 | 0 |
| 2 | p3_20 | **0.329** | 0.543 | +6 | 0.314 | 0.057 | 0 |
| 3 | p13_06 | 0.214 | 0.521 | +4 | 0.271 | 0.064 | 0 |
| 4 | p3_19 | 0.200 | **0.386** | +5 | 0.450 | 0.057 | 0 |
| 5 | p13_23 | 0.179 | 0.464 | +6 | 0.329 | 0.036 | 0 |
| 6 | p10_27 | **0.336** | 0.543 | +5 | 0.321 | 0.050 | 0 |
| 7 | p10_05 | 0.179 | **0.407** | +6 | 0.414 | 0.036 | 0 |
| 8 | p12_41 | 0.186 | 0.464 | +4 | 0.357 | 0.057 | 0 |
| 9 | p10_22 | 0.207 | 0.493 | +5 | 0.264 | 0.064 | 0 |
| 10 | p12_29 | 0.200 | 0.500 | +4 | 0.300 | 0.029 | 0 |
| | **median** | **0.200** | **0.493** | **+5** | 0.318 | 0.057 | 0 |

Two unambiguously good properties: **no cysteine anywhere in the ten**, so there is no free-thiol
or disulfide-scrambling liability to manage; and glycine stays between 2.9% and 7.9%.

The set is **alanine-rich and hydrophobic-rich**, and two designs are extreme on both axes
together: rank 2 (`p3_20`) at 32.9% alanine and 54.3% AVILMFWY, and rank 6 (`p10_27`) at 33.6%
and 54.3%. One third of a designed binder being alanine is a lot of sequence doing very little.

The set is also **net positive**, +4 to +6 on every design, which is a deliberate inversion and
is discussed in section 4.

---

## 2. The campaign around them

The ten are not outliers within this campaign. They are typical of it.

| set | n | %Ala median (range) | AVILMFWY median (range) | cysteine-free |
|---|---|---|---|---|
| the ten submitted | 10 | 0.200 (0.179 - 0.336) | 0.493 (0.386 - 0.543) | 10/10 |
| seven rounds' best-50, pooled | 350 | 0.221 (0.064 - 0.450) | 0.471 (0.329 - 0.657) | 350/350 |
| p14 charge band | 142 | 0.214 (0.064 - 0.429) | 0.486 (0.393 - 0.636) | — |

Round by round, from the 50 best of each:

| round | %Ala median | %Ala max | AVILMFWY median | AVILMFWY max |
|---|---|---|---|---|
| p2 | 0.214 | 0.357 | 0.457 | 0.571 |
| p3 | 0.207 | 0.386 | 0.443 | 0.586 |
| p10 | 0.243 | 0.421 | 0.482 | 0.643 |
| p12 | 0.243 | 0.450 | 0.489 | 0.657 |
| p13 | 0.214 | 0.371 | 0.471 | 0.600 |
| p14 | 0.214 | 0.314 | 0.482 | 0.564 |
| p15 | 0.214 | 0.436 | 0.486 | 0.650 |

Composition is flat across seven rounds. It never moved because nothing ever pushed on it —
ProteinMPNN was given no alanine bias and no composition term entered any gate. The alanine
content of this campaign is a default of the sampler, not a choice anyone made.

**Every design in the campaign is 140 aa**, so length cannot be compared to anything.

---

## 3. The EGFR envelope, the caveat attached to it, and the comparison

### The caveat, from the document that defines the bounds

`submission/egfr-ph-selective-binders/analysis/composition_envelope.md` states where its four
bounds come from and what they are not. Quoting it:

> They are the **observed range of a small set of all-alpha EGFR binders from the previous
> Adaptyv competition round that were experimentally confirmed sub-micromolar.** Each bound is
> the minimum or maximum actually seen in that set …

> - **Not a fitted model.** No regression, no cross-validation, no held-out test. It is a
>   bounding box around a handful of successes.
> - **Not a physical limit.** A design outside it is not predicted to fail; it is outside the
>   range where this particular class of binder has been observed to succeed.
> - **Not derived from a large n.** The confirmed sub-micromolar all-alpha set is on the order
>   of ten designs. Four bounds from ten observations will be tight in some directions by chance.

And the companion document `analysis/what_predicts_binding.md`, on the same repository's own
3027-design / 32-target experimental data, is blunter still:

> Directions reverse between targets — binders are significantly *less* hydrophobic on
> DERF7/FGF-R1/MDM2 and significantly *more* on MZB1/PD-L1. **Composition filters must be derived
> from the target's own cohort, never imported.**

**GLP-1R is a different target, in a different structural class** — a class B GPCR engaged at an
acidic extracellular domain, not an EGFR domain III ligand pocket — and **no GLP-1R binder
cohort of experimentally confirmed binders exists anywhere in this repository** to derive a
target-specific envelope from. So the only quantitative bounds available are precisely the ones
their own author says not to import.

The comparison below is therefore run and reported, because a reader will otherwise run it
themselves and reach a worse-informed conclusion. It is **not** a pass/fail judgement on these
designs, and the word "violation" is avoided deliberately.

### The comparison

| axis | EGFR bound | EGFR centroid | the ten submitted | outside |
|---|---|---|---|---|
| alanine fraction | ≤ 0.104 | 0.085 | median **0.200** (0.179 - 0.336) | **10 / 10** |
| hydrophobic AVILMFWY | 0.329 - 0.42 | 0.376 | median **0.493** (0.386 - 0.543) | **8 / 10** |
| net charge at pH 7 | −19.3 … −5.0 | — | **not computable for these ten** — see below | — |
| pI | 4.60 - 5.00 | — | **not computable for these ten** — see below | — |

On alanine the gap is large and uniform: the median is **1.92× the ceiling** and **2.35× the
centroid**, and no design is close — the lowest is 0.179, still 1.7× the ceiling. Across the 350
round-best sequences, **346 of 350** are above 0.104.

On hydrophobicity, 8 of 10 sit above the 0.42 upper bound; the two inside are rank 4
(`p3_19`, 0.386) and rank 7 (`p10_05`, 0.407), both near the centroid. Across the 350 round-best
sequences only **41 of 350** fall inside 0.329 - 0.42.

If the EGFR envelope were imported wholesale and applied as a filter, it would reject this
entire campaign. That is itself a reason to doubt its transferability rather than the campaign,
but it is also exactly what a reader should know before weighting the comparison either way.

### Why the charge and pI axes cannot be evaluated for the ten

`metrics/campaign_seq_properties.csv` carries pI and the Henderson–Hasselbalch charge at pH 7
for 4167 designs — the complete refold sets of p2 (1200), p3 (480), p10 (90), p12 (1197) and
p13 (1200). It does **not** cover p14, p15, or the ten submitted sequences, and it carries no
sequence column, so the ten cannot be joined into it.

Substituting each design's backbone parent was tested and rejected: the ten sequences differ
from their same-named p-round parents at **43 to 64 of the 130 designed positions**, so the
parents' pI and pH-7 charge are not the designs' own. No number is reported for those two axes
rather than a near-enough one.

What **can** be said, at campaign level, from the 4167 designs that are covered:

| axis | campaign median | inside the EGFR bound |
|---|---|---|
| net charge at pH 7 | −8.02 (p10 −13.99, p90 −3.02) | 3199 / 4167 (**76.8%**) |
| pI | 4.95 (p10 4.62, p90 5.58) | 2032 / 4167 (**48.8%**) |
| both together | — | 2019 / 4167 (48.5%) |

Note the irony: on the two axes that *can* be checked campaign-wide, the pre-inversion designs
sit largely **inside** the EGFR envelope, and on the two that can be counted from sequence they
sit far outside it. The envelope is not pointing in one consistent direction here.

### One cross-table trap

The `f_hydrophobic` column in `campaign_seq_properties.csv` uses the alphabet **AVILMFWC** —
verified by reproducing the receptor's reported 0.4743 exactly from its sequence, where AVILMFWY
gives 0.4890. It includes cysteine and **excludes tyrosine**. The EGFR envelope is defined on
**AVILMFWY**. The two columns are not interchangeable, and every AVILMFWY figure in this
document is counted from sequence rather than read from that column.

---

## 4. The one composition decision anyone actually made: charge inversion

Composition drifted by default in this campaign with a single exception — the p14 round
deliberately flipped the sign of the designed region's charge, and the ten submitted designs
carry that inversion.

Formal net charge, (K+R) − (D+E), from the 50 best of each round:

| round | whole chain, median | designed 11-140, median | designs with a positive body |
|---|---|---|---|
| p2 | −10 | −8 | 0 / 50 |
| p3 | −7.5 | −5.5 | 2 / 50 |
| p10 | −7 | −5 | 3 / 50 |
| p12 | −6 | −4 | 3 / 50 |
| p13 | −7 | −5 | 5 / 50 |
| **p14** | **+6** | **+10** | **50 / 50** |
| p15 | +6 | +8 | 50 / 50 |
| the ten | +5 | +7 | 10 / 10 |

(The `q_des` column of `p14_band_q4to9.csv` is exactly this quantity: it equals the computed
formal charge of residues 11-140 on **142 of 142** rows.)

**The target supports the direction.** From `metrics/receptor_properties.csv`, the face being
designed against is acidic and the rest of the receptor is not:

| region | length | pI | charge at pH 7 |
|---|---|---|---|
| receptor, full 1-409 | 409 | 7.93 | +2.45 |
| **ECD 1-117 — the designable face** | 117 | **4.88** | **−5.18** |
| TMD 118-409 | 292 | 8.96 | +7.39 |
| epitope A86-A89 ±5 | 14 | 4.37 | −1.25 |
| epitope A27/A42/A86 ±5 | 33 | 4.43 | −3.26 |
| epitope A27/A42/A86/A95 ±5 | 42 | 4.39 | −3.26 |
| epitope A24/A25/A27/A42 ±5 | 25 | 6.20 | −0.25 |

Every epitope neighbourhood the campaign targeted is acidic or neutral, and the ECD as a whole
carries −5.18 at pH 7. A positively charged binder is **charge-complementary to the face it is
meant to bind**, which is a real argument for the inversion and is independent of the EGFR
envelope — whose net-charge window (−19.3 … −5.0) was derived against a different receptor
surface entirely, and would have pushed these designs in the wrong direction.

**What the inversion bought was modest.** Designed-region salt bridges have median 9 across the
142 rows of `p14_band_q4to9.csv` (9 across the 107 of them that pass the gate), against
`sb_designed` medians of 7.0 (p12, p13), 7.5 (p10) and 8.0 (p3) in the pre-inversion baseline
tables `analysis/physchem/baseline_p*_best_interactions.csv`, which live in the campaign tree and are not shipped here. That is a shift of one to two
bridges, between two tables built on different design sets under different selection — it is not
a controlled comparison and should not be read as a measured effect.

Two further reasons to hold it loosely. The band table's `sb_des` is taken to mean the same
quantity as the baseline tables' `sb_designed`, which is an inference from column naming: no
file states it. In the baseline tables `salt_bridges` = `sb_designed` + `sb_anchor` exactly on
all 200 rows, whereas in the band table `sb_des` **exceeds** `sb` on 80 of 142 rows, so `sb`
there must be the deduplicated count (the baselines' `unique_salt_bridges`, median 6.0-6.5)
rather than the total. The two tables are being read against each other on a column-name match,
not a documented definition.

**What it cost, before being repaired, was not modest.** p14 also let ProteinMPNN redesign the
10-residue agonist anchor, and the result is in `README.md` section 4: 0 of 600 p14 anchors
satisfy the agonist consensus, including 600/600 wrong at the invariant G and T positions.
p15 re-imposed the consensus on positions 1-10 while keeping the charge-inverted body, which is
why the ten submitted designs combine a compliant anchor with a positive body.

---

## 5. Interface composition

`metrics/interface_properties.csv`, 250 designs — **the 50 best of each of p2, p3, p10, p12 and
p13 only**. It was never recomputed for p14, p15 or the ten, so it describes the campaign
*before* the charge inversion and cannot be used to characterise what is being submitted.

Medians across all 250:

| | f_hydrophobic | f_positive | f_negative | net charge | residues |
|---|---|---|---|---|---|
| binder interface | 0.378 | 0.156 | 0.163 | 0 | 50 |
| …of which the designed region | 0.400 | — | — | **+2** | 40 |
| …of which the anchor | — | — | — | **−2** | **10 in all 250** |
| receptor interface | 0.372 | — | — | +1 | 79 |

Three things worth recording.

**The anchor contributes all ten of its residues to the interface in all 250 designs**, and it
contributes a net −2 while the designed region contributes +2. The acidic anchor and the basic
designed surface are doing opposite electrostatic jobs, which is consistent with the anchor
being a fixed agonist motif and the designed surface being free to complement the acidic ECD
around it.

**Interface charge is not consistently signed even before the inversion**: 115 of 250 binder
interfaces are net negative, 46 neutral and 89 positive, against receptor interfaces that are
142 positive, 58 neutral and 50 negative. The campaign was not building a uniformly
charge-complementary interface at this stage; p14 made that explicit.

**About 44% of the predicted receptor interface is outside the ECD.** The median fraction of
receptor interface residues lying within residues 1-117 is **0.562** (range 0.487 - 0.636);
246 of 250 are at or above half. The hotspots are all in the ECD and the designs were built for
it, but in the predicted complexes a substantial minority of contacts land in 118-409. This is
recorded as a limitation in `README.md` section 5.

---

## 6. What this document does not establish

- **No composition figure here is tied to a binding outcome.** There is no experimental cohort
  for GLP-1R anywhere in this repository, so nothing can be said about which compositions bind
  this target. Every statement above is descriptive.
- **The EGFR envelope comparison is reported, not applied.** Its own documentation says the
  bounds are a bounding box around roughly ten successes on another target, are not a fitted
  model, and must not be imported. Nothing was filtered, ranked, penalised or excluded on it.
- **Campaign-wide alanine content is not recoverable.** `campaign_seq_properties.csv` has no
  sequence and no per-amino-acid column. Every alanine figure here comes from the 350 round-best
  sequences, the 142 in the p14 band, or the ten submitted — never from the full 4167.
- **pI and pH-7 titration charge are unavailable for the ten submitted designs** and are
  reported only at campaign level, for the five pre-inversion rounds the physchem table covers.
- **The 33% alanine of ranks 2 and 6 is flagged, not resolved.** It is outside anything observed
  in the one reference set available, and this campaign has no evidence either way about whether
  it matters on GLP-1R. If one axis of this submission deserves a re-design round, it is this one.
