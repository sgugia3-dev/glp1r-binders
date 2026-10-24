# Three predictions of the same input: what held, what did not

150 designs, each predicted three times by Boltz-2 from **byte-identical input**, differing only
in random seed. Campaign runs `p15_1099477` (replicate 1), `p18s2_1106925` (replicate 2) and
`p18s3_1106926` (replicate 3). One row per design, columns prefixed `p15_`, `p18s2_`, `p18s3_`.

Source: `metrics/seed_consistency.csv`, 150 data rows, 32 columns. Every number below is
recomputed from that file; the stored `n_seeds_gate` and `n_seeds_full` columns reproduce
exactly on all 150 rows under the stated gate rule, so the file is internally consistent.

**What this measures:** sampling variance inside one predictor. Nothing more. The replicates do
not vary the model, the MSA, the receptor construct or the method, so they cannot say whether
Boltz-2 is right about these designs — only how stable its answer is.

---

## 1. Read the headline number with its selection effect first

**132 of 150 designs (88.0%) hold the acceptance gate in all three replicates.**

That figure is not an unbiased reproducibility rate, and should not be quoted as one. The 150
designs **were selected on replicate 1**: they are the top 150 by composite score of the 187
p15 refolds that passed the gate (verified against `metrics/p15_RANKED.csv` — the 150 names are
exactly the top 150 by `score`, cutting at 0.4897). So replicate 1 passes 150 of 150 **by
construction**.

The defensible statement is narrower and more useful: **of 150 designs that one prediction
passed, 18 were not confirmed by two more.** Equivalently, a single pass carries roughly a
one-in-eight chance of not replicating.

The two later replicates, which carry no such selection, pass at similar rates to each other:

| replicate | passes the gate | would have been a false positive |
|---|---|---|
| 1 (`p15`) | 150 / 150 — **by construction** | 18 |
| 2 (`p18s2`) | 138 / 150 (92.0%) | 6 |
| 3 (`p18s3`) | 135 / 150 (90.0%) | 3 |

("False positive" = passes that replicate but does not hold all three.)

---

## 2. Which criteria held, and which did not

The gate is three conditions on an unrestrained co-fold: **anchor RMSD < 2.0 Å**, **≥ 3 hotspots
contacted** (measured from the designed region with the 10-residue anchor excluded), and
**receptor RMSD < 5.0 Å**.

| criterion | holds in all 3 replicates | holds in ≥ 1 |
|---|---|---|
| receptor RMSD < 5.0 Å | **150 / 150** | 150 |
| ≥ 3 hotspots contacted | **148 / 150** | 150 |
| anchor RMSD < 2.0 Å | **133 / 150** | 150 |
| all three together (the gate) | **132 / 150** | 150 |
| the *full* hotspot set (3 of 3 or 4 of 4, per lineage) | **137 / 150** | 150 |

**One criterion does essentially all of the failing.** Across 450 predictions there are 27
replicate-level gate failures: **24 are the anchor RMSD and 3 are the hotspot count. Receptor
RMSD never fails once** — the worst receptor RMSD anywhere in the 450 is 4.857 Å, inside a 5.0 Å
gate. That criterion is not discriminating anything here and should not be read as having been
passed in a meaningful sense.

**And the anchor failures are all marginal, not collapses.** The 24 failing anchor RMSDs run
**2.011 Å to 2.606 Å** — every one within 0.61 Å of the threshold. Nine are below 2.10 Å.
Moving the line would move the result a long way:

| anchor RMSD threshold | designs holding the gate in all three |
|---|---|
| 2.00 Å (as gated) | 132 / 150 |
| 2.25 Å | 142 / 150 |
| 2.50 Å | 146 / 150 |
| 3.00 Å | 148 / 150 |

The worst-of-three anchor RMSD over all 150 designs has median 1.435 Å and maximum 2.606 Å.
**No design in this set loses its anchor.** What the 2.0 Å gate is separating is a tail of
designs whose anchor wobbles by a few tenths of an Ångström between seeds, not a population that
binds in one prediction and misses in another.

**The full-hotspot losses are different in kind and more evenly spread.** There are 18
replicate-level losses of the full set, distributed 6 / 7 / 5 across replicates 1, 2 and 3 —
replicate 1 is not better here, because the gate it was selected on requires ≥ 3 hotspots, not
all of them. Fifteen of the 18 are four-hotspot designs dropping to three; three are
three-hotspot designs dropping to two (which also fails the gate). By lineage: p13 (10), p3 (5),
p12 (3), p10 (0).

### Per lineage

| lineage | n | hotspots in its set | gate 3/3 | full set 3/3 |
|---|---|---|---|---|
| p10 | 38 | 3 | 34 (89.5%) | **38 (100%)** |
| p12 | 42 | 3 | 34 (81.0%) | 40 (95.2%) |
| p13 | 55 | 4 | **52 (94.5%)** | 47 (85.5%) |
| p3 | 15 | 4 | 12 (80.0%) | 12 (80.0%) |

The four-hotspot lineages hold their full set less often, which is what you would expect from
having one more thing to hold; p10 and p12 are not better designs for showing 3 of 3.

---

## 3. Spread per metric

Per design, range = max − min across the three replicates. Values pooled over all 450
predictions for the left-hand columns.

| metric | median value | range over 450 | median per-design spread | 90th pct | worst |
|---|---|---|---|---|---|
| binder pLDDT | 0.907 | 0.852 - 0.936 | **0.005** | 0.010 | 0.016 |
| ipTM | 0.926 | 0.882 - 0.961 | 0.009 | 0.020 | 0.041 |
| ipSAE | 0.893 | 0.850 - 0.930 | **0.016** | 0.028 | 0.058 |
| receptor RMSD | 3.284 Å | 1.783 - 4.857 Å | 0.385 Å | 1.064 Å | 1.777 Å |
| anchor RMSD | 1.204 Å | 0.367 - 2.606 Å | **0.478 Å** | 0.837 Å | 1.521 Å |
| hotspots contacted | 3 | 2 - 4 | **0** | 0 | 1 |
| interface salt bridges | 7 | 3 - 13 | 2 | 3 | 5 |
| interface hydrophobic contacts | 81 | 40 - 123 | 11 | 20 | 37 |

Two things stand out.

**The anchor RMSD spread is comparable to the whole gate margin.** Median spread 0.478 Å against
a median value of 1.204 Å and a 2.0 Å threshold — so a typical design sits about 0.8 Å from the
line with 0.5 Å of seed noise on it. That is exactly the geometry that produces 18 marginal
failures, and it is a property of the measurement, not of the designs.

**The interaction counts are the noisiest thing reported.** Salt bridges move by a median of 2
(of a median 7) and hydrophobic contacts by a median of 11 (of a median 81) between predictions
of identical input. The `seed1_sb` / `seed1_hc` columns in `SUBMISSION.csv` should be read as
one draw from a distribution, not as a count.

### Between-replicate correlation

Pearson *r* and Spearman ρ over the three replicate pairs, 150 designs:

| metric | Pearson r (min-max, mean) | Spearman ρ (mean) |
|---|---|---|
| binder pLDDT | +0.945 - +0.954 (**+0.949**) | +0.940 |
| ipTM | +0.795 - +0.859 (+0.830) | +0.831 |
| hydrophobic contacts | +0.715 - +0.763 (+0.732) | +0.709 |
| salt bridges | +0.570 - +0.612 (+0.589) | +0.598 |
| anchor RMSD | +0.443 - +0.622 (+0.523) | +0.513 |
| **ipSAE** | +0.309 - +0.450 (**+0.401**) | +0.399 |

**ipSAE is the least reproducible metric in the set, and it is the one the ranking leans on.**
Binder pLDDT reproduces at r ≈ 0.95; ipSAE at r ≈ 0.40. A re-seed reshuffles the ipSAE order
substantially. This compounds a separate, already-documented problem with ipSAE on this target —
it barely separates on-target from off-target at all (median 0.882 at 4/4 hotspots against 0.871
at 0 hotspots over 1200 folds; 0.894 against 0.870 over 480 folds, from the round READMEs). A
metric that neither discriminates the thing you care about nor reproduces across a reseed cannot
carry a ranking on its own, and in this campaign it is not asked to: the geometric hotspot check
is what establishes engagement, and ipSAE only orders designs that check has already cleared.

---

## 4. Could a single prediction have predicted which designs were unstable?

**No.** This is the clearest result in the document.

Taking replicate 1 alone and asking how well each of its metrics separates the 132 designs that
hold the gate in all three from the 18 that do not (AUC, oriented so that higher is better for
each metric's own direction; 0.5 is chance):

| replicate-1 metric | AUC | stable median | unstable median |
|---|---|---|---|
| ipSAE | **0.487** | 0.894 | 0.894 |
| binder pLDDT | **0.468** | 0.907 | 0.908 |
| hydrophobic contacts | 0.516 | 80.0 | 76.5 |
| salt bridges | 0.534 | 7.0 | 7.0 |
| ipTM | 0.548 | 0.927 | 0.927 |
| receptor RMSD (lower better) | 0.643 | 3.223 Å | 3.427 Å |
| anchor RMSD (lower better) | **0.662** | 1.136 Å | 1.371 Å |

The three confidence metrics are **at chance**. Their medians are identical to three decimal
places in two cases. A design that would fail a reseed looks exactly like one that would not.

The only terms with any signal are the two RMSDs, and the signal is weak (0.662, 0.643) and
goes in the obvious direction: a design already closer to the anchor-RMSD threshold in
replicate 1 is somewhat likelier to cross it in replicate 2. That is a restatement of the
threshold, not an independent predictor, and it is nowhere near strong enough to screen on —
the stable and unstable distributions overlap heavily (stable median 1.136 Å, unstable 1.371 Å,
against a 2.0 Å line).

**For the full hotspot set the picture is worse: nothing predicts it at all.** Separating the
137 designs that keep their full set in all three from the 13 that do not:

| replicate-1 metric | AUC |
|---|---|
| ipSAE | 0.469 |
| binder pLDDT | 0.477 |
| anchor RMSD (lower better) | 0.469 |
| replicate-1 hotspot count | 0.523 |
| ipTM | 0.595 |

Every one is at or below chance except ipTM at 0.595, which on 13 positives is not a result.

---

## 5. What this implies

1. **Run the replicates.** There is no cheaper substitute. Nothing available from one
   prediction — not ipSAE, not pLDDT, not ipTM, not contact counts — identifies the 12% of
   designs a reseed will not confirm. The only way to know is to fold again.

2. **Three is enough to detect the problem, and not enough to characterise it.** Three draws
   locate a design's spread crudely; they cannot distinguish a design that fails one seed in
   three from one that fails one in ten. The per-design `range_*` columns in `SUBMISSION.csv`
   are a floor on the true spread, not an estimate of it.

3. **The anchor-RMSD gate is where the whole replicate result lives, and it is a knife-edge.**
   24 of the 27 failures are this criterion, all within 0.61 Å of the line, and relaxing it to
   2.5 Å would move 132/150 to 146/150. The gate is defensible — the anchor is the designed
   agonist motif and holding it matters — but the 132/150 headline is as much a statement about
   where the threshold was placed as about the designs.

4. **Receptor RMSD < 5.0 Å contributed nothing.** It never fired across 450 predictions. It
   should either be tightened to a value that discriminates or dropped from the gate's
   description, because as stated it implies a filter that is not filtering.

5. **Do not drop `glp1r-08` on the strength of its replicate record.** It holds the gate in 1 of
   3 replicates — but its failures are anchor RMSD 2.398 Å and 2.188 Å, both inside the marginal
   band that 24 of 24 failures occupy; it keeps all three of its hotspots in every replicate;
   and no replicate-1 metric distinguished it in advance. Removing it would mean screening on a
   criterion this document shows is at chance for exactly this purpose. It is submitted with
   `n_seeds_gate = 1` recorded in the file, so a reader who disagrees can drop it in one filter.

6. **The ranking should not be re-derived from ipSAE.** At r ≈ 0.40 between replicates, an
   ipSAE-ordered list is substantially a different list on the next seed. The `mean_ipsae` and
   `range_ipsae` columns are reported together for this reason: the spread (median 0.0124 across
   the ten submitted) is a meaningful fraction of the total spread of the cohort's means
   (0.8832 - 0.9180, i.e. 0.0348).
