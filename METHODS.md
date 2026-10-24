# Methods

Everything below was run on a CentOS 7 / glibc 2.17 cluster (IBiSCo, V100-SXM2-32GB) with no
container runtime, which constrains every tool version; `protocol/envs/` carries the exact pins
and `protocol/envs/tool-commits.txt` the tool commits. Every script named here is in
`protocol/scripts/`; the originals live in `challenge2026/glp1r2026/scripts/` and the job logs in
`challenge2026/glp1r2026/logs/`.

The pipeline produces 140-residue de novo miniproteins that present the GLP-1 N-terminal
pharmacophore (residues 7-16) on a designed scaffold, and that are predicted to dock on the GLP-1R
extracellular domain with the anchor in its native pocket and the designed body engaging a chosen
hotspot patch. It has three stages, run as one SLURM job per shard so the GPU allocation is never
released between them:

```
  AF3 GLP-1R:GLP-1 complex (5 models)
        |  RFdiffusion3, motif scaffolding on peptide residues 1-10   GPU
        v
  backbones  --> screen: max CA deviation <= 2.0 A, <= 1 chainbreak
        |  soluble ProteinMPNN, 16 sequences each, anchor fixed       CPU
        v
  sequences  --> 3 per backbone advanced
        |
  Boltz-2 co-fold against the 409-aa receptor, cached receptor MSA    GPU
        |
        v
  measure: anchor RMSD vs AF3, hotspot contact, interface contacts, ipSAE
        |
        v
  gate -> rank -> deduplicate -> 50 per round
```

No restraints are used at any point. There is no restrained-then-released cycle here: the fold is
unconstrained from the first pass, and what makes the anchor land in its pocket is that RFD3 built
the backbone around it, not a contact term in Boltz.

## 0. Target preparation

| item | value |
|---|---|
| target | human GLP-1R, the 409-aa chain of `protocol/inputs/fold_glp1r_glp_job_request.json` |
| template | AlphaFold3 server job `glp1r_glp`, dialect `alphafoldserver` v3, modelSeed 1777850401 |
| template files | `fold_glp1r_glp_model_0.cif` .. `_model_4.cif` (receptor chain A 409 aa, peptide chain B 30 aa) |
| reference peptide | GLP-1 `HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR` |
| motif transplanted | peptide residues 1-10, `HAEGTFTSDV` (= GLP-1 His7-Val16) |
| receptor MSA | `glp1r2026/inputs/msa/glp1r_receptor_A.a3m`, 9929 sequences, cached from the AF3 server |
| designable face | the ECD; every hotspot used lies in receptor residues 24-95 |

The binder chain is always folded with `msa: empty` because it is de novo. **The receptor MSA is not
optional.** Measured directly (job 1088573, `protocol/scripts/h4_msatest.sbatch`), the same binder
scored ipTM 0.841 with the cached a3m and ipTM 0.316 without it; Boltz 2.2.1 reads the AF3 a3m
file directly. The MSA is 5.3 MB and is not copied into this submission
(`protocol/inputs/MSA_NOTE.txt` records its path).

## 1. Motif definition and hotspot patches

The anchor is fixed, not searched: RFD3 is given the AF3 peptide's first ten residues as a motif
segment and asked to grow a 130-residue scaffold off its C-terminus. What varies between rounds is
the **hotspot patch** on the receptor that the diffusion model is oriented toward. Hotspots are
specified at atom level in the RFD3 input JSONs (`protocol/inputs/rfd3/`), and
`infer_ori_strategy: hotspots` places the binder from them:

| round input | hotspot atoms | residue identities |
|---|---|---|
| `prod_r3.json` | A24:OE1,NE2  A25:NH1,NH2  A27:CD1  A42:OG1 | Gln24, Arg25, Leu27, Thr42 |
| `w2_r86_89.json` | A86:CD2  A87:CD1  A88:CD1  A89:NE2 | Leu86, Trp87, Leu88, Gln89 |
| `w3_r27_42_86.json` | A27:CD1  A42:OG1  A86:CD2 | Leu27, Thr42, Leu86 |
| `w4_r27_42_86_95.json` | A27:CD1  A42:OG1  A86:CD2  A95:CD2 | Leu27, Thr42, Leu86, Leu95 |

Ser26 (`A26:OG`) appears in the earlier `glp1r_ecd.json` and was **dropped**: 0 of 30 designs
reached it, median distance 6.8 A (recorded in the `h5_epitope.py` header). A95:CD2 was chosen by
burial with `protocol/scripts/pick_atom.py` — 30 heavy-atom neighbours within 8 A, the most exposed
atom in the set, against A86:CD2 41, A42:OG1 58 and A27:CD1 78. Patch geometry was measured before
launching with `protocol/scripts/hot_geom.py`: A27-A42 6.1 A, A42-A86 11.9 A, A27-A86 16.5 A, so
A27+A42 are effectively one site and A86 a second site 12-16 A away; the anchor sits 37-47 A from
all of them.

## 2. Environments

Three, and they cannot be merged. Pins from `protocol/envs/`.

| env | python | used for | key pins |
|---|---|---|---|
| `rfd3` | 3.12.14 | RFdiffusion3, soluble ProteinMPNN, **and every measurement script** | torch 2.6.0, rc-foundry 0.2.0, atomworks 2.2.1, biotite 1.4.0, numpy 2.5.3, scipy 1.16.3, hydra-core 1.3.7 |
| `boltz2` | 3.12.14 | Boltz-2 co-folding, every round | torch 2.6.0, boltz 2.2.1, numpy 1.26.4, numba 0.61.0, scipy 1.13.1, pytorch-lightning 2.5.0 |
| `ftmap` | 3.11.16 | physicochemical analysis | numpy 2.3.5, scipy 1.17.1, pandas, biotite |

ipSAE is a checkout, not a package: `tools/ipsae/ipsae.py` at commit
`6174cf9e71cb1bd660cc805856a18c4871a6dec3`.

Every job begins `source challenge2026/scripts/env.sh`, which is not optional. It redirects
`HF_HOME`, `TORCH_HOME`, `PIP_CACHE_DIR` and `XDG_CACHE_HOME` to scratch, points `C26_WEIGHTS` at
scratch, sets `CONDA_SOLVER=classic` (the shipped libmamba solver is broken here), disables core
dumps, and raises `ulimit -u` from the cluster default of 100 to 8192.

## 3. Backbone generation — RFdiffusion3

`rc-foundry` 0.2.0 (`rfd3`), one V100 per shard, inside a retry loop of 3 attempts:

```
rfd3 design out_dir=$DIFF inputs=$G/inputs/rfd3/<round>.json ckpt_path=$CKPT \
     n_batches=1 diffusion_batch_size=4 seed=$((BASE + SH)) \
     inference_sampler.step_scale=3 inference_sampler.gamma_0=0.2 \
     prevalidate_inputs=True skip_existing=True
```

Contig, identical in every scaffolding round:

```
B1-10,130-130,/0,A1-409
```

Chain B residues 1-10 are the AF3 peptide motif; `130-130` is a fixed 130-residue designed
extension; `/0` is the chain break; `A1-409` is the receptor held fixed. Binder length is therefore
exactly **140 aa in every design** — 10 anchor + 130 designed — which is why no length sweep appears
in the gates. `is_non_loopy: true` and `dialect: 2` on every input. Each array task uses its own
diffusion seed (`seed` is a real config field), so shards cannot duplicate designs; the bases were
7000 (p2, p3), 9000 (p10), 12000 (p12) and 15000 (p13), plus the shard index.

Five AF3 models are supplied as separate input entries, so one shard at the default
`diffusion_batch_size=4` yields 20 backbones and, at `--per-backbone 3`, 60 folds. p10 (the smoke
test) was launched smaller, at `BATCH=3 PERBB=2`: 15 backbones and 30 folds per shard.

RFD3 renumbers: the binder comes back as chain A starting at 1 and the receptor as chain B. Each
design's JSON carries `diffused_index_map` (`B<i> -> A<j>`), and **the motif position is read from
that map per design, never assumed** (`h2_mpnn.py`).

## 4. Backbone screen

Before any sequence design, `h2_mpnn.py` reads each design's own RFD3 metrics JSON and keeps only
backbones with

```
max_ca_deviation <= 2.0 A   AND   n_chainbreaks <= 1
```

`max_ca_deviation` is how far the scaffolded motif moved from the transplanted coordinates, so this
discards backbones that did not actually hold the pharmacophore. Skipped backbones are listed by
name and deviation in the stage log.

## 5. Sequence design — soluble ProteinMPNN

Driven as a Python API (`MPNNInferenceEngine`, which ships inside the `rfd3` env), not the CLI, so
the checkpoint loads once per shard. Soluble weights `soluble_v_48_020.pt` with
`is_legacy_weights=True`, `device="cpu"`.

```
--nseq 16   --temperature 0.1   --seed $((BASE + SH))
--acid-pos 3   --max-ca-dev 2.0   --max-breaks 1
omit: ["CYS"]   (global, every position)
```

Two requirements are **enforced by omission and then asserted**, because a logit bias cannot
guarantee an identity:

1. **The 10-residue anchor must survive byte-identical.** `--fixed_residues` is mutually exclusive
   with the chain selectors in this engine, so instead every *designable* position is enumerated in
   `designed_residues` and the motif is compared character-for-character against the input structure
   afterwards. A mismatch is exit code 7, not a warning.
2. **Binder position 3 must be Asp or Glu.** With C-terminal-only extension the motif occupies
   binder positions 1-10, so position 3 falls *inside* the fixed anchor — the native Glu satisfies
   the constraint and the script verifies that rather than silently forcing nothing. (The same
   argument, applied to a design where position 3 is free, omits the other 18 amino acids there.)

16 sequences are sampled per backbone and the first **3** are advanced to folding
(`h3_make_yamls.py --per-backbone 3`). MPNN scores are recorded but never filtered on.

## 6. Co-folding — Boltz-2

`h3_make_yamls.py` writes one YAML per sequence, binder as chain A with `msa: empty` and receptor as
chain B with the cached a3m. The manifest is written **one level above** the YAML directory, because
`boltz predict` parses every file in the directory it is given and dies on a stray `.json`.

```
boltz predict $Y --out_dir $OUT --cache $SHMC --model boltz2 \
      --write_full_pae --override --accelerator gpu --devices 1 \
      --num_workers 2 --preprocessing-threads 2 --diffusion_samples 1 --seed 1
```

(`--num_workers 1 --preprocessing-threads 1` in the 2-CPU jobs: p15, p16, p17, p18.)

Operational details that are load-bearing:

- the Boltz checkpoint is copied to `/dev/shm/boltzcache_$SLURM_JOB_ID` per job and `mols` /
  `boltz2_aff.ckpt` symlinked beside it, so the folding stage does not read weights off Lustre;
- `NUMBA_CACHE_DIR` is exported per shard, because concurrent shards otherwise race on Boltz's
  `@njit(cache=True)` cache;
- folding runs under `timeout ${LEFT}s` with `DEADLINE_S` set to ~87% of the wall limit, so the job
  stops accepting new folds and still measures whatever finished rather than being killed mid-write.
  Return code 124 is an expected outcome, not a failure.

Chain assignment is **never** taken from the label. `h3_measure.py`, `h5_epitope.py`,
`h7_interactions.py` and `build_q.py` all identify the receptor as the longer chain and then confirm
the binder by length and sequence identity against the manifest, because Boltz renumbers chains
from 1 and the binder's label depends on YAML order.

## 7. Measurement

Four independent measurements per fold, all geometric except ipSAE.

**7a. Anchor placement — `h3_measure.py`.** Superpose the refold onto the AF3 model using the
**receptor CA atoms only** (Kabsch), then, with that single transform, measure the RMSD of the ten
motif residues against where AF3 put the peptide's first ten. Superposing on the receptor and
scoring the motif keeps the two independent: an identical pose gives ~0 A, a binder that dragged the
pharmacophore elsewhere gives a large number. Reported as `motif_rmsd_vs_af3` and `receptor_rmsd`.

**7b. Epitope engagement — `h5_epitope.py`.** Minimum heavy-atom distance from the **designed region
only** (binder 11-140, the anchor explicitly excluded) to each requested hotspot residue; contact
cutoff **5.0 A**; `n_hotspots_hit` is how many are within it. `--name-prefix` restricts scoring to
one lineage, needed from p14 onward where one prediction directory mixes backbones from campaigns
with different hotspot sets.

> **A high confidence score does not mean the binder hit the intended epitope.** On this target
> ipSAE barely separates on-target from off-target. Measured over the 1200 p2 folds: median ipSAE
> 0.882 at 4/4 hotspots against 0.871 at 0/4 — a separation of 0.011. Over the 480 p3 folds, 0.894
> against 0.870. In the final p15 gated set the ordering even inverts (4 hits median 0.888, 3 hits
> median 0.894). Epitope adherence is therefore established by the geometric check in 7b and by
> nothing else; ipSAE is used only to order designs that check has already verified.

**7c. Interface contacts — `h7_interactions.py`.** Distance-only proxies, stated as such:

| contact | definition |
|---|---|
| salt bridge | charged N of K/R/H within **4.0 A** of carboxylate O of D/E |
| hydrogen bond | N/O donor-acceptor pair within **3.5 A**, **no angle term** — over-counts, comparative only |
| hydrophobic | side-chain carbon of A/V/L/I/M/F/W/P within **4.5 A** of the same on the other side |

Only binder-to-receptor pairs are counted. Reported separately for the whole binder, the designed
region 11-140 (which contacts the ECD) and the anchor 1-10 (which contacts the TMD).

**7d. ipSAE.** Computed explicitly from the Boltz PAE `.npz` plus the model CIF at cutoffs **15/15**,
after folding, in its own pass:

```
python tools/ipsae/ipsae.py $pae $cif 15 15
```

Both the optimistic `max` row (`ipsae`) and the conservative direction (`ipsae_cons`) are recorded.
The gap between them is large for some designs — in the p15 gated set `ipsae` ranges 0.857-0.918
while `ipsae_cons` ranges 0.314-0.896 — so one value alone is not enough. ipSAE writes its output
next to the **input** files, not the working directory.

## 8. Gate, score and packaging

**Gate first, then rank** (`rank_campaign.py`), applied to every fold, advancing every passer rather
than a fixed top-N:

```
motif RMSD vs AF3 < 2.0 A   AND   >= 3 hotspots contacted   AND   receptor RMSD < 5.0 A
```

**Score** is the mean of three min-max normalised terms computed over the gated set only: motif RMSD
(inverted), ipSAE at 15/15, and binder-chain mean pLDDT.

Binder pLDDT is read per-residue from the Boltz `plddt_*.npz` and sliced to the binder chain by
walking the CIF's CA order and checking the slice length against the manifest's `binder_len`.
**`complex_plddt` is recorded but never ranked on**: the complex is 409 receptor residues to 140
binder residues, so it is three-quarters a constant.

**Packaging** (`best_dir.py`): order by hotspot count first and composite score second (`--order
tier`, which matters because the score itself ignores hotspot count), then greedy-cluster at **70%
identity over the designed region** keeping the best-ordered member, and take the top **50**. Each
round directory gets `metrics.csv`, `candidates.fasta`, `designed_region.fasta`,
`structures/` (mmCIF complexes) and `binder_only/` (binder PDB with pLDDT in the B-factor column),
plus a `README.txt` generated from that round's own numbers — every quantity in it is an argument or
derived from the data, so a round README always describes its own directory.

## 9. The rounds as run

Each row is one SLURM array; 20 shards is both the `gpus` QOS limit (`MaxJobsPU=20`) and the number
of reachable free GPUs, so throttles are planned at that.

| round | run tag | what changed | shards | refolds | gated | full epitope in the packaged 50 |
|---|---|---|---|---|---|---|
| p2 | `p2_1088883` | hotspots A24/25/27/42 | 20 | 1200 | 869 (72%) | 50/50 at 4/4 |
| p3 | `p3_1088996` | hotspots A86-89 (the Trp87 patch) | 8 | 480 | 309 (64%) | 50/50 at 4/4 |
| p10 | `p10_1089071` | combined patch A27/42/86, smoke test | 3 | 90 | 64 (71%) | 50/50 at 3/3 |
| p12 | `p12_1089085` | p10 at production scale | 20 | 1197 | 779 (65%) | 50/50 at 3/3 |
| p13 | `p13_1089119` | p12 plus Leu95 | 20 | 1200 | 822 (68%) | 50/50 at 4/4 |
| p14 | `p14_1092297` | charge-complementary redesign, anchor designable | 20 | 600 | 432 (72%) | 50/50 at own set* |
| p15 | `p15_1099477` | anchor repaired under the agonist consensus | 10 | 210 | 187 (89%) | 50/50 at own set* |
| p18 s2 | `p18s2_1106925` | seed-2 replicate of the top 150 p15 | 10 | 150 | 138 (92%) | — |
| p18 s3 | `p18s3_1106926` | seed-3 replicate of the top 150 p15 | 10 | 150 | 135 (90%) | — |

\* p14 and p15 mix lineages with different hotspot-set sizes in one directory. `best_dir.py` was
originally called with `--hotspot-count 4` for both, which counted only designs reaching 4 of 4 and so
reported p15 as 34/50 — a p10 or p12 design can contact at most 3 and was never counted. Both
`p*_best/README.txt` files have since been corrected to count each design against its own lineage's
set, giving 50/50 for both rounds; the figures above are the corrected ones. Per-lineage
full-engagement rates are in the round logs (p15: p10 100%, p12 95%, p13 89%, p3 88%).

Campaign-wide hotspot-count histograms are in each `p*_best/README.txt`. p2 and p3 use the same
diffusion seed base and the same MPNN seed base (37 + shard), so the two patches are directly
comparable; p10 used MPNN seed 71 + shard, p12 131 + shard, p13 191 + shard.

## 10. p14 — charge-complementary redesign (`h6_mpnn_charge.py`)

Rounds p2-p13 fixed the anchor and left the designed surface to ProteinMPNN's default composition,
which came out **acidic**: measured over the 200 packaged input complexes, median q(11-140) was -6
(p3), -5 (p10), -4 (p12), -5 (p13). Against an acidic epitope that is the wrong sign.

**Grounding** (`contact_charge.py` over the same 200 complexes, printed in `ctcharge_1090959.out`):

```
binder 11-140 contacts the ECD : 44 res, net charge -3.0  -> wants BASIC
binder 1-10   contacts the TMD : 27 res, net charge +3.0  -> wants ACIDIC
```

Two changes from `h2_mpnn.py`:

1. **The anchor becomes designable.** All 140 positions are designed, with position 1 forced to HIS
   and position 3 to ASP/GLU by omission. The native `HAEGTFTSDV` already satisfies both, so this
   relaxes the other eight positions rather than fighting the native sequence. Positions 2 and 4-10
   carry a mild acidic bias (0.8).
2. **Charge is modulated per position from measured geometry, not globally.** Each binder position
   within 5 A of the receptor gets a bias chosen by its partner composition: partner acidic -> K/R
   up, partner basic -> D/E down-weighted, partner apolar -> A/V/L/I/M/F, partner polar ->
   N/Q/S/T/Y/H. Interior positions are left unbiased so MPNN can still build a core. Only
   `bias_per_residue` is used, never the global `bias`, because per-residue entries override the
   global one and mixing them makes the net effect at a position ambiguous. Exposed non-interface
   positions (burial < 85 neighbours within 8 A of CB) carry the bulk charge bias; the interface set
   alone (~36 of 130) is too small to move the net charge.

Parameters as run:

```
--nseq 32 --keep 3 --temperature 0.15 --seed $((500 + SH)) \
--bias-weight 1.6 --target-charge 5.5 --charge-window 1.5 --charge-weight 1.2 \
--charge-basis full --bulk-charge-bias 0.85 --charge-rounds 3
```

Sequences are **selected, not merely sampled**: many are drawn per backbone and ranked on
salt-bridge potential against the actual partners, with a steep penalty outside the charge window,
so the charge target is met by selection rather than by an aggressive bias.

**Result.** 600 refolds from 200 input backbones, whose inputs were 100% fully engaged: motif
RMSD < 2 A in 486 (81%), full epitope retained in 457 (76%), >= 3 hotspots in 526 (88%), gate passed
by 432 (72%), 10 complete misses. Retention differed sharply by lineage — p10 87%, p12 83%, p13 80%,
p3 55% full-epitope retention — so the p3 (Trp87) lineage is the one that loses most when its
surface is redesigned.

**The cost.** Designing the anchor freely destroyed the agonist pharmacophore. Measured with
`consensus.py`: **0 of 600 p14 anchors satisfied the agonist consensus, and 100% broke the invariant
Gly4 and Thr5.** The 39 distinct anchors p14 produced (`HAEEEFFAEV`, `HPEEEFFDEV`, ...) are listed
in `p14_best/README.txt`.

## 11. p15 — anchor repair under the agonist consensus (`h8_anchor.py`)

p15 repairs what p14 broke without giving back the surface p14 won. **Positions 11-140 are held at
their p14 identities and the change is asserted, not assumed** — a candidate whose designed region
differs from the parent at any position is rejected (`designed_region_changed`), so the validated
charge-complementary surface is preserved exactly and only the 10-residue pharmacophore moves.

**The consensus** (`consensus.py`), framed so that consensus position 1 = GLP-1 His7 = binder
position 1, with non-canonical residues (X = Aib / Acb / a-Me-Leu, which MPNN cannot produce)
replaced by the canonical alternatives from the counts:

| pos | allowed | top | conservation across the 23 agonists |
|---|---|---|---|
| 1 | H, Y, F | H | variable, 82.6% H |
| 2 | G, A, S | G | variable; X=Aib leads at 34.8% but is non-canonical |
| 3 | E, D | E | acid class 65.2%; Q excluded to keep the acid requirement |
| 4 | G | G | **invariant, 100%** |
| 5 | T | T | **invariant, 100%** |
| 6 | F | F | **invariant, 100%** |
| 7 | T, I | T | near-invariant, 95.7% T |
| 8 | S | S | **invariant, 100%** |
| 9 | D | D | **invariant, 100%** |
| 10 | Y, V, L, K | Y | variable; hydrophobic class 56.5% |

Each position is restrained **by omission**, which forces an identity; the five invariant positions
are fully determined and the variable ones carry a bias of 1.0 toward the consensus top residue.
Native GLP-1 `HAEGTFTSDV` satisfies every restraint, so this is a repair, not an exploration.

**Inputs** were selected by `pick_inputs.py` from the p14 low-charge band: gated, full epitope for
their own lineage, and q(11-140) in +4..+9. That band was characterised first with
`band_quality.py` — 142 folds from 101 backbones, 107/142 (75%) gated, 94/142 (66%) gated and fully
engaged from 71 backbones. Deduplicating by backbone gave **71 distinct inputs** (p13 22, p12 20,
p10 15, p3 14), score 0.323-0.916, q(11-140) 5-8. Parameters: `--nseq 16 --keep 3 --temperature 0.2
--seed $((900 + SH)) --top-bias 1.0`.

**Result**, verified on the refolded structures rather than the design-stage log (`p15_check.py`):

- **210/210 anchors satisfy all ten restraints**, including every invariant position, against 0/600
  for p14;
- 22 distinct anchors; identity to native GLP-1 median 9/10, min 7, max 10; the commonest is native
  `HAEGTFTSDV` (37 designs);
- 190/210 (90%) anchor RMSD < 2 A, median 1.35 A; 195/210 (93%) full epitope; **187/210 (89%) pass
  the gate** — the highest gate rate of any round;
- median 7 unique salt bridges, 30 H-bond pairs, 80 hydrophobic contacts; q(11-140) median +8,
  inherited from p14 and unchanged by design;
- gated-set confidence: ipSAE median 0.892, binder pLDDT median 0.903, ipTM median 0.924, pDockQ
  median 0.634.

## 12. Validation

Three independent checks, all unconstrained, all on the full p15 cohort.

### 12a. Three-seed replicates (`p18_reseed.sbatch`, `mk_reseed.py`, `consistency.py`)

The top 150 p15 designs by composite score were re-folded at Boltz seeds 2 and 3. Inputs are
**byte-identical** to the originals — same binder sequence, same receptor, same cached MSA, same
manifest names — so the only difference between runs is `--seed`, which makes this a replicate test
and not a re-design. Both replicate arms ran as 10 shards each.

Spread across the three seeds, per design, median over 150 designs:

| metric | median range | median SD | worst range |
|---|---|---|---|
| ipSAE | 0.016 | 0.007 | 0.058 |
| binder pLDDT | 0.005 | 0.002 | 0.016 |
| ipTM | 0.009 | 0.004 | 0.041 |
| anchor RMSD | 0.48 A | 0.21 A | 1.52 A |
| receptor RMSD | 0.38 A | 0.16 A | 1.78 A |
| unique salt bridges | 2.0 | 0.8 | 5.0 |
| hydrophobic contacts | 11 | 5 | 37 |

Does each criterion hold in **every** seed, not on average?

| criterion | all 3 seeds | any seed | fragile |
|---|---|---|---|
| receptor RMSD < 5 A | 150 (100%) | 150 | 0 |
| binder pLDDT >= 0.80 | 150 (100%) | 150 | 0 |
| ipSAE >= 0.50 | 150 (100%) | 150 | 0 |
| >= 3 hotspots | 148 (99%) | 150 | 2 |
| full epitope contact | 137 (91%) | 150 | 13 |
| anchor RMSD < 2 A | 133 (89%) | 150 | 17 |
| salt bridges >= 5 | 135 (90%) | 150 | 15 |

The combined gate holds in all three seeds for **132 of 150 (88%)**; 18 designs pass sometimes and
fail sometimes. Per-seed the gate rate is 150/150 at seed 1 (by construction — these were selected
there), 138 (92%) at seed 2 and 135 (90%) at seed 3.

**The confidence metrics reproduce; the geometry does not, as well.** Spearman between seed pairs:

| metric | p15~s2 | p15~s3 | s2~s3 |
|---|---|---|---|
| binder pLDDT | 0.928 | 0.941 | 0.950 |
| ipTM | 0.857 | 0.832 | 0.803 |
| hydrophobic contacts | 0.692 | 0.738 | 0.698 |
| salt bridges | 0.573 | 0.628 | 0.593 |
| anchor RMSD | 0.467 | 0.506 | 0.566 |
| **ipSAE** | **0.462** | **0.318** | **0.418** |

ipSAE ranks almost as poorly as the geometric quantities here. Its *values* are tight (median range
0.016 across seeds) but its *ordering* is not, because the whole cohort sits in a 0.06-wide band.
Ranking a cohort this uniform by ipSAE is close to ranking it at random; the gate, not the score, is
what separates designs.

### 12b. Monomer control (`p16_mono.sbatch`, `mk_monomer_yamls.py`, `measure_monomer.py`)

All 210 p15 binder sequences were folded **alone**, no receptor chain, `msa: empty`, seed 1, and
compared residue-for-residue with the same chain inside its complex. Superposition is on binder CA
atoms only.

- monomer pLDDT (whole chain): median 0.880 (p10 0.827, p90 0.919); **200/210 >= 0.80, 0/210 < 0.70**
- designed region 11-140: median 0.901; anchor 1-10: median 0.584
- RMSD monomer vs bound, by window (`refine_mono.py`, `p16_monomer_windows.csv`):

| window | p10 | median | p90 | within 3 A |
|---|---|---|---|---|
| all 140 | 1.41 | 10.28 | 17.58 | 65 (31%) |
| designed 11-140 | 0.76 | 5.81 | 11.74 | 86 (41%) |
| **domain 30-140** | **0.24** | **0.42** | **1.77** | **203 (97%)** |
| domain 40-140 | 0.19 | 0.29 | 0.72 | 207 (99%) |
| core 50-120 | 0.16 | 0.24 | 0.45 | 207 (99%) |

**The whole-chain RMSD is misleading and the windowed one is the real answer.** The folded domain
(residues 30-140) is essentially identical bound and unbound — median 0.42 A, 203/210 within 3 A,
and the distribution is unimodal with 174 designs under 1 A. What moves is the N-terminal 10-30
segment that carries the anchor into its pocket, which is exactly the part expected to be
receptor-ordered: its monomer pLDDT is 0.584 against 0.901 for the body. Quoting the all-140 median
of 10.28 A as "the binder does not hold its fold" would be wrong. By lineage, domain-within-3 A is
100% for p10, p12 and p13 and 83% for p3.

### 12c. Self-association control (`p17_dimer.sbatch`, `mk_dimer_yamls.py`, `dimer_analysis.py`)

The same 210 sequences folded as **two copies of themselves**, no receptor. The question is not
whether Boltz draws two chains together — with two chains and nothing else in the box it generally
will, and the median A-B contact count here is 732 — but how confident that interface is relative to
the same binder's receptor interface.

| | median |
|---|---|
| binder-receptor ipTM (from the p15 complexes) | 0.923 |
| binder-binder ipTM | 0.262 |
| difference | **+0.661 in favour of the receptor** |

`ipTM(self) >= 0.80`: **1 of 210**. `>= 0.85`: 0. `< 0.60`: 200 of 210 (95%). By lineage the median
self-ipTM is 0.243 (p3), 0.334 (p10), 0.259 (p12), 0.247 (p13). Self-association is not a liability
in this cohort.

## 13. Packaging the deliverable

`build_q.py` takes a list of design names (`qmap.txt`), pulls each design's structure from all three
replicate runs (`p15_1099477`, `p18s2_1106925`, `p18s3_1106926`), converts mmCIF to PDB and joins
the per-replicate metrics from `seed_consistency.csv`. `pack10.py` writes the final directory.

PDB rather than mmCIF because 3Dmol's mmCIF reader silently truncated these files in the browser.
Chain A = binder, residues 1-140 (1-10 is the anchor); chain B = GLP-1R, residues 1-409; both start
at 1; the B-factor column carries per-residue pLDDT on a 0-100 scale. **The chain counts are
verified, not assumed**: a structure whose chains are not exactly 140 and 409 is reported as a
problem rather than written.

Cohort figures for the ten shipped designs, median across the ten: mean ipSAE 0.902 (median spread
across replicates 0.012), mean binder pLDDT 0.906, mean anchor RMSD 0.77 A (median spread 0.46 A).
**9 of 10 hold the full gate in all three replicates and 9 of 10 keep full hotspot contact in all
three.**

> **The hotspot set differs by lineage and the counts are not comparable across lineages.** p3 uses
> A86-A89 (4), p10 and p12 use A27/A42/A86 (3), p13 uses A27/A42/A86/A95 (4). A design showing 3 of
> 3 is not weaker than one showing 4 of 4. `n_hotspots_total` in `SUBMISSION.csv` and
> `metrics/seed_consistency.csv` records which set applies.

## 14. Reproducing this

```bash
source challenge2026/scripts/env.sh     # ALWAYS first: caches to scratch, ulimit -u 100 -> 8192
conda activate rfd3                     # backbones, MPNN, all measurement
conda activate boltz2                   # co-folding

# one scaffolding round (20 shards, diffusion + MPNN + fold + measure in one job)
sbatch --array=0-19%20 protocol/scripts/p13_prod.sbatch
sbatch protocol/scripts/p13_finish.sbatch        # ipSAE -> rank_campaign -> best_dir

# charge redesign on the packaged backbones, then anchor repair
sbatch --array=0-19%20 protocol/scripts/p14_redesign.sbatch
sbatch protocol/scripts/p14_finish.sbatch
sbatch --array=0-9%10  protocol/scripts/p15_anchor.sbatch
sbatch protocol/scripts/p15_rank.sbatch

# validation
SEED=2 TOPN=150 sbatch --array=0-9%10 protocol/scripts/p18_reseed.sbatch
SEED=3 TOPN=150 sbatch --array=0-9%10 protocol/scripts/p18_reseed.sbatch
sbatch protocol/scripts/p16_mono.sbatch          # monomer control
sbatch protocol/scripts/p17_dimer.sbatch         # self-association control
python protocol/scripts/consistency.py analysis \
       p15_1099477=<run>,p18s2_1106925=<run>,p18s3_1106926=<run> analysis/seed_consistency.csv
```

## 15. Things that will bite

- **`scripts/env.sh` first, always.** `ulimit -u` is 100 on login *and* compute nodes; conda's
  download thread pool dies well before that.
- **Pre-create output directories, and retry.** Lustre `mkdir` from compute nodes intermittently
  fails with EIO. Every job wraps `mkdir -p` in a 5-attempt loop with a 20 s backoff, and the
  diffusion stage retries 3 times.
- **End every sbatch with an explicit status line.** A trailing `echo` makes SLURM record COMPLETED
  on a failed job; these scripts print `=== exit: $rc ===` and exit on it.
- **Export a per-shard `NUMBA_CACHE_DIR`.** Concurrent Boltz shards otherwise race on the numba
  cache.
- **ipSAE writes next to its INPUT**, not the working directory, and must not be run with
  `>/dev/null` alone without checking that the `*_15_15.txt` file appeared.
- **Never identify a chain by its label.** Boltz renumbers from 1 and the binder's label depends on
  YAML order; RFD3 returns the binder as chain A and the receptor as chain B, the opposite of the
  Boltz co-fold. Every script here resolves chains by length and confirms by sequence.
- **RFD3 names designs from the input JSON key, which doubles the run prefix**
  (`prod_r3_prod_r3_m1_...`). `best_dir.py` strips one copy; anything joining on names must handle it.
- **`--fixed_residues` is mutually exclusive with the chain selectors** in this MPNN engine. Enumerate
  the designable positions in `designed_residues` instead, and verify the fixed region afterwards.
- **A bias cannot force an identity.** Every hard constraint in this pipeline — position 1 = His,
  position 3 = acid, the ten consensus positions — is implemented by omitting the other residues and
  then asserted on the output.
- **The epitope check must exclude the anchor.** The anchor is in its pocket by construction; scoring
  it would report engagement for a design whose entire designed body missed.
- **Quote the windowed monomer RMSD, not the whole-chain one.** See 12b: the same cohort reads as
  10.28 A or 0.42 A depending on whether the receptor-ordered N-terminal segment is included.
