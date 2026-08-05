# The NB-ARC region model

`module/reference/NB-ARC_region.hmm` defines the NB-ARC region. InterPro signatures still contribute
evidence for classification, but the **boundaries** of the region come from this model.

The model is locked. `module/reference/NB-ARC_region.LOCK` records its checksum, and
`tests/run_tests.sh` fails if the file changes. Rebuilding it changes every NB-ARC coordinate the
tool has ever emitted, so it must be a deliberate, versioned act.

---

## Why not Pfam PF00931

PF00931 is the obvious candidate and is the wrong one, for two measured reasons.

**It covers about half the region.** PF00931 has `LENG 175`. The NB-ARC regions this tool defines have
a median length of **356 aa** (n = 508, p25 = 336, p75 = 361, max = 424). The NB-ARC is
NBD-HD1-WHD; aligning to PF00931 discards the HD1 and WHD parts that the definition exists to
capture.

**It sits on a detection cliff.** PF00931's curated gathering threshold is `GA 23.5`. Against it:

| Sequence | Best-domain score |
|---|---:|
| ZAR1 | 197 |
| Rx | 143 |
| MLA1 | 112 |
| Roq1 | 83 |
| RPS4 | 64 |
| divergent TIR-NLR outgroup | **26.0 – 26.1** |

The outgroup passes by 2.5 bits. Anything slightly more divergent falls below GA and disappears
entirely. A model built from full-length NB-ARC regions removes that cliff.

---

## How it was built

```
1. seed  = RefPlantNLR NB-ARC regions >= 300 aa (467 of 508) + the outgroup regions
2. MSA   = mafft, once
3. HMM   = hmmbuild  ->  module/reference/NB-ARC_region.hmm     (LENG 382)
4. lock  = sha256 into NB-ARC_region.LOCK
```

Regions under 100 aa were excluded from the seed; a fragment would distort the model. `mafft` ran
once, at build time, and is **not** a runtime dependency. `NB-ARC_region.seed.aln` is the alignment
that produced the model, so the build is auditable without re-running anything.

Everything after step 3 is profile alignment (`hmmalign`), which is deterministic and does not shift
when the query set changes. That is what makes coordinates and clade assignments stable across
proteomes.

---

## Calibrated thresholds

Two thresholds, both set from data rather than chosen:

| Threshold | Value | Basis |
|---|---:|---|
| `NBARC_HMM_PROTEIN_GA` | 40 | the lowest-scoring genuine RefPlantNLR NLR scores 40.7 |
| `NBARC_HMM_DOMAIN_GA` | 20 | RRS1's partial second NB-ARC scores 36; a domain floor below that keeps it |

The protein threshold decides whether a protein is NLR-like at all. The domain threshold decides
whether an individual match becomes a region. They are separate because a genuine partial domain
scores well below a genuine full-length one.

`completeness` is reported per region as `(hmm_to - hmm_from + 1) / 382`, so a partial is visible in
the output rather than silently indistinguishable from a complete one.

---

## What the model changed

Applied across RefPlantNLR and two proteomes:

- **478 / 249 / 1003** NB-ARC boundaries replaced by HMM envelope coordinates
- **2 / 9 / 13** proteins gained an NB-ARC region that signatures had missed
- **1 / 8 / 16** candidate matches rejected below threshold

Classification agreement with the 481 published RefPlantNLR subclasses did not move: 453 exact,
23 honestly ambiguous, 1 wrong.

### Partial and duplicated NB-ARCs

Partial NB-ARCs are biologically real, so there is no length filter. Two rules protect them:

**Adjacent envelopes are trimmed, never merged.** Where two HMM envelopes overlap, the earlier one
ends at `next_start - 1`. Without this, a tandem pair fuses into one region. Verified survivors:

| Protein | NB-ARC regions | Nature |
|---|---:|---|
| `Bph1/9-1`, `Bph1/9-2`, `Bph1/9-7`, `Bph1/9-9`, `Pib`, `RGA2` | 2 | genuine partial duplication |
| `Prv^1`, `Prv^2`, `RPP2A` | 2 | genuine tandem |
| `RRS1-R` | 2 | second copy partial, invisible to signatures — recovered by the model |

**A region needs domain evidence to exist.** Motif hits may extend a region but never create one.
This removes eight spurious 15 aa "NB-ARC regions", each traceable to a single stray FIMO hit — for
example the `RNBS-C` hit inside an LRR repeat region in `Pi2`, `Pi9`, `PigmR`, `PigmS` and `Piz-t`.
Length was the wrong filter for this: a stray motif is noise at any length, and a genuine partial is
real at any length.

`tests/run_tests.sh --full` asserts all four cases.

---

## Sequences that are NB-ARC-like but not NLRs

The model finds P-loop NTPases outside the NLR family, as any NB-ARC model will. Of 23 candidates
detected by the model but not by any NLR signature:

- **9** carry genuine independent NLR evidence (scores 56.6 – 132.3)
- **6** are CDC48 / AAA+ ATPases (scores 12.4 – 15.8) — well below the protein threshold
- **8** remain unclear

The score separation between the two confident groups is what `NBARC_HMM_PROTEIN_GA = 40` sits in.
Support lists in `module/R/03_regions.R` (`NBARC_NLR_SUPPORT`, `NBARC_EXCLUDE_SUPPORT`) name the
signature classes that confirm or veto an HMM-only call.

---

## Rebuilding it

Do not, unless a release intends to. If it must happen:

1. Regenerate the seed from the current NB-ARC regions.
2. `mafft` and `hmmbuild` as above.
3. Recalibrate both thresholds against RefPlantNLR — they are properties of the model, not constants.
4. Update `NB-ARC_region.LOCK` with the new checksum.
5. Re-run `tests/run_tests.sh --full`. Treat any movement in 453 / 23 / 1 as a finding to explain,
   not a number to update.
