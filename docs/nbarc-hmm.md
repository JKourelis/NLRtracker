# The NB-ARC region model

`module/reference/NB-ARC_region.hmm` defines NB-ARC region coordinates. InterPro signatures supply domain and lineage classification evidence, while this profile sets the boundaries of every `NB-ARC` region in `Domains.tsv` and `sequences/NBARC.fasta`.

The model (`LENG 382`) is locked. `module/reference/NB-ARC_region.LOCK` records its SHA-256 checksum, and `tests/run_tests.sh` fails if the HMM file changes.

## 1. Why a full-length model is used instead of `PF00931`

Pfam `PF00931` (`NB-ARC`) is shorter than a complete plant NB-ARC and sits close to its gathering threshold (`GA 23.5`) on divergent NLRs:

- **Partial coverage:** `PF00931` had model length `175` aa in Pfam `38.2` (`InterPro 109.0`) and shortened to `152` aa in Pfam `39.0` (`InterPro 111.0`), whereas full NB-ARC regions (NBD, HD1, and WHD segments) have a RefPlantNLR median length of **356 aa** (`n = 508`, p25 = `336`, p75 = `361`, max = `424`).
- **Divergent receptor sensitivity:** Under Pfam `38.2`, the divergent TIR-NLR outgroup (`TNP`) scored `26.0` to `26.1` against `GA 23.5` (compared to `72` to `89` against `NB-ARC_region.hmm`). Under Pfam `39.0` (`PF00931.30`), `TN2` drops below gathering threshold so `Status` calls `degenerate NLR`, while `NB-ARC_region.hmm` still recovers its `NB-ARC` region (`168–374`).

## 2. Model construction and thresholds

1. Seed: `467` RefPlantNLR NB-ARC regions of at least `300` aa plus `4` `TNP` outgroup regions (`471` total).
2. Alignment: `mafft --auto` once at build time (`module/reference/NB-ARC_region.seed.aln`). `mafft` is not a runtime dependency.
3. Profile: `hmmbuild --amino -n NB-ARC_region` -> `NB-ARC_region.hmm` (`LENG 382`).

Two calibrated score cutoffs separate NLR detection from domain boundary placement:

| Parameter | Score | Purpose | Empirical basis |
|---|---:|---|---|
| `NBARC_HMM_PROTEIN_GA` | 40 | Accepts an HMM-only protein as NLR-like | Minimum RefPlantNLR score is `40.7`; CDC48 and AAA+ ATPases peak at `15.8` |
| `NBARC_HMM_DOMAIN_GA` | 20 | Keeps partial NB-ARC copies inside an NLR | Recovers `RRS1` partial second NB-ARC (`1068–1187`, score `36`) |

HMM-only candidates (hits not covered by a member-database NB-ARC signature) also require corroborating NLR evidence from `NBARC_NLR_SUPPORT` and zero hits in `NBARC_EXCLUDE_SUPPORT` (`AAA+`/`CDC48` signatures in `module/R/03_regions.R`). Completeness is reported as `(hmm_to - hmm_from + 1) / 382`.

## 3. Handling duplicated, partial, and split NB-ARCs

### Adjacent envelopes are trimmed rather than merged

When two HMM envelopes overlap on a protein (`Bph1/9-1` partial copy `138–388` and complete copy `381–768`), the earlier envelope ends at `next_start - 1`. Verified multi-NB-ARC controls asserted by `tests/run_tests.sh --full`:

| Protein | NB-ARC regions | Biological nature |
|---|---:|---|
| `Bph1/9-1`, `Bph1/9-2`, `Bph1/9-7`, `Bph1/9-9`, `Pib`, `RGA2` | 2 | Genuine partial duplication |
| `Prv^1`, `Prv^2`, `RPP2A` | 2 | Genuine tandem NB-ARC |
| `RRS1-R` | 2 | Second partial NB-ARC recovered at score `36` |
| `Pi2` | 1 | Stray LRR-embedded motif excluded from forming a second NB-ARC |

### Evidence-based joining of interrupted envelopes (`2.0.3`)

`hmmsearch` occasionally splits a single NB-ARC across an insertion (for example `Rpp4C4`, which hits model positions `27–148` at `142–286` and `147–373` at `312–554`). `join_split_nbarc()` in `module/R/01_parse.R` runs inside `build_regions()` right before `apply_hmm_nbarc()` and joins consecutive envelopes `A` then `B` when all five conditions hold:

1. The model continues: `B$hmm_from >= A$hmm_to - 1`.
2. The model does not restart: `B$hmm_from > A$hmm_from` (separating insertions from tandem copies).
3. Protein coordinates do not overlap: `B$start > A$end`.
4. No unrelated domain is contained inside the gap `[A$end + 1, B$start - 1]` (ignoring `NBARC`, `NACHT`, `PLOOP`, `nbarc-motif`, `Winged-helix`, and `nacht-winged-helix`).
5. At least one corroboration rule fires: a canonical FIMO NB-ARC motif run spans the gap (strongest), a signature-derived `NBARC` or `NACHT` domain spans the gap, or the total uncorroborated backstop gap across the chain is at most `30` residues.

Joined chains must also stay within `382` model positions (`max_hmm_to - min_hmm_from + 1 <= 382`). On RefPlantNLR, `Rpp4C4` joins on motif continuity into one `NB-ARC` region, while all tandem and partial duplications stay separate.

### Motif hits extend regions without creating them

Isolated FIMO motif hits cannot create an `NB-ARC` region on their own, preventing stray 15-aa artefacts inside LRR repeats (`Pi2`, `Pi9`, `PigmR`, `PigmS`, `Piz-t`).

## 4. Reference phylogeny scaffold

`module/reference/NB-ARC_region.reference.fasta` and `module/reference/NB-ARC_region.reference.itol` hold `500` NB-ARC regions extracted under InterPro `111.0` from the `481` published RefPlantNLR proteins (`496` regions, including `Pb1`, `RXL`, and all tandem/partial pairs such as `RRS1-R_(1)` and `RRS1-R_(2)`) plus the `4` `TNP` outgroup proteins. `module/nlr_clades.tsv` defines `23` active clades (`1` outgroup plus `22` anchored clades).

## 5. Rebuilding the model

Keep the model locked unless a major release explicitly intends to shift coordinates:

1. Regenerate the seed alignment from curated full-length NB-ARC regions.
2. Run `mafft` and `hmmbuild`, recalibrating `NBARC_HMM_PROTEIN_GA` and `NBARC_HMM_DOMAIN_GA` on RefPlantNLR.
3. Update `module/reference/NB-ARC_region.LOCK` with the new SHA-256 digest.
4. Run `tests/run_tests.sh --full` and explain any shift in the `457 / 23 / 1` RefPlantNLR benchmark.
