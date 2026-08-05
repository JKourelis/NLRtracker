# NLRtracker

Extract and functionally annotate NLR immune receptors from plant proteomes.

NLRtracker classifies NLRs from protein or transcript sequences using the core features of the
[RefPlantNLR](https://doi.org/10.1371/journal.pbio.3001124) dataset — a curated collection of
experimentally validated plant disease resistance proteins from the NLR family.

**Version 2.0.0** runs on **InterProScan 6.0.1 / InterPro 109.0**, which makes it portable to
macOS and Windows. See [Platform support](#platform-support) and [Upgrading from 1.x](#upgrading-from-1x).

---

## What is new in 2.0.0

| | |
|---|---|
| **InterProScan 6** | Nextflow + containers, no local Java install. This is what removes the Linux-only constraint. |
| **InterPro 109.0** | The bundled entry list was frozen at InterPro 5.53-87.0 from 2021 — 40 038 entries against 54 191 today. |
| **9 new NLR signatures** | Three NB-ARC winged-helix Pfams, four NLR-specific LRR Pfams, an MLKL N-terminal CDD, a 4H-bundle superfamily. All curated against RefPlantNLR and a full proteome — see [`docs/signatures.md`](docs/signatures.md). |
| **CCG10-NLR is assignable** | Previously it could only ever appear inside the ambiguous string `"CC-NLR or CCR-NLR or CCG10-NLR"`. |
| **Arbitrated regions** | TIR, CC, RPW8, R1, LRR, C-JID and a general ID region are now first-class and **do not overlap**. |
| **Per-region FASTA** | `TIR.fasta`, `CC.fasta`, `LRR.fasta`, `CJID.fasta`, `ID.fasta`, … alongside `NBARC.fasta`. |
| **Figures are automatic** | No separate `refplantnlR` step needed. |
| **Signature logic is data** | A ~90-line `case_when` chain became three TSV tables that regenerate on each InterPro release. |
| **Organised output** | Results at the top level; sequences, annotations, figures and raw tool output in subdirectories. |
| **Provenance** | `run_info.yaml` records every version, the applications, and the batch size. |

Measured against the 481 published RefPlantNLR subclass assignments: **94.2% exact**, a further
4.8% ambiguous-but-containing-the-truth, **one wrong call**.

---

## Requirements

| | | |
|---|---|---|
| [Nextflow](https://www.nextflow.io/) | **≥ 25.10.4** | InterProScan 6.0.1 refuses to start below this |
| Container runtime | Docker, SingularityCE or Apptainer | |
| R | ≥ 4.1 | see the package list below |
| [HMMER](http://hmmer.org/) | ≥ 3.3 | `hmmsearch` for the NB-ARC and C-JID models. `hmmalign` as well, with `-p` |
| [MEME Suite](https://meme-suite.org/) | ≥ 5.2 | `fimo` |
| [VeryFastTree](https://github.com/citiususc/veryfasttree) | ≥ 4.0 | only with `-p`. [FastTree](http://www.microbesonline.org/fasttree/) 2.1 also works |
| Disk | ~40 GB | the InterProScan reference data |

You do not need a local InterProScan installation. Nextflow fetches the workflow and the reference
data on first use.

R packages, all from CRAN except `Biostrings`:

```r
install.packages(c("dplyr", "tidyr", "readr", "stringr", "tibble", "purrr",
                   "ggplot2", "scales", "jsonlite", "ape"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Biostrings")
```

`ape` is needed only for `-p`. Everything else is needed for every run.

`mafft` is **not** a requirement. It built the shipped reference models once, and those models are
in the repository, so a run never calls it. `module/reference/*.seed.aln` records what was aligned.

### Install

1. Install Nextflow, a container runtime, R, HMMER and the MEME Suite.
2. Install the R packages listed above.

3. Clone this repository.
4. Download the InterProScan reference data. Choose a location with 40 GB free:

   ```bash
   nextflow run ebi-pf-team/interproscan6 -r 6.0.1 -profile docker,test \
     --datadir /path/to/interproscan6/data --interpro 109.0
   ```

5. Tell NLRtracker where that data lives. Either export the variable or pass `-D`:

   ```bash
   export INTERPROSCAN_DATADIR=/path/to/interproscan6/data
   ```

   The default is `~/Bioinformatics/interproscan6/data`.

---

## Usage

```bash
./NLRtracker.sh -s proteins.fasta -o out_dir
```

That produces the full output tree described below, figures included.

### Options

```
Required:
  -s FILE     Amino acid (or nucleotide, with -t n) sequence FASTA
  -o DIR      Output directory

Reuse existing tool output (skips that step):
  -i FILE     InterProScan 6 TSV result
  -f FILE     FIMO result (.gff)
  -X FILE     HMMER result (hmmsearch --domtblout)

Reference data:
  -m FILE     meme.xml for FIMO            [module/meme.xml]
  -x FILE     HMM for hmmsearch            [module/abe3069_Data_S1.hmm]
  -d FILE     InterPro entry list          [module/InterPro_109.0.list]
  -D DIR      InterProScan data directory  [$INTERPROSCAN_DATADIR]

Execution:
  -t p|n      Sequence type                [p]
  -c INT      CPUs per InterProScan task   [2]
  -b INT      InterProScan batch size      [auto]
  -w INT      Max concurrent tasks         [all available processors]
  -P NAME     Container profile: docker | singularity | apptainer  [docker]
  -p          Also classify phylogenetically against RefPlantNLR
  -L          Disable the InterPro Matches API (fully local, reproducible)
  -k          Keep the Nextflow work directory
  -h          Help
```

### Reuse a previous run

The InterProScan step dominates the runtime. To re-run only the classification after changing a
signature table, feed the previous tool output back in:

```bash
./NLRtracker.sh -s proteins.fasta -o out2 \
  -i out/raw/interproscan/proteins.normalised.fasta.tsv \
  -f out/raw/fimo/fimo.gff \
  -X out/raw/hmmer/CJID.domtblout
```

### Nucleotide input

```bash
./NLRtracker.sh -s transcripts.fasta -o out_dir -t n
```

C-JID detection is skipped for nucleotide input, because the HMM is a protein model.

---

## Output

```
out_dir/
├── NLRtracker.tsv          classification: Status, Subclass (putative), support, architecture
├── Domains.tsv             every annotation, per protein — the table figures are drawn from
├── summary.tsv             counts per class, per region, and per signature
├── run_info.yaml           provenance: tool versions, InterPro release, applications, batch size
├── figures/
│   ├── NLR.pdf             domain architecture, paginated at 50 proteins
│   ├── NLR-associated.pdf
│   └── by_subclass/*.pdf
├── sequences/
│   ├── NLR.fasta / NLR.lst
│   ├── NLR-associated.fasta / NLR-associated.lst
│   ├── NBARC.fasta / NBARC_deduplicated.fasta
│   └── TIR.fasta  CC.fasta  RPW8.fasta  R1.fasta  LRR.fasta  CJID.fasta  ID.fasta
├── annotations/
│   ├── NLR.gff3 / NLR-associated.gff3
│   └── iTOL.txt / iTOL_dedup.txt
└── raw/
    ├── interproscan/  fimo/  hmmer/
    └── logs/          incl. dropped_regions.tsv
```

### `NLRtracker.tsv`

| Column | Meaning |
|---|---|
| `seqname` | protein identifier |
| `Status` | `NLR`, `degenerate NLR`, `CCX`, `RPW8`, `TX`, `CJID`, `MLKL`, `NACHT`, `likely non-plant NLR` |
| `Subclass (putative)` | `CC-NLR`, `TIR-NLR`, `CCR-NLR`, `CCG10-NLR`, or an explicit ambiguity such as `CC-NLR or CCG10-NLR` |
| `subclass_support` | **how** the call was reached — see below |
| `Domain` | full architecture, e.g. `(CC)(NBARC)(LRR)` |
| `Domain architecture simplified` | one letter per region, e.g. `CNL` |

### `subclass_support` — read this before trusting a subclass

The support tier records the kind of evidence behind each call, and on RefPlantNLR it predicts
reliability:

| Tier | Meaning | Wrong calls in 481 |
|---|---|---:|
| `domain` | a diagnostic N-terminal domain was found | 0 |
| `probable` | a gated rule fired; correct when positive, silent when negative | 0 |
| `lineage` | an NB-ARC winged-helix marker only, so the call is reported as an explicit ambiguity | 0 |
| `supporting` | a low-coverage signature corroborated an existing call | 0 |
| `motif` | MEME motifs only | 0 |
| *(empty)* | no subclass evidence at all | 1 |

A subclass carrying `(probable)` is trustworthy when present and proves nothing when absent.

---

## Regions and how their boundaries are decided

Raw signature coordinates overlap: an LRR signature routinely extends into the NB-ARC, and a TIR
signature into its N-terminus. Version 1.x emitted those coordinates unchanged, so extracted
domain sequences needed manual trimming and figures drew domains on top of one another.

Version 2.0 arbitrates them. NB-ARC is the anchor — a bare P-loop anchors the same way — and it is
never trimmed:

1. Build the NB-ARC region. Merge the P-loop, HD1, winged-helix and NB-ARC motif evidence.
2. Truncate every N-terminal region — TIR, CC, RPW8, R1 — to end before the NB-ARC starts.
3. Truncate the LRR to start after the NB-ARC ends.
4. Truncate the C-JID to start after the LRR ends.
5. Subtract every canonical region from the integrated-domain (`ID`) regions.

Regions trimmed out of existence are dropped **and reported** in
`raw/logs/dropped_regions.tsv`, grouped by the rule that removed them. A TIR lying entirely inside
an NB-ARC call is a signal worth seeing, not something to hide.

The result is a disjoint region set, which is asserted on every run. Per-region FASTA files are
therefore directly usable for alignment and phylogeny without further trimming.

### The `ID` region

Every integrated domain contributes an `ID` region, giving one consistent handle for "this NLR has
an integrated domain" across `Domains.tsv`, the figures, `ID.fasta` and the GFF3. The specific
domain name is kept in a parallel `DOMAIN` row, so figures still label it `WRKY`, `HMA`, `PKin`
and so on.

---

## The signature tables

All biological logic lives in three generated TSV files, not in R code:

| File | Contents |
|---|---|
| `module/nlr_signatures.tsv` | 231 signatures → database, InterPro entry, domain class, subclass evidence, hit counts, provenance |
| `module/nlr_classes.tsv` | 41 classes → region, merge group, precedence, simplified-architecture letter, colour, shape |
| `module/nlr_subclass_rules.tsv` | rules needing more than one signature |

To update after a new InterPro release:

1. Download the new release with InterProScan.
2. Recover the 1.x script, which holds the original signature-to-class mapping. It is not
   duplicated in the working tree, because git already has it:

   ```bash
   git show v1.0.3:module/NLRtracker.R > /tmp/NLRtracker.legacy.R
   ```

3. Run the generator against the new `entries.json`:

   ```bash
   python3 module/build_signature_table.py \
     --entries  $INTERPROSCAN_DATADIR/interpro/110.0/entries.json \
     --legacy   /tmp/NLRtracker.legacy.R \
     --outdir   module \
     --interpro-version 110.0
   ```

4. Read the generator report. Check the `ABSENT` list — those signatures no longer exist and can
   never match.
5. Download the new entry list to `module/InterPro_<version>.list`.
6. Update `INTERPRO_VERSION` in `NLRtracker.sh`.
7. Re-run `tests/run_tests.sh --full` and compare the classification with the previous release.

`docs/signatures.md` explains, per domain class, which signatures belong to it, how they overlap,
and what a hit means for the classification.

---

## Performance

InterProScan 6 parallelises by splitting the input into batches and running every analysis once per
batch. The upstream default batch size of 5000 is tuned for whole proteomes, and on smaller inputs
it collapses to a single batch that leaves most of the machine idle.

NLRtracker sizes the batch from the input instead. Measured on 24 physical cores, **with the
Matches API disabled** so these are local-compute numbers:

| Input | Upstream default | Auto-sized | Speedup |
|---|---:|---:|---:|
| 481 sequences | 777 s | **198 s** | **3.92×** |
| 6 000 proteins | 1369 s | **438 s** | **3.13×** |

### What that speedup does and does not apply to

Both runs used `--no-matches-api`, so every sequence was computed locally. **That is the case the
speedup applies to, and it is the normal case when annotating a new genome**, where nothing is
precalculated.

With the Matches API **on** — the default — InterProScan retrieves precomputed annotations for any
input protein already present in UniParc. For a published proteome most sequences come back from
the lookup, wall clock becomes network-bound rather than compute-bound, and batching matters
correspondingly less. Do not expect 3–4× on a well-known proteome; do expect it on novel sequences.

Two practical consequences:

1. Benchmark with `-L` (`--no-matches-api`), or you are timing EBI's latency rather than your
   machine.
2. Record which path a published result took. The same input can be served either way depending on
   release timing, so a run is not reproducible without that.

Override with `-b` if you need a specific batch size, and `-c` / `-w` for threads per task and
total concurrency.

> **Batch size is recorded in `run_info.yaml`, and it matters.** Two runs of the same input at
> different batch sizes differ by roughly 0.1% of InterProScan rows — 44 of 34 732 at proteome
> scale — because SMART's prefilter is not normalised for database size upstream, and CATH-Gene3D
> shifts with its sub-batch split. Keep the batch size constant across runs you intend to compare.

---

## Platform support

InterProScan 6 removes the Java and Linux-only constraint, so NLRtracker now runs on all three
platforms. The containers are `linux/amd64` only, which matters on Apple Silicon.

| Platform | Status |
|---|---|
| **Linux** | Native. This is the tested platform. |
| **Intel Mac** | Native `linux/amd64` containers. Install HMMER, the MEME Suite and VeryFastTree from Homebrew or Bioconda. |
| **Apple Silicon** | Works under Docker Desktop's `amd64` emulation, and is materially slower. Raise the Docker Desktop memory limit — the default 2 GB is below what the InterProScan tasks request. |
| **Windows** | Run inside WSL2, which behaves as Linux. `NLRtracker.sh` is a bash script and does not run in PowerShell or `cmd`. |

Every component by platform:

| Component | Linux | macOS | Windows |
|---|---|---|---|
| `NLRtracker.sh` | yes | yes | WSL2 only — it is bash, and it reads `lscpu` and `/proc` to size batches |
| InterProScan 6 (Nextflow + container) | yes | yes | yes, inside WSL2 |
| `module/*.R` | yes | yes | yes |
| HMMER, MEME Suite | yes | yes | no practical native build; use WSL2 |
| VeryFastTree | yes | yes | WSL2 |

The R code and the reference data are portable. The shell wrapper is what ties a run to a
POSIX system, and WSL2 satisfies it.

> **Automatic batch sizing** reads the **physical** core count, from `lscpu` on Linux and
> `sysctl hw.physicalcpu` on macOS. If neither answers, it uses `getconf _NPROCESSORS_ONLN`, which
> counts logical CPUs and therefore over-counts on a machine with simultaneous multithreading. The
> batch size affects speed only, never results. Set `-b` to remove the guess.

---

## Upgrading from 1.x

**The output layout changed.** Results stay at the top level; sequences, annotations, figures and
raw tool output moved into subdirectories. Scripts that read `out/NLR.fasta` need
`out/sequences/NLR.fasta`, and `out/NLR.gff3` becomes `out/annotations/NLR.gff3`.

**`NBARC.fasta` boundaries changed.** The new NB-ARC winged-helix Pfams extend the NB-ARC over the
winged-helix subdomain, so extracted NB-ARC sequences are generally longer and more complete.
Alignments and trees built from 1.x output are not directly comparable.

**Region coordinates changed** for TIR, CC, RPW8, LRR and C-JID, because they are now arbitrated
rather than raw signature coordinates.

**`-i` now expects an InterProScan 6 TSV, not an InterProScan 5 GFF3.** The formats are
incompatible: InterProScan 5 put the accession in `Name=` with a `signature_desc=`; InterProScan 6
puts the description in `Name=` and the accession in `Alias=`. The parser rejects a version-5 file
rather than producing wrong results.

**`refplantnlR` needs 0.2.0 or later.** Earlier versions abort on the new region names with
*"Insufficient values in manual scale"*.

---

## Citation

If you use NLRtracker, cite:

- **RefPlantNLR** — Kourelis J, Sakai T, Adachi H, Kamoun S (2021). *RefPlantNLR is a
  comprehensive collection of experimentally validated plant disease resistance proteins from the
  NLR family.* PLoS Biology 19(10): e3001124.
- **InterPro** — Blum M, et al. (2025). *InterPro: the protein sequence classification resource in
  2025.* Nucleic Acids Research 53(D1): D444–D456.
- **MEME Suite** — Bailey TL, et al. (2015). *The MEME Suite.* Nucleic Acids Research 43(W1):
  W39–W49.
- **NLR Annotator** motifs — Steuernagel B, et al. (2020). *The NLR-Annotator tool enables
  annotation of the intracellular immune receptor repertoire.* Plant Physiology 183(2): 468–482.
- **C-JID HMM** — Ma S, et al. (2020). *Direct pathogen-induced assembly of an NLR immune receptor
  complex to form a holoenzyme.* Science 370(6521): eabe3069.

---

## Version history

### 2.0.0

- InterProScan 6.0.1 via Nextflow; InterPro 109.0
- Runs on macOS and Windows (WSL2)
- 9 new NLR signatures, curated against RefPlantNLR and a full proteome
- CCG10-NLR assignable for the first time, via a gated rule reported as `(probable)`
- Arbitrated, non-overlapping regions: TIR, CC, RPW8, R1, NB-ARC, PLOOP, LRR, C-JID, ID
- Per-region FASTA extraction
- Automatic domain-architecture figures, paginated and per subclass
- `subclass_support` column recording the evidence tier behind each call
- Signature logic externalised into three regenerable TSV tables
- Auto-sized InterProScan batching — 3–4× faster than the upstream default
- `summary.tsv` and `run_info.yaml`
- Reorganised output tree (**breaking**)
- Hardened driver: `set -euo pipefail`, validated arguments, quoted paths

### 1.0.3 and earlier

See [`docs/CHANGELOG-1.x.md`](docs/CHANGELOG-1.x.md).
