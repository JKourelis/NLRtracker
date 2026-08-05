# NLRtracker

Find and annotate NLR immune receptors in a plant proteome.

Give NLRtracker a protein FASTA file. It tells you which proteins are NLRs, which subclass each one
belongs to, where every domain starts and stops, and how confident that subclass call is. It also
extracts each domain as its own FASTA file and draws the domain architecture.

Version 2.0.1 runs on InterProScan 6.0.1 and InterPro 109.0, so it now works on macOS and on Windows
under WSL2 as well as on Linux.

---

## What you get out

Nine kinds of output, from one command.

| | |
|---|---|
| A classification per protein | NLR or not, the subclass, and the evidence tier behind the call |
| Every domain with coordinates | non-overlapping, so you can use them without trimming |
| One FASTA per domain type | `NBARC.fasta`, `TIR.fasta`, `CC.fasta`, `LRR.fasta`, `CJID.fasta`, `ID.fasta` and more |
| Domain architecture figures | drawn automatically, paginated, and split by subclass |
| A self-contained HTML report | one file, opens in any browser, no internet needed |
| GFF3 annotations | for a genome browser |
| iTOL strips | for annotating a tree |
| A summary table | counts per class, per region and per signature |
| A provenance file | every tool version, every setting, so a result can be reproduced |

Optionally, with `-p`, it also places each protein in the RefPlantNLR phylogeny and tells you where
the tree and the domain architecture disagree.

### How well does it work

We tested it against the 481 NLRs in RefPlantNLR, whose subclasses were published independently.

NLRtracker gets 453 of them exactly right. For 23 more it reports an honest ambiguity that contains
the right answer, for example "CC-NLR or CCG10-NLR" when the evidence genuinely cannot separate the
two. One call in 481 is wrong.

That last number is the one worth watching. Reporting an ambiguity is not the same as being wrong,
and the tool is built to tell you which situation you are in.

---

## Can you run it?

You need seven things. Six are one-line installs, and the seventh is 40 GB of disk.

| | | |
|---|---|---|
| [Nextflow](https://www.nextflow.io/) | 25.10.4 or later | earlier versions cannot start InterProScan 6.0.1 |
| A container runtime | Docker, SingularityCE or Apptainer | |
| R | 4.1 or later | package list below |
| [HMMER](http://hmmer.org/) | 3.3 or later | `hmmsearch`, plus `hmmalign` if you use `-p` |
| [MEME Suite](https://meme-suite.org/) | 5.2 or later | for `fimo` |
| [VeryFastTree](https://github.com/citiususc/veryfasttree) | 4.0 or later | only if you use `-p`. [FastTree](http://www.microbesonline.org/fasttree/) 2.1 also works |
| Disk | about 40 GB | the InterProScan reference data |

You do not need to install InterProScan itself. Nextflow fetches it and its reference data the first
time you run.

Every R package comes from CRAN except `Biostrings`:

```r
install.packages(c("dplyr", "tidyr", "readr", "stringr", "tibble", "purrr",
                   "ggplot2", "scales", "jsonlite", "ape"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Biostrings")
```

`ape` is needed only for `-p`. Everything else is needed for every run.

`mafft` is not a requirement. It built the reference models once, and those models ship with the
repository, so a run never calls it. `module/reference/*.seed.aln` records what was aligned.

### Install

1. Install Nextflow, a container runtime, R, HMMER and the MEME Suite.
2. Install the R packages listed above.
3. Clone this repository.
4. Choose a location with 40 GB free.
5. Download the InterProScan reference data into it:

   ```bash
   nextflow run ebi-pf-team/interproscan6 -r 6.0.1 -profile docker,test \
     --datadir /path/to/interproscan6/data --interpro 109.0
   ```

6. Tell NLRtracker where that data lives:

   ```bash
   export INTERPROSCAN_DATADIR=/path/to/interproscan6/data
   ```

   The default location is `~/Bioinformatics/interproscan6/data`.

---

## Run it

```bash
./NLRtracker.sh -s proteins.fasta -o out_dir
```

That produces the full output tree described below, figures included.

Add `-p` to place every protein in the RefPlantNLR phylogeny as well:

```bash
./NLRtracker.sh -s proteins.fasta -o out_dir -p
```

For nucleotide input, set the sequence type:

```bash
./NLRtracker.sh -s transcripts.fasta -o out_dir -t n
```

C-JID detection is skipped for nucleotide input, because that model is a protein model.

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

### Re-run the classification without re-running InterProScan

The InterProScan step dominates the runtime. To change a signature table and see the effect, feed the
previous tool output back in:

```bash
./NLRtracker.sh -s proteins.fasta -o out2 \
  -i out/raw/interproscan/proteins.normalised.fasta.tsv \
  -f out/raw/fimo/fimo.gff \
  -X out/raw/hmmer/CJID.domtblout
```

---

## Output files

```
out_dir/
├── NLRtracker.tsv          classification: Status, Subclass (putative), support, architecture
├── Domains.tsv             every annotation, per protein. The table the figures are drawn from
├── summary.tsv             counts per class, per region, and per signature
├── report.html             self-contained interactive report
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
├── phylogeny/              only with -p
│   ├── NB-ARC.input.fasta  NB-ARC.aln  NB-ARC.tree
│   ├── clade_assignment.tsv       one row per protein: clade, parent chain, nearest
│   │                              RefPlantNLR relatives and their distances
│   ├── clade_summary.tsv          counts and recovery per reference clade
│   ├── discrepancies.tsv          architecture against phylogeny, categorised
│   ├── unplaced.tsv               proteins outside every reference clade
│   └── excluded_low_occupancy.tsv sequences too short to place reliably
└── raw/
    ├── interproscan/  fimo/  hmmer/
    └── logs/          including dropped_regions.tsv
```

### The columns in `NLRtracker.tsv`

| Column | Meaning |
|---|---|
| `seqname` | protein identifier |
| `Status` | `NLR`, `degenerate NLR`, `CCX`, `RPW8`, `TX`, `CJID`, `MLKL`, `NACHT`, `likely non-plant NLR` |
| `Subclass (putative)` | `CC-NLR`, `TIR-NLR`, `CCR-NLR`, `CCG10-NLR`, or an explicit ambiguity such as `CC-NLR or CCG10-NLR` |
| `subclass_support` | how the call was reached, and how much to trust it |
| `Domain` | full architecture, for example `(CC)(NBARC)(LRR)` |
| `Domain architecture simplified` | one letter per region, for example `CNL` |

---

## A subclass call is only as good as its evidence

Two proteins can both be called `CC-NLR` for very different reasons. One has a recognisable coiled-coil
domain at its N-terminus. The other has no N-terminal domain at all, and the call rests on a marker
inside the NB-ARC that happens to correlate with the CC lineage. Those are not equally trustworthy,
and a single subclass column cannot tell them apart.

The `subclass_support` column records which kind of evidence was used. On RefPlantNLR it predicts
reliability directly.

| Tier | What was found | Wrong calls in 481 |
|---|---|---:|
| `domain` | a diagnostic N-terminal domain | 0 |
| `probable` | a gated rule fired | 0 |
| `lineage` | an NB-ARC winged-helix marker only, so the call is reported as an explicit ambiguity | 0 |
| `supporting` | a low-coverage signature corroborated an existing call | 0 |
| `motif` | MEME motifs only | 0 |
| *(empty)* | no subclass evidence at all | 1 |

Read `(probable)` as trustworthy when it appears and meaningless when it does not. The gated rule
that produces it is designed to stay silent rather than guess.

The single wrong call in 481 sits in the bottom row, where there was no subclass evidence at all. That
is the pattern worth knowing: when the tool has nothing to go on, it says so.

---

## Why regions do not overlap

Raw signature coordinates overlap constantly. An LRR signature routinely runs into the NB-ARC, and a
TIR signature runs into the protein's N-terminus. Version 1.x reported those coordinates unchanged,
which meant extracted domain sequences needed trimming by hand, and figures drew domains on top of one
another.

Version 2.0 settles the boundaries instead. The NB-ARC is the anchor, a bare P-loop anchors the same
way, and the anchor is never trimmed.

1. Build the NB-ARC region. Merge the P-loop, HD1, winged-helix and NB-ARC motif evidence.
2. Truncate every N-terminal region, meaning TIR, CC, RPW8 and R1, to end before the NB-ARC starts.
3. Truncate the LRR to start after the NB-ARC ends.
4. Truncate the C-JID to start after the LRR ends.
5. Subtract every canonical region from the integrated-domain regions.

A region trimmed out of existence is dropped and then reported in `raw/logs/dropped_regions.tsv`,
grouped by the rule that removed it. A TIR lying entirely inside an NB-ARC call is a signal worth
seeing, not something to hide.

Every run asserts that the result is disjoint. That is what makes the per-region FASTA files usable
directly for alignment and phylogeny, with no further trimming.

### The integrated-domain region

An integrated domain is a domain fused into an NLR that does not belong to the NLR core, such as a
WRKY transcription factor domain or a heavy-metal-associated domain. Those are how many NLRs
recognise their targets, so they matter.

Every integrated domain contributes an `ID` region. That gives one consistent handle for "this NLR has
an integrated domain" across `Domains.tsv`, the figures, `ID.fasta` and the GFF3. The specific name is
kept in a parallel `DOMAIN` row, so a figure still labels it `WRKY`, `HMA` or `PKin`.

An integrated domain can sit on either side of the NB-ARC, so truncating it from one end would be
wrong. Subtracting the canonical regions can split one integrated domain into pieces, which is the
correct answer when a domain genuinely straddles a boundary.

---

## The biology lives in three tables, not in the code

Version 1.x encoded the rule for every signature as a 90-line branching statement inside the R script.
Updating InterPro meant editing that statement, and there was no way to see what any signature did
without reading code.

All of it is now data.

| File | Contents |
|---|---|
| `module/nlr_signatures.tsv` | 231 signatures, each with its database, InterPro entry, domain class, subclass evidence, hit counts and provenance |
| `module/nlr_classes.tsv` | 41 classes, each with its region, merge group, precedence, architecture letter, colour and shape |
| `module/nlr_subclass_rules.tsv` | the rules that need more than one signature |

Three consequences follow. An InterPro update becomes a data refresh rather than a code change.
Anyone can see what a signature does by reading one line. And the same tables can drive a different
implementation, which is what makes the tool portable beyond R.

`docs/signatures.md` explains, per domain class, which signatures belong to it, how they overlap, and
what a hit means for the classification.

### Update the tables after a new InterPro release

1. Download the new release with InterProScan.
2. Recover the 1.x script from git. It holds the original signature-to-class mapping:

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

4. Read the generator report.
5. Check the `ABSENT` list. Those signatures no longer exist and can never match.
6. Download the new entry list to `module/InterPro_<version>.list`.
7. Update `INTERPRO_VERSION` in `NLRtracker.sh`.
8. Run `tests/run_tests.sh --full`.
9. Compare the classification with the previous release.

---

## The speedup applies to a new genome, not a published one

InterProScan 6 splits the input into batches and runs every analysis once per batch. Its default batch
size of 5000 suits a whole proteome. On a smaller input it collapses to a single batch, which leaves
most of the machine idle.

NLRtracker sizes the batch from the input instead. Measured on 24 physical cores, with the Matches API
switched off so that every sequence was computed locally:

| Input | Upstream default | Auto-sized | Speedup |
|---|---:|---:|---:|
| 481 sequences | 777 s | 198 s | 3.92× |
| 6 000 proteins | 1369 s | 438 s | 3.13× |

Now the part that decides whether those numbers apply to you.

By default InterProScan looks each sequence up in the InterPro Matches API and computes only what is
missing. For a published proteome most sequences come back from that lookup, so the wall clock is set
by network latency rather than by your processor, and the batch size barely matters.

For a new genome nothing is precalculated, everything is computed locally, and the batch size decides
the runtime. That is the case the table above measures, and it is the case that actually hurts.

Two practical consequences:

1. Benchmark with `-L`, which disables the lookup. Otherwise you are timing EBI's latency.
2. Record which path a published result took. The same input can be served either way depending on
   release timing.

> Two runs of the same input at different batch sizes differ by roughly 0.1% of InterProScan rows, or
> 44 rows out of 34 732 at proteome scale. SMART's prefilter is not normalised for database size
> upstream, and CATH-Gene3D shifts with its sub-batch split. Keep the batch size constant across runs
> you intend to compare. `run_info.yaml` records the value used.

> Automatic batch sizing reads the physical core count, from `lscpu` on Linux and
> `sysctl hw.physicalcpu` on macOS. If neither answers it falls back to `getconf _NPROCESSORS_ONLN`,
> which counts logical processors and so over-counts on a machine with simultaneous multithreading.
> The batch size affects speed only, never results. Set `-b` to remove the guess.

---

## Platform support

InterProScan 6 removed the Java and Linux-only requirement, so NLRtracker runs on all three
platforms. The containers are `linux/amd64` only, which matters on Apple Silicon.

| Platform | Status |
|---|---|
| Linux | Native. This is the tested platform. |
| Intel Mac | Native `linux/amd64` containers. Install HMMER, the MEME Suite and VeryFastTree from Homebrew or Bioconda. |
| Apple Silicon | Runs under Docker Desktop's `amd64` emulation, and is much slower. Raise the Docker Desktop memory limit, because the default 2 GB is below what the InterProScan tasks request. |
| Windows | Run inside WSL2, which behaves as Linux. `NLRtracker.sh` is a bash script, so it does not run in PowerShell or `cmd`. |

Every component, by platform:

| Component | Linux | macOS | Windows |
|---|---|---|---|
| `NLRtracker.sh` | yes | yes | WSL2 only. It is bash, and it reads `lscpu` and `/proc` to size batches |
| InterProScan 6 (Nextflow and a container) | yes | yes | yes, inside WSL2 |
| `module/*.R` | yes | yes | yes |
| HMMER, MEME Suite | yes | yes | no practical native build. Use WSL2 |
| VeryFastTree | yes | yes | WSL2 |

The R code and the reference data are portable. The shell wrapper is the part that ties a run to a
POSIX system, and WSL2 satisfies it.

---

## What changed in 2.0.0

| | |
|---|---|
| InterProScan 6 | Runs through Nextflow and containers, with no local Java install. This is what removed the Linux-only constraint. |
| InterPro 109.0 | The bundled entry list had been frozen at InterPro 5.53-87.0 from 2021, holding 40 038 entries against 54 191 today. |
| 9 new NLR signatures | Three NB-ARC winged-helix Pfam entries, four NLR-specific LRR Pfam entries, an MLKL N-terminal CDD entry, and a four-helix-bundle superfamily. All curated against RefPlantNLR and a full proteome. See [`docs/signatures.md`](docs/signatures.md). |
| CCG10-NLR can now be assigned | It could previously appear only inside the ambiguous string `CC-NLR or CCR-NLR or CCG10-NLR`. |
| Regions are arbitrated | TIR, CC, RPW8, R1, LRR, C-JID and a general integrated-domain region no longer overlap. |
| One FASTA per region | `TIR.fasta`, `CC.fasta`, `LRR.fasta`, `CJID.fasta`, `ID.fasta` and more, alongside `NBARC.fasta`. |
| Figures are automatic | No separate `refplantnlR` step is needed. |
| An HTML report | One self-contained file, with no external assets. |
| Optional phylogeny | `-p` places each protein in the RefPlantNLR phylogeny and flags disagreement with the architecture. |
| Signature logic became data | A 90-line branching statement became three TSV tables that regenerate on each InterPro release. |
| Organised output | Results stay at the top level. Sequences, annotations, figures and raw tool output moved into subdirectories. |
| Provenance | `run_info.yaml` records every version, the applications used, and the batch size. |

---

## Upgrading from 1.x

Five things changed that will break an existing script.

**The output layout moved.** A script reading `out/NLR.fasta` now needs `out/sequences/NLR.fasta`, and
`out/NLR.gff3` becomes `out/annotations/NLR.gff3`.

**`NBARC.fasta` boundaries moved.** The new winged-helix Pfam entries extend the NB-ARC over the
winged-helix subdomain, so extracted NB-ARC sequences are generally longer and more complete.
Alignments and trees built from 1.x output are not directly comparable.

**Region coordinates moved** for TIR, CC, RPW8, LRR and C-JID, because they are now arbitrated rather
than taken raw from the signatures.

**`-i` now expects an InterProScan 6 TSV file, not an InterProScan 5 GFF3.** The two formats are
incompatible: InterProScan 5 put the accession in `Name=` with a separate `signature_desc=`, while
InterProScan 6 puts the description in `Name=` and the accession in `Alias=`. The parser rejects a
version-5 file rather than producing wrong results.

**`refplantnlR` needs 0.2.0 or later.** Earlier versions stop on the new region names with
*"Insufficient values in manual scale"*.

---

## Citation

If you use NLRtracker, please cite:

- **RefPlantNLR**. Kourelis J, Sakai T, Adachi H, Kamoun S (2021). *RefPlantNLR is a comprehensive
  collection of experimentally validated plant disease resistance proteins from the NLR family.* PLoS
  Biology 19(10): e3001124.
- **InterPro**. Blum M, et al. (2025). *InterPro: the protein sequence classification resource in
  2025.* Nucleic Acids Research 53(D1): D444–D456.
- **MEME Suite**. Bailey TL, et al. (2015). *The MEME Suite.* Nucleic Acids Research 43(W1): W39–W49.
- **NLR Annotator** motifs. Steuernagel B, et al. (2020). *The NLR-Annotator tool enables annotation
  of the intracellular immune receptor repertoire.* Plant Physiology 183(2): 468–482.
- **C-JID model**. Ma S, et al. (2020). *Direct pathogen-induced assembly of an NLR immune receptor
  complex to form a holoenzyme.* Science 370(6521): eabe3069.

---

## Version history

### 2.0.1

Two fixes to the optional phylogeny step. Nothing else changed, and classification is
untouched.

**Tip names in `phylogeny/NB-ARC.tree` are no longer corrupted.** `ape::write.tree` sanitises
labels silently, so `Bph19-1_(1)` came out as `Bph19-1_-1`. Names in the published tree could
therefore not be joined back to `NLRtracker.tsv`, and it affected exactly the proteins carrying
duplicated NB-ARC domains: Bph19, Pib, RGA2 and Prv. The tree is now published by substituting
quoted names into VeryFastTree's own output, so the file is what the tree program wrote with the
identifiers put back.

**Tree building now uses one thread by default.** VeryFastTree 4.0.5 is not reproducible above
one thread, and on this data it also returns a much worse tree. Each thread count run twice on
the 507-tip reference scaffold:

| Threads | Splits reported bad | Reproducible | Time |
|---:|---:|---|---:|
| 1 | 6 of 414 | yes | 2.71 s |
| 2 | 401 of 414 | no | 1.96 s |
| 4 | 378 of 414 | no | 1.42 s |
| 8 | 378 of 414 | no | 1.47 s |

The whole saving was 1.2 seconds. Pass `--threads` to `module/NLRtracker.R` to raise it if a
very large query set needs it, and record the value you used.

If you built a tree with 2.0.0, rebuild it.

### 2.0.0

- InterProScan 6.0.1 through Nextflow, and InterPro 109.0
- Runs on macOS and on Windows under WSL2
- 9 new NLR signatures, curated against RefPlantNLR and a full proteome
- CCG10-NLR assignable for the first time, through a gated rule reported as `(probable)`
- Arbitrated, non-overlapping regions: TIR, CC, RPW8, R1, NB-ARC, PLOOP, LRR, C-JID, ID
- One FASTA file per region
- Automatic domain-architecture figures, paginated and split by subclass
- A self-contained HTML report
- Optional phylogenetic classification against RefPlantNLR, with `-p`
- A `subclass_support` column recording the evidence tier behind each call
- Signature logic moved out of the code into three regenerable TSV tables
- Auto-sized InterProScan batching, 3 to 4 times faster than the upstream default on local compute
- `summary.tsv` and `run_info.yaml`
- Reorganised output tree. This is a breaking change.
- Hardened driver: `set -euo pipefail`, validated arguments, quoted paths

### 1.0.3 and earlier

See [`docs/CHANGELOG-1.x.md`](docs/CHANGELOG-1.x.md).

---

## Licence

MIT. See [`LICENSE`](LICENSE).
