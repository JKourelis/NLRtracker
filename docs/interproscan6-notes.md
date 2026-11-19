# InterProScan 6 notes

Working notes on InterProScan `6.0.3.0` (`InterPro 111.0`) behaviour observed while testing NLRtracker.

Official documentation: <https://interproscan6.readthedocs.io/latest/>.

## 1. TSV column layout

InterProScan `6` emits a 15-column TSV without a header. Two columns differ from early `6.x` documentation:

| Column | Content | Note |
|---:|---|---|
| 1 | Protein accession | |
| 2 | Sequence MD5 | |
| 3 | Sequence length | |
| 4 | Analysis | `Pfam`, `CDD`, `SUPERFAMILY`, `PRINTS`, `SMART`, `CATH-Gene3D`, `ProSiteProfiles` |
| 5 | Signature accession | |
| 6 | Signature description | Populated even for unintegrated signatures |
| 7 | Start | |
| 8 | Stop | |
| 9 | Score or E-value | |
| 10 | Match status / source | `module/R/01_parse.R` does not require `"T"` here |
| 11 | Date | |
| 12 | InterPro accession | `-` when unintegrated |
| 13 | InterPro description | `-` when unintegrated |
| 14 | GO annotations | `-` when absent |
| 15 | Pathway annotations | `-` when absent |

`module/R/01_parse.R` reads TSVs with `quote = ""` so bare double quotes inside InterPro entry names (such as `IPR036390`, `"Winged helix" DNA-binding domain`) do not shift columns or corrupt quotes.

## 2. InterProScan 5 GFF3 incompatibility

NLRtracker `1.x` parsed InterProScan `5` GFF3 files, whose `Name=` and `Alias=` attributes map differently from InterProScan `6`. Passing an InterProScan `5` output to `-i` triggers an explicit validation error in `read_interproscan()` rather than producing silent misannotations.

## 3. Applications and `prosite-patterns`

NLRtracker runs seven applications by default:

```text
pfam,cathgene3d,superfamily,prints,smart,cdd,prosite-profiles
```

In InterProScan `6.0.1`, `prosite-patterns` failed with `ps_scan.pl: command not found`. Although `6.0.3.0` executes `prosite-patterns` cleanly, NLRtracker keeps it disabled because none of the `231` curated NLR signatures relies on a PROSITE pattern (`prosite-profiles` covers all required PROSITE entries).

## 4. Matches API, offline `-L` mode, and `run_info.yaml`

By default, InterProScan checks the remote EBI Matches API and runs containers only for uncached sequences. Pass `-L` to force fully local execution (`noMatchesApi = true` and `skipInterproVersionCheck = true`), which guarantees reproducibility across InterPro updates without depending on remote lookup state.

`run_info.yaml` records `nlrtracker_version`, `interproscan_version`, `interpro_version`, `applications`, `interproscan_batch_size`, `seqtype`, `input`, input/classified/region/signature counts, `R_version`, and `run_date`.

## 5. Nextflow isolation and symlinked caches

`NLRtracker.sh` launches `nextflow` inside `<outdir>/raw/interproscan/` with absolute paths and physically resolved (`pwd -P`) `datadir`, `outdir`, `fasta`, and `NXF_HOME`. That prevents concurrent runs from colliding on `.nextflow/` session locks in the working directory and ensures Docker bind-mounts succeed when `~/.nextflow` is a symlink. Setting `subBatchSize` equal to `batchSize` in `nlrtracker.config` also avoids Nextflow `splitFasta` serialisation races on large inputs.

## 6. `module/json_to_tsv.py`

`module/json_to_tsv.py` converts InterProScan `6` JSON into the 15-column TSV format using string number parsing (`parse_float=str`, `parse_int=str`) and location-level E-value fallback for `CDD` matches.
