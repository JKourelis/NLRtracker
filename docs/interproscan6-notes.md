# InterProScan 6 notes

Working notes on InterProScan 6.0.1, kept because two of them are places where the official
documentation is wrong and the parsers here depend on the real behaviour.

Official documentation: <https://interproscan6.readthedocs.io/latest/>. Note `latest`. Several
older manifest URLs still in circulation return 404.

---

## 1. The TSV columns

InterProScan 6 writes a 15-column TSV. Two columns do not match the published description.

| # | Content | Note |
|---:|---|---|
| 1 | Protein accession | |
| 2 | Sequence MD5 | |
| 3 | Sequence length | |
| 4 | Analysis | `Pfam`, `CDD`, `SUPERFAMILY`, … |
| 5 | Signature accession | |
| 6 | Signature description | present even when the signature is unintegrated |
| 7 | Start | |
| 8 | Stop | |
| 9 | Score / E-value | |
| 10 | **Match source** | documented as the boolean `T`. It is not. |
| 11 | Date | |
| 12 | InterPro accession | `-` when unintegrated |
| 13 | InterPro description | `-` when unintegrated |
| 14 | GO annotations | optional |
| 15 | **Pathway annotations** | the column order around 14/15 differs from the documentation |

Consequences for this repository:

- Do not test column 10 for `T`. `module/R/01_parse.R` ignores it.
- An unintegrated signature gives `-` in columns 12 and 13 but still carries its description in
  column 6, so a signature must never be discarded for lacking an InterPro accession. Several NLR
  signatures are unintegrated.

## 2. IPS5 GFF3 is not a drop-in substitute

NLRtracker 1.x parsed the InterProScan 5 GFF3. The `Name=` and `Alias=` attribute semantics are
effectively swapped between the two versions, so an IPS5 GFF3 read as if it were IPS6 output yields
signature accessions in the wrong field. This is why 2.x reads the TSV and why the 1.x fixtures were
removed from `sample_data/` rather than kept alongside.

## 3. `prosite-patterns` is broken in 6.0.1

The `prosite-patterns` application fails with `ps_scan.pl: command not found`. It is therefore not
in the `APPLICATIONS` list in `NLRtracker.sh`. `prosite-profiles` works and is used.

## 4. The Matches API changes what is computed, not what is reported

By default InterProScan looks sequences up in the InterPro Matches API and computes only the
remainder. `-L` disables the lookup and forces a fully local run.

Results are equivalent. Runtime is not, and neither is reproducibility: a lookup result depends on
the state of a remote service. Use `-L` for anything that must be reproducible, and see
[`performance.md`](performance.md) for what that costs.

`run_info.yaml` records whether the API was used.

## 5. `json_to_tsv.py`

`module/json_to_tsv.py` converts InterProScan 6 JSON to the 15-column TSV, so existing JSON output
can be reused without recomputing. It reproduces the real TSV byte for byte, which required two
non-obvious choices:

- Parse numbers as strings (`parse_float=str, parse_int=str`). Round-tripping through a float
  changes the formatting of E-values.
- For CDD, fall back to the location-level E-value when the match level has none.
