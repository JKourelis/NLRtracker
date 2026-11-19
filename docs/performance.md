# Performance and batch sizing

InterProScan accounts for nearly the entire runtime of an NLRtracker job. This page records measured timings so the wrapper defaults in `NLRtracker.sh` can be verified directly.

Test hardware: 24 physical cores (`48` logical threads, `2` threads per core), NVMe storage, Docker.

## Batch size tuning

InterProScan 6 divides input proteins into batches of `batchSize` sequences and fans each batch across the seven requested analyses. Nextflow schedules tasks automatically; `NLRtracker.sh` tunes only batch granularity.

InterProScan's upstream default of `batchSize = 5000` collapses small datasets into one batch, leaving most CPU cores idle. Measured wall-clock times with the InterPro Matches API disabled (`-L`):

| Input | Batch size | Wall clock | Speedup |
|---|---|---:|---:|
| 481 sequences (`RefPlantNLR`) | `5000` (upstream default) | 776 s | 1.00x |
| 481 sequences (`RefPlantNLR`) | `64` (auto-sized) | 197 s | 3.93x |
| 6,000 sequences | `5000` (upstream default) | 1,368 s | 1.00x |
| 6,000 sequences | `250` (auto-sized) | 437 s | 3.13x |

Even at 6,000 sequences, where the upstream default already creates two batches, smaller batches cut runtime more than threefold by keeping all worker slots balanced instead of waiting on long-running tail batches.

### Automatic calculation in `NLRtracker.sh`

```text
cores  = physical_cores()           (/proc/cpuinfo -> lscpu -> sysctl hw.physicalcpu -> getconf)
slots  = max(1, cores / cpus_per_task)
target = clamp(n_sequences / (2 * slots), min = 16, max = 5000)
```

- **Physical cores rather than logical threads:** HMMER (`hmmsearch`, `hmmscan`) and `pfsearch` are memory-bandwidth and AVX-512 bound, so counting simultaneous multithreading pairs over-allocates work and reduces throughput.
- **Two work units per slot:** Overlapping two batches per slot hides slow batches behind concurrent work.
- **Bounds:** `16` prevents container startup overhead from dominating tiny FASTAs, and `5000` avoids exceeding the upstream default on massive proteomes. Override with `-b INT` when comparing runs.

## Local compute versus InterPro Matches API

- **Local compute (`-L` or unpublished sequences):** Every sequence executes in `SCAN_SEQUENCES` containers locally, so auto-sized batches deliver the 3–4x speedup measured above.
- **Published proteomes (default API lookup enabled):** InterProScan retrieves cached matches over HTTPS first and computes only unmatched proteins locally. Network latency controls wall-clock time unless many sequences are absent from the remote lookup.

## SMART batch-size sensitivity

Changing `-b` slightly shifts `SMART` hit counts (roughly `0.1%` of total InterProScan rows at proteome scale) because `PREFILTER_SMART` upstream does not pass `-Z` database normalisation across chunks, making E-values sensitive to chunk size. `SMART` supplies no NLR-defining subclass marker, so classification barely moves; for strict bit-for-bit reproducibility across runs, pass both `-L` and an explicit `-b INT`.
