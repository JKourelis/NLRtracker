# Performance

The InterProScan step dominates every run. This page records what was measured, so the defaults
in `NLRtracker.sh` can be checked rather than trusted.

Hardware: 24 physical cores (48 logical, 2 threads per core), NVMe storage, Docker.

---

## Batch size is the only tuning parameter that matters

InterProScan 6 splits the input into batches, then fans each batch out across the requested
applications. Nextflow schedules the resulting tasks. Concurrency is therefore already handled
upstream, and `NLRtracker.sh` does not try to manage it. What the wrapper sets is the **granularity**
of the work units.

The upstream default is `batchSize = 5000`. That value is sensible for a large proteome and wrong
for a small input, because 481 sequences become a **single** batch. One batch means one task per
application, and most of the machine stays idle.

| Input | Batch size | Wall clock | Speedup |
|---|---:|---:|---:|
| 481 sequences (RefPlantNLR) | 5000 (upstream default) | 777 s | — |
| 481 sequences | 64 (auto-sized) | 198 s | **3.92×** |
| 6000 sequences | 5000 | 2 118 s | — |
| 6000 sequences | 250 (auto-sized) | 677 s | **3.13×** |

The speedup does not collapse as the input grows. At 6000 sequences the default already produces
two batches, yet auto-sizing still gives 3.13×. The cause is **load balance**, not task count: batches
finish at different times, and coarse batches leave long idle tails.

### How the wrapper chooses it

```
cores  = physical cores          (lscpu on Linux, sysctl hw.physicalcpu on macOS)
slots  = cores / cpus_per_task
target = n_sequences / (2 * slots)      clamped to [16, 5000]
```

Twice as many work units as slots, so that a slow batch is overlapped by others rather than
becoming the critical path. Physical cores rather than logical ones: the `hmmsearch`/`pfsearch`
family is AVX-512 and memory-bandwidth bound, where simultaneous multithreading returns far less
than 2×.

The floor of 16 exists because container startup would otherwise dominate. The ceiling of 5000
never exceeds the upstream default.

Set `-b` to override the calculation. The value used is recorded in `run_info.yaml`.

---

## What the speedup does and does not apply to

**It applies to local compute.** All four measurements above ran with the InterPro Matches API
disabled, so every sequence was computed locally. The executed-task lists confirm it: only
`SCAN_SEQUENCES` ran, with no `LOOKUP` or `SCAN_REMAINING` stage.

**It does not apply to a proteome that is already in the API.** With the API enabled — the default —
a published proteome returns mostly from lookup. Wall clock is then network-bound, and the batch
size matters much less.

So the 3–4× figure is the **new-genome** case: an unpublished assembly, or any run with `-L`. That
is also the case where runtime actually hurts, which is why the default is auto-sized rather than
left at 5000.

---

## A known upstream inconsistency: SMART and `-b`

Changing the batch size changes SMART results slightly.

`PREFILTER_SMART` is not normalised with `-Z` upstream, so its E-values depend on the size of the
database chunk it is given, and that chunk size follows the sub-batch split. The effect is small and
confined to SMART, which contributes no NLR-defining signature. Everything else is unaffected.

Two consequences:

- Keep the batch size constant across runs you intend to compare. `run_info.yaml` records it.
- For a fully reproducible run, pass `-b` explicitly as well as `-L`.

This has not been reported upstream yet.
