#!/usr/bin/env bash
#
# NLRtracker - extract and functionally annotate NLRs from plant proteomes
#
# Runs InterProScan 6 (via Nextflow + containers), FIMO and HMMER, then classifies
# NLRs in R.  See README.md for requirements and docs/ for the measurements
# behind the parallelisation defaults.
#
set -euo pipefail

VERSION="2.0.0"
INTERPROSCAN_REVISION="6.0.1"
INTERPRO_VERSION="109.0"

# The seven analyses NLRtracker classifies on. Renamed for InterProScan 6:
#   Gene3D -> cathgene3d, ProSiteProfiles -> prosite-profiles.
APPLICATIONS="pfam,cathgene3d,superfamily,prints,smart,cdd,prosite-profiles"

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
module_dir="${here}/module"

# ------------------------------------------------------------------ defaults
fasta=""
outdir=""
seqtype="p"
interpro_result=""
fimo_result=""
hmmer_result=""
meme_xml="${module_dir}/meme.xml"
hmm_file="${module_dir}/abe3069_Data_S1.hmm"
entry_list=""
datadir="${INTERPROSCAN_DATADIR:-${HOME}/Bioinformatics/interproscan6/data}"
container_profile="docker"
cpus_per_task=2
batch_size=""          # empty means auto-size, see choose_batch_size()
max_workers=""         # empty means let Nextflow use all available processors
phylogeny=0
use_matches_api=1
keep_work=0

# ------------------------------------------------------------------ logging
log()  { printf '[%s] %s\n'  "$(date '+%H:%M:%S')" "$*"; }
warn() { printf '[%s] WARNING: %s\n' "$(date '+%H:%M:%S')" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n'   "$(date '+%H:%M:%S')" "$*" >&2; exit 1; }

usage() {
  cat <<EOM
NLRtracker ${VERSION}  (InterProScan ${INTERPROSCAN_REVISION}, InterPro ${INTERPRO_VERSION})

Usage: $(basename "$0") -s <proteins.fasta> -o <outdir> [options]

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
  -d FILE     InterPro entry list          [module/InterPro_<ver>.list]
  -D DIR      InterProScan data directory  [\$INTERPROSCAN_DATADIR or
                                            ~/Bioinformatics/interproscan6/data]

Execution:
  -t p|n      Sequence type                [p]
  -c INT      CPUs per InterProScan task   [${cpus_per_task}]
  -b INT      InterProScan batch size      [auto, see below]
  -w INT      Max concurrent tasks         [all available processors]
  -P NAME     Container profile: docker | singularity | apptainer  [docker]
  -p          Also classify phylogenetically against RefPlantNLR.
              Needs hmmalign, VeryFastTree (or FastTree) and the R package ape.
  -L          Disable the InterPro Matches API (fully local, reproducible)
  -k          Keep the Nextflow work directory
  -h          This help

Batch size is auto-sized from the input, because the InterProScan default of
5000 collapses to a single work unit on small inputs and leaves most of the
machine idle. Measured on 481 sequences: 777 s at the default against 198 s
auto-sized. See docs/performance.md.

Note: changing the batch size changes SMART results slightly, because
PREFILTER_SMART is not normalised with -Z upstream. The value used is recorded
in <outdir>/run_info.yaml.
EOM
  exit "${1:-2}"
}

while getopts ":s:o:i:f:X:m:x:d:D:t:c:b:w:P:pLkh" opt; do
  case "${opt}" in
    s) fasta="${OPTARG}" ;;
    o) outdir="${OPTARG}" ;;
    i) interpro_result="${OPTARG}" ;;
    f) fimo_result="${OPTARG}" ;;
    X) hmmer_result="${OPTARG}" ;;
    m) meme_xml="${OPTARG}" ;;
    x) hmm_file="${OPTARG}" ;;
    d) entry_list="${OPTARG}" ;;
    D) datadir="${OPTARG}" ;;
    t) seqtype="${OPTARG}" ;;
    c) cpus_per_task="${OPTARG}" ;;
    b) batch_size="${OPTARG}" ;;
    w) max_workers="${OPTARG}" ;;
    P) container_profile="${OPTARG}" ;;
    p) phylogeny=1 ;;
    L) use_matches_api=0 ;;
    k) keep_work=1 ;;
    h) usage 0 ;;
    :) die "option -${OPTARG} needs an argument" ;;
    \?) die "unknown option -${OPTARG}" ;;
  esac
done

# ------------------------------------------------------------------ validate
[[ -n "${fasta}"  ]] || { warn "-s is required"; usage; }
[[ -n "${outdir}" ]] || { warn "-o is required"; usage; }
[[ -f "${fasta}"  ]] || die "input FASTA not found: ${fasta}"
[[ "${seqtype}" == "p" || "${seqtype}" == "n" ]] || die "-t must be p or n (got '${seqtype}')"
[[ "${cpus_per_task}" =~ ^[0-9]+$ && "${cpus_per_task}" -ge 1 ]] \
  || die "-c must be a positive integer (got '${cpus_per_task}')"
[[ -z "${batch_size}"  || "${batch_size}"  =~ ^[0-9]+$ ]] || die "-b must be an integer"
[[ -z "${max_workers}" || "${max_workers}" =~ ^[0-9]+$ ]] || die "-w must be an integer"
case "${container_profile}" in docker|singularity|apptainer) ;; *)
  die "-P must be docker, singularity or apptainer (got '${container_profile}')" ;; esac

if [[ -z "${entry_list}" ]]; then
  entry_list="${module_dir}/InterPro_${INTERPRO_VERSION}.list"
fi
for f in "${entry_list}" "${module_dir}/nlr_signatures.tsv" \
         "${module_dir}/nlr_classes.tsv" "${module_dir}/nlr_subclass_rules.tsv"; do
  [[ -f "${f}" ]] || die "reference file missing: ${f}"
done
[[ -f "${meme_xml}" ]] || die "meme.xml not found: ${meme_xml}"

need() { command -v "$1" >/dev/null 2>&1 || die "$1 not found in PATH. See README.md"; }
need Rscript
[[ -n "${fimo_result}" ]] || need fimo
if [[ "${seqtype}" == "p" && -z "${hmmer_result}" ]]; then
  need hmmsearch
  [[ -f "${hmm_file}" ]] || die "HMM not found: ${hmm_file}"
fi
if [[ -z "${interpro_result}" ]]; then
  need nextflow
  [[ -d "${datadir}" ]] || die "InterProScan data directory not found: ${datadir}
Set it with -D, or with the INTERPROSCAN_DATADIR environment variable."
fi
# -p adds three dependencies. Check them now: the InterProScan step runs first and
# takes far longer, so a missing tree program must not surface at the end of it.
if (( phylogeny )); then
  need hmmalign
  command -v VeryFastTree >/dev/null 2>&1 || command -v FastTree >/dev/null 2>&1 \
    || die "-p needs VeryFastTree or FastTree in PATH. See README.md"
  Rscript -e 'if (!requireNamespace("ape", quietly = TRUE)) quit(status = 1)' \
    || die "-p needs the R package ape. Install it with install.packages(\"ape\")"
  [[ -f "${module_dir}/nlr_clades.tsv" ]] || die "reference file missing: ${module_dir}/nlr_clades.tsv"
  [[ -f "${module_dir}/reference/NB-ARC_region.hmm" ]] \
    || die "reference file missing: ${module_dir}/reference/NB-ARC_region.hmm"
fi

# ------------------------------------------------------------------ layout
mkdir -p "${outdir}"/{figures,sequences,annotations,raw/logs}
outdir="$(cd "${outdir}" && pwd)"
fasta="$(cd "$(dirname "${fasta}")" && pwd)/$(basename "${fasta}")"

log_file="${outdir}/raw/logs/nlrtracker.log"
exec > >(tee -a "${log_file}") 2>&1

log "NLRtracker ${VERSION}"
log "input        ${fasta}"
log "output       ${outdir}"

# One-line-per-sequence FASTA. Headers are truncated at the first whitespace,
# which is what the downstream identifiers assume.
work_fasta="${outdir}/raw/input.normalised.fasta"
awk '{ if ($1 ~ /^>/) print "\n"$1; else printf "%s", $1 }' "${fasta}" \
  | sed -e '1d' > "${work_fasta}"
n_seq="$(grep -c '^>' "${work_fasta}" || true)"
[[ "${n_seq}" -gt 0 ]] || die "no sequences found in ${fasta}"
log "sequences    ${n_seq}"

# ------------------------------------------------------------------ batching
# Give the Nextflow scheduler more work units than it has slots. Concurrency
# itself is Nextflow's job and is deliberately not computed here.
choose_batch_size() {
  local n="$1" cpus="$2" cores slots target
  cores="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"
  # Physical cores, not SMT threads: the hmmsearch/pfsearch family is
  # AVX-512 and memory-bandwidth bound, where SMT returns far less than 2x.
  if command -v lscpu >/dev/null 2>&1; then
    local tpc
    tpc="$(lscpu 2>/dev/null | awk -F: '/^Thread\(s\) per core/ {gsub(/ /,"",$2); print $2}')"
    [[ "${tpc}" =~ ^[0-9]+$ && "${tpc}" -ge 1 ]] && cores=$(( cores / tpc ))
  elif command -v sysctl >/dev/null 2>&1; then
    # macOS reports physical cores directly.
    local phys
    phys="$(sysctl -n hw.physicalcpu 2>/dev/null || true)"
    [[ "${phys}" =~ ^[0-9]+$ && "${phys}" -ge 1 ]] && cores="${phys}"
  fi
  (( cores < 1 )) && cores=1
  slots=$(( cores / cpus )); (( slots < 1 )) && slots=1
  target=$(( n / (2 * slots) ))
  (( target < 16 ))   && target=16      # container startup would dominate
  (( target > 5000 )) && target=5000    # never exceed the upstream default
  echo "${target}"
}

if [[ -z "${batch_size}" ]]; then
  batch_size="$(choose_batch_size "${n_seq}" "${cpus_per_task}")"
  log "batch size   ${batch_size} (auto)"
else
  log "batch size   ${batch_size} (set with -b)"
fi

# ------------------------------------------------------------------ 1. InterProScan
ips_tsv="${interpro_result}"
if [[ -z "${ips_tsv}" ]]; then
  log "running InterProScan ${INTERPROSCAN_REVISION} / InterPro ${INTERPRO_VERSION}"
  ips_out="${outdir}/raw/interproscan"
  mkdir -p "${ips_out}"

  # Numeric and boolean parameters MUST go through a config file. Nextflow
  # 26.x passes command-line params as String, and InterProScan calls
  # splitFasta(Path, int, boolean), so `--batchSize 250` crashes SPLIT_FASTA.
  run_config="${ips_out}/nlrtracker.config"
  {
    echo "// generated by NLRtracker ${VERSION} - do not edit"
    echo "params {"
    echo "    batchSize    = ${batch_size}"
    echo "    cpus         = ${cpus_per_task}"
    [[ -n "${max_workers}" ]] && echo "    maxWorkers   = ${max_workers}"
    [[ "${use_matches_api}" -eq 0 ]] && echo "    noMatchesApi = true"
    echo "}"
  } > "${run_config}"

  nf_args=(run ebi-pf-team/interproscan6 -r "${INTERPROSCAN_REVISION}"
           -profile "${container_profile}"
           -c "${run_config}"
           -work-dir "${ips_out}/work"
           --input "${work_fasta}"
           --outdir "${ips_out}"
           --datadir "${datadir}"
           --interpro "${INTERPRO_VERSION}"
           --applications "${APPLICATIONS}"
           --formats tsv,json
           -with-report "${ips_out}/nextflow_report.html"
           -with-trace  "${ips_out}/trace.txt"
           -resume)
  [[ "${seqtype}" == "n" ]] && nf_args+=(--nucleic)

  nextflow -log "${ips_out}/nextflow.log" "${nf_args[@]}"

  ips_tsv="$(find "${ips_out}" -maxdepth 1 -name '*.tsv' -type f | head -1)"
  [[ -n "${ips_tsv}" ]] || die "InterProScan produced no TSV in ${ips_out}"
  log "interproscan ${ips_tsv} ($(wc -l < "${ips_tsv}") rows)"
else
  log "reusing InterProScan result ${ips_tsv}"
  [[ -f "${ips_tsv}" ]] || die "InterProScan result not found: ${ips_tsv}"
fi

# ------------------------------------------------------------------ 2. FIMO
if [[ -z "${fimo_result}" ]]; then
  log "running FIMO"
  fimo --oc "${outdir}/raw/fimo" "${meme_xml}" "${work_fasta}" >/dev/null 2>&1 \
    || die "FIMO failed. See ${outdir}/raw/fimo"
  fimo_result="${outdir}/raw/fimo/fimo.gff"
  [[ -f "${fimo_result}" ]] || die "FIMO produced no fimo.gff"
  log "fimo         ${fimo_result} ($(wc -l < "${fimo_result}") rows)"
else
  log "reusing FIMO result ${fimo_result}"
  [[ -f "${fimo_result}" ]] || die "FIMO result not found: ${fimo_result}"
fi

# ------------------------------------------------------------------ 3. HMMER
if [[ "${seqtype}" == "p" ]]; then
  if [[ -z "${hmmer_result}" ]]; then
    log "running hmmsearch (C-JID)"
    hmmer_result="${outdir}/raw/hmmer/CJID.domtblout"
    mkdir -p "$(dirname "${hmmer_result}")"
    hmmsearch --domtblout "${hmmer_result}" "${hmm_file}" "${work_fasta}" >/dev/null \
      || die "hmmsearch failed"
    log "hmmer        ${hmmer_result} ($(grep -vc '^#' "${hmmer_result}" || true) rows)"
  else
    log "reusing HMMER result ${hmmer_result}"
    [[ -f "${hmmer_result}" ]] || die "HMMER result not found: ${hmmer_result}"
  fi
else
  warn "nucleotide input: C-JID detection is skipped (protein input only)"
  hmmer_result=""
fi

# ------------------------------------------------------------------ 3b. NB-ARC HMM
# The locked NB-ARC_region.hmm defines the NB-ARC region boundaries. It is a strict superset
# of signature-based detection and recovers partial copies that no member database finds --
# RRS1's second NB-ARC among them. See docs/nbarc-hmm.md.
nbarc_hmm_out=""
nbarc_hmm="${module_dir}/reference/NB-ARC_region.hmm"
if [[ -f "${nbarc_hmm}" ]] && command -v hmmsearch >/dev/null 2>&1; then
  log "running hmmsearch (NB-ARC region)"
  nbarc_hmm_out="${outdir}/raw/hmmer/NBARC.domtblout"
  mkdir -p "$(dirname "${nbarc_hmm_out}")"
  # -E is permissive on purpose; the calibrated score thresholds are applied in R
  hmmsearch --cpu "${cpus_per_task}" -E 1e-3 --domtblout "${nbarc_hmm_out}" \
    "${nbarc_hmm}" "${work_fasta}" >/dev/null \
    || die "NB-ARC hmmsearch failed"
  log "nbarc hmm    ${nbarc_hmm_out} ($(grep -vc '^#' "${nbarc_hmm_out}" || true) domains)"
else
  warn "NB-ARC HMM not available; falling back to signature-defined NB-ARC boundaries"
fi

# ------------------------------------------------------------------ 4. classify
log "classifying"
Rscript "${module_dir}/NLRtracker.R" \
  --interproscan   "${ips_tsv}" \
  --fimo           "${fimo_result}" \
  --hmmer          "${hmmer_result}" \
  --nbarc-hmm      "${nbarc_hmm_out}" \
  --fasta          "${work_fasta}" \
  --entry-list     "${entry_list}" \
  --signatures     "${module_dir}/nlr_signatures.tsv" \
  --classes        "${module_dir}/nlr_classes.tsv" \
  --subclass-rules "${module_dir}/nlr_subclass_rules.tsv" \
  --outdir         "${outdir}" \
  --seqtype        "${seqtype}" \
  --version        "${VERSION}" \
  --interproscan-version "${INTERPROSCAN_REVISION}" \
  --interpro-version     "${INTERPRO_VERSION}" \
  --applications   "${APPLICATIONS}" \
  --batch-size     "${batch_size}" \
  --phylogeny      "$( ((phylogeny)) && echo true || echo false )" \
  --input-name     "$(basename "${fasta}")"

# ------------------------------------------------------------------ tidy
if [[ "${keep_work}" -eq 0 && -d "${outdir}/raw/interproscan/work" ]]; then
  log "removing the Nextflow work directory (-k keeps it)"
  rm -rf "${outdir}/raw/interproscan/work"
fi
rm -f "${outdir}/raw/input.normalised.fasta.tmp"

log "done"
log "  classification  ${outdir}/NLRtracker.tsv"
log "  domains         ${outdir}/Domains.tsv"
log "  summary         ${outdir}/summary.tsv"
log "  figures         ${outdir}/figures/"
log "  sequences       ${outdir}/sequences/"
log "  annotations     ${outdir}/annotations/"
if ((phylogeny)); then log "  phylogeny       ${outdir}/phylogeny/"; fi
log "  provenance      ${outdir}/run_info.yaml"
