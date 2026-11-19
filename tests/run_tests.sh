#!/usr/bin/env bash
#
# NLRtracker regression tests.
#
# These are positive controls, in the laboratory sense: the answers are already known, so any
# change in them means something broke. Run this before every release, and after any edit to
# module/nlr_signatures.tsv, module/nlr_classes.tsv or module/nlr_subclass_rules.tsv -- those
# are data files, and a signature added to improve one class can silently break another.
#
#   tests/run_tests.sh            # fast: sample_data only, no external tools needed
#   tests/run_tests.sh --full     # also re-checks the 481 published RefPlantNLR subclasses
#
# --full needs the RefPlantNLR fixtures under dev/. They are not shipped, because they are
# 3 MB of intermediate output; regenerate them with the command printed if they are absent.
#
set -uo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${here}"
export R_LIBS_USER="${R_LIBS_USER:-${HOME}/R/x86_64-pc-linux-gnu-library/4.6}"

FULL=0
[[ "${1:-}" == "--full" ]] && FULL=1

pass=0; fail=0
ok()   { printf '  \033[32mPASS\033[0m  %s\n' "$1"; pass=$((pass+1)); }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }
check() { # check <description> <expected> <actual>
  if [[ "$2" == "$3" ]]; then ok "$1 = $3"; else bad "$1: expected $2, got $3"; fi
}

tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT

# ---------------------------------------------------------------- 1. reference data integrity
echo
echo "1. Reference data"
lock="module/reference/NB-ARC_region.LOCK"
if [[ -f "${lock}" ]]; then
  want="$(awk '/^sha256_hmm/{print $2}' "${lock}")"
  got="$(sha256sum module/reference/NB-ARC_region.hmm | awk '{print $1}')"
  check "NB-ARC HMM checksum matches the LOCK file" "${want}" "${got}"
else
  bad "module/reference/NB-ARC_region.LOCK is missing"
fi
check "signature table rows" "231" "$(grep -vc '^#' module/nlr_signatures.tsv | awk '{print $1-1}')"
check "class table rows"     "41"  "$(grep -vc '^#' module/nlr_classes.tsv    | awk '{print $1-1}')"

# every signature must reference a class that exists, or regions vanish silently
orphan="$(awk -F'\t' 'NR>1 && !/^#/{print $6}' module/nlr_signatures.tsv | sort -u > "${tmp}/a"
          awk -F'\t' 'NR>1 && !/^#/{print $1}' module/nlr_classes.tsv    | sort -u > "${tmp}/b"
          comm -23 "${tmp}/a" "${tmp}/b" | tr '\n' ' ')"
check "no signature references an undefined class" "" "${orphan// /}"

# ---------------------------------------------------------------- 2. sample_data smoke test
echo
echo "2. Sample data (no external tools required)"
Rscript module/NLRtracker.R \
  --interproscan sample_data/sample.interproscan6.tsv \
  --fimo         sample_data/sample.fimo.gff \
  --hmmer        sample_data/sample.CJID.domtblout \
  --nbarc-hmm    sample_data/sample.NBARC.domtblout \
  --fasta        sample_data/sample.fasta \
  --entry-list   module/InterPro_111.0.list \
  --signatures   module/nlr_signatures.tsv \
  --classes      module/nlr_classes.tsv \
  --subclass-rules module/nlr_subclass_rules.tsv \
  --outdir "${tmp}/sample" --seqtype p --version test \
  --interproscan-version 6.0.3.0 --interpro-version 111.0 \
  --applications pfam --batch-size 16 --input-name sample > "${tmp}/sample.log" 2>&1
rc=$?
check "pipeline exits cleanly" "0" "${rc}"
if [[ ${rc} -ne 0 ]]; then sed -n '$p' "${tmp}/sample.log"; fi

exp() { awk -F'\t' -v k="$1" '$1==k{print $2}' sample_data/sample.expected.tsv; }
if [[ -f "${tmp}/sample/NLRtracker.tsv" ]]; then
  check "classified sequences" "$(exp n_classified)" "$(($(wc -l < "${tmp}/sample/NLRtracker.tsv")-1))"
  check "regions"              "$(exp n_regions)" \
        "$(awk -F'\t' 'NR>1 && $2=="REGION"' "${tmp}/sample/Domains.tsv" | wc -l | tr -d ' ')"
  # every declared output must exist
  for f in NLRtracker.tsv Domains.tsv summary.tsv run_info.yaml report.html \
           sequences/NLR.fasta sequences/NBARC.fasta annotations/NLR.gff3 figures/NLR.pdf; do
    [[ -s "${tmp}/sample/${f}" ]] && ok "wrote ${f}" || bad "missing or empty: ${f}"
  done
  # a warning-free run: warnings that are always present are warnings nobody reads
  w="$(grep -ciE '^warning|parsing issue' "${tmp}/sample.log" || true)"
  check "no R warnings" "0" "${w}"
fi

# ---------------------------------------------------------------- 3. region invariants
echo
echo "3. Region invariants"
if [[ -f "${tmp}/sample/Domains.tsv" ]]; then
  # arbitration must leave regions disjoint; the pipeline asserts this, so catch the message
  ov="$(grep -c 'residual region overlap' "${tmp}/sample.log" || true)"
  check "no residual region overlaps" "0" "${ov}"
  # motif-only regions must never be created (a single stray FIMO hit is not a domain)
  mo="$(awk -F'\t' 'NR>1 && $2=="REGION" && ($5-$4+1)<30' "${tmp}/sample/Domains.tsv" | wc -l | tr -d ' ')"
  check "no region under 30 aa" "0" "${mo}"
fi

# ---------------------------------------------------------------- 4. full RefPlantNLR check
if [[ ${FULL} -eq 1 ]]; then
  echo
  echo "4. RefPlantNLR: 481 published subclasses"
  D="${REFPLANTNLR_FIXTURE_DIR:-dev}"
  need=("${D}/RefPlantNLR.interproscan6-111.0.tsv" "${D}/RefPlantNLR.fimo.gff"
        "${D}/RefPlantNLR.CJID.domtblout" "${D}/RefPlantNLR.NBARC.domtblout"
        "${D}/RefPlantNLR.norm.fasta")
  missing=0
  for f in "${need[@]}"; do [[ -f "${f}" ]] || { missing=1; echo "    absent: ${f}"; }; done
  if [[ ${missing} -eq 1 ]]; then
    echo "    Regenerate with:"
    echo "      curl -sSL -o dev/RefPlantNLR.fasta \\"
    echo "        'https://journals.plos.org/plosbiology/article/file?type=supplementary&id=10.1371/journal.pbio.3001124.s013'"
    echo "      ./NLRtracker.sh -s dev/RefPlantNLR.fasta -o dev/out_refplantnlr"
    bad "RefPlantNLR fixtures absent - skipped"
  else
    Rscript module/NLRtracker.R \
      --interproscan "${D}/RefPlantNLR.interproscan6-111.0.tsv" --fimo "${D}/RefPlantNLR.fimo.gff" \
      --hmmer "${D}/RefPlantNLR.CJID.domtblout" --nbarc-hmm "${D}/RefPlantNLR.NBARC.domtblout" \
      --fasta "${D}/RefPlantNLR.norm.fasta" --entry-list module/InterPro_111.0.list \
      --signatures module/nlr_signatures.tsv --classes module/nlr_classes.tsv \
      --subclass-rules module/nlr_subclass_rules.tsv --outdir "${tmp}/ref" --seqtype p \
      --version test --interproscan-version 6.0.3.0 --interpro-version 111.0 \
      --applications pfam --batch-size 64 --input-name RefPlantNLR > "${tmp}/ref.log" 2>&1
    check "pipeline exits cleanly" "0" "$?"
    check "sequences classified" "481" "$(($(wc -l < "${tmp}/ref/NLRtracker.tsv")-1))"

    # The specific boundary and lineage cases established for this release.
    nb() { awk -F'\t' -v P="$1" 'NR>1 && $1==P && $2=="REGION" && $3=="NB-ARC"' \
             "${tmp}/ref/Domains.tsv" | wc -l | tr -d ' '; }
    check "RRS1-R keeps both NB-ARC domains (HMM recovers the partial)" "2" "$(nb RRS1-R)"
    check "Bph1/9-1 partial duplication not merged"                     "2" "$(nb Bph19-1)"
    check "Prv^1 tandem NB-ARC not merged"                              "2" "$(nb 'Prv^1')"
    check "Pi2 stray-motif region not created"                          "1" "$(nb Pi2)"
    check "Rpp4C4 split NB-ARC envelopes joined on motif evidence"      "1" "$(nb Rpp4C4)"
    check "CCG10-NLR (probable) calls"                                 "11" \
          "$(awk -F'\t' 'NR>1 && $3 ~ /^CCG10-NLR/' "${tmp}/ref/NLRtracker.tsv" | wc -l | tr -d ' ')"
    ov="$(grep -c 'residual region overlap' "${tmp}/ref.log" || true)"
    check "no residual region overlaps" "0" "${ov}"
    w="$(grep -ciE '^warning|parsing issue' "${tmp}/ref.log" || true)"
    check "no R warnings" "0" "${w}"

    # Agreement with the published subclasses, if Table S1 is available.
    # Set REFPLANTNLR_TABLE_S1 to the published Table S1 spreadsheet to enable this check.
    S1="${REFPLANTNLR_TABLE_S1:-}"
    if [[ -f "${S1}" ]] && command -v python3 >/dev/null 2>&1; then
      read -r ex am wr < <(python3 tests/score_refplantnlr.py "${tmp}/ref/NLRtracker.tsv" "${S1}")
      check "exact agreement with published subclass" "457" "${ex}"
      check "honestly ambiguous (contains the truth)" "23"  "${am}"
      check "wrong"                                   "1"   "${wr}"
    else
      echo "    Table S1 not found - skipped the published-subclass comparison"
    fi
  fi
fi

# ---------------------------------------------------------------- summary
echo
printf '%s\n' "-----------------------------------------"
printf '  %d passed, %d failed\n' "${pass}" "${fail}"
[[ ${fail} -eq 0 ]] || exit 1
