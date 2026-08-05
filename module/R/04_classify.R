## 04_classify.R -- assign Status and Subclass (putative) from the domain architecture.
##
## Two evidence streams are kept separate and combined here:
##
##   regions   what the protein is built from, after arbitration (03_regions.R)
##   evidence  per-signature subclass hints from nlr_signatures.tsv, plus the
##             multi-signature rules in nlr_subclass_rules.tsv
##
## The distinction matters most for the NB-ARC winged-helix Pfams. PF23559 and PF23282
## are region components of the NB-ARC, never N-terminal regions, but they are the
## strongest lineage discriminators available: PF23559 covers 99% of CC-NLRs and PF23282
## 94% of TIR-NLRs, with no cross-talk. See docs/signatures.md.

## Ordered architecture string, e.g. "(CC)(NBARC)(LRR)".
architecture_string <- function(domains) {
  domains %>%
    arrange(seqname, start, end) %>%
    group_by(seqname) %>%
    summarise(Domain = paste0("(", domain_class, ")", collapse = ""), .groups = "drop")
}

## Which subclass hints does each sequence carry?
collect_evidence <- function(annotation, cfg) {
  ev <- cfg$signatures %>%
    filter(!is.na(subclass_evidence), subclass_evidence != "-") %>%
    select(signature, subclass_evidence)
  if (nrow(ev) == 0) {
    return(tibble::tibble(seqname = character(), evidence = character()))
  }
  annotation %>%
    inner_join(ev, by = c("Signature" = "signature")) %>%
    distinct(seqname, subclass_evidence) %>%
    rename(evidence = subclass_evidence)
}

has_evidence <- function(ev, seqname, tag) {
  seqname %in% ev$seqname[ev$evidence == tag]
}

## Multi-signature rules: require_any present AND forbid_any absent.
apply_subclass_rules <- function(annotation, cfg) {
  if (nrow(cfg$rules) == 0) {
    return(tibble::tibble(seqname = character(), rule = character(),
                          subclass = character(), qualifier = character()))
  }
  sig_by_seq <- annotation %>% distinct(seqname, Signature)
  purrr::map_dfr(seq_len(nrow(cfg$rules)), function(i) {
    r <- cfg$rules[i, ]
    req <- stringr::str_split(r$require_any, ",", simplify = TRUE) %>% as.character() %>% trimws()
    fbd <- if (is.na(r$forbid_any) || r$forbid_any == "-") character() else
      stringr::str_split(r$forbid_any, ",", simplify = TRUE) %>% as.character() %>% trimws()
    hit <- sig_by_seq %>% filter(Signature %in% req) %>% pull(seqname) %>% unique()
    bad <- sig_by_seq %>% filter(Signature %in% fbd) %>% pull(seqname) %>% unique()
    tibble::tibble(seqname = setdiff(hit, bad), rule = r$rule,
                   subclass = r$subclass, qualifier = r$qualifier)
  })
}

classify <- function(domains, annotation, cfg) {
  arch <- architecture_string(domains)
  ev   <- collect_evidence(annotation, cfg)
  rules <- apply_subclass_rules(annotation, cfg)

  ## Presence flags, computed on the architecture string so the semantics match the
  ## legacy behaviour exactly.
  d <- arch %>%
    mutate(
      has_nbarc   = grepl("\\(NBARC\\)", Domain),
      has_nacht   = grepl("\\(NACHT\\)|\\(nacht-winged-helix\\)|\\(FISNA\\)", Domain),
      has_ploop   = grepl("\\(PLOOP\\)", Domain),
      has_tir     = grepl("\\(TIR\\)", Domain),
      has_cc      = grepl("\\(CC\\)|\\(R1\\)", Domain),
      has_rpw8    = grepl("\\(RPW8\\)", Domain),
      has_cbln    = grepl("\\(CblN\\)", Domain),
      has_mlkl    = grepl("\\(MLKL\\)", Domain),
      has_cjid    = grepl("\\(CJID\\)", Domain),
      has_lrr     = grepl("\\(LRR\\)", Domain),
      has_disease = grepl("\\(disease\\)", Domain),
      has_ccmotif = grepl("\\(cc-motif\\)", Domain),
      has_rnbsd   = grepl("\\(rnbs-d\\)", Domain),
      has_othermotif = grepl("\\(other-motif\\)", Domain),
      has_linkermhd  = grepl("\\(linker-MHD\\)", Domain),
      has_nonplant   = grepl("\\(HET\\)|\\(NN\\)", Domain),
      has_deathdom   = grepl("\\(NN\\)|\\(CARD\\)|\\(PYD\\)|\\(BIR\\)", Domain),
      ev_cc      = purrr::map_lgl(seqname, has_evidence, ev = ev, tag = "CC-NLR"),
      ev_tir     = purrr::map_lgl(seqname, has_evidence, ev = ev, tag = "TIR-NLR"),
      ev_cc_sup  = purrr::map_lgl(seqname, has_evidence, ev = ev, tag = "CC-NLR:supporting"),
      ev_tir_sup = purrr::map_lgl(seqname, has_evidence, ev = ev, tag = "TIR-NLR:supporting")
    )

  ## ---------------------------------------------------------------- Status
  d <- d %>%
    mutate(Status = dplyr::case_when(
      has_nonplant & (has_nacht | has_nbarc | has_ploop) ~ "likely non-plant NLR",
      has_nacht                                          ~ "NACHT",
      has_deathdom & has_ploop                           ~ "NACHT",
      has_nbarc                                          ~ "NLR",
      (has_cc | has_rpw8 | has_tir) & has_ploop          ~ "degenerate NLR",
      (has_cc | has_rpw8 | has_tir | has_lrr) &
        (has_rnbsd | has_linkermhd)                      ~ "degenerate NLR",
      has_disease                                        ~ "degenerate NLR",
      has_cc                                             ~ "CCX",
      has_rpw8                                           ~ "RPW8",
      has_tir                                            ~ "TX",
      has_cjid                                           ~ "CJID",
      has_mlkl                                           ~ "MLKL",
      TRUE                                               ~ NA_character_
    ))

  ## ------------------------------------------------------- Subclass (putative)
  ## Ordered most to least specific. The gated CCG10 rule is applied after the TIR
  ## and RPW8 routes so a diagnostic N-terminal call always wins, and before the
  ## generic CC routes because CCG10-NLRs have no detectable CC domain at all.
  ccg10 <- rules %>% filter(subclass == "CCG10-NLR")
  d <- d %>%
    mutate(
      is_nlr = Status %in% c("NLR", "degenerate NLR"),
      ccg10_hit = seqname %in% ccg10$seqname,
      `Subclass (putative)` = dplyr::case_when(
        Status %in% c("TX", "CJID", "RPW8", "MLKL", "CCX", "NACHT") ~ Status,

        ## diagnostic N-terminal domains
        is_nlr & has_tir  & (has_othermotif | has_linkermhd)  ~ "TIR-NLR",
        is_nlr & has_cjid                                     ~ "TIR-NLR",
        is_nlr & has_tir                                      ~ "TIR-NLR",
        is_nlr & has_rpw8                                     ~ "CCR-NLR",
        is_nlr & has_cbln                                     ~ "CCR-NLR",

        ## diagnostic lineage evidence from the NB-ARC winged-helix Pfams
        is_nlr & ev_tir & !has_cc                             ~ "TIR-NLR",

        ## CCG10: an RPS2-like LRR with no detectable CC domain
        is_nlr & ccg10_hit                                    ~ "CCG10-NLR",

        is_nlr & has_cc                                       ~ "CC-NLR",
        is_nlr & has_ccmotif                                  ~ "CC-NLR",

        ## Lineage evidence without a detectable CC domain cannot separate CC-NLR
        ## from CCG10-NLR: PF23559 covers 343/347 CC-NLRs but also 8/19 CCG10-NLRs,
        ## and CCG10-NLRs have no CC domain by definition. Measured on RefPlantNLR,
        ## 7 of 23 such calls are in fact CCG10, so reporting a bare "CC-NLR" here
        ## would be wrong 30% of the time. Report the ambiguity instead.
        is_nlr & (ev_cc | ev_cc_sup) & !has_cc   ~ "CC-NLR or CCG10-NLR",
        is_nlr & ev_cc                                        ~ "CC-NLR",
        is_nlr & ev_cc_sup                                    ~ "CC-NLR",
        is_nlr & ev_tir_sup                                   ~ "TIR-NLR",
        is_nlr & has_rnbsd                    ~ "CC-NLR or CCR-NLR or CCG10-NLR",
        Domain %in% c("(TIR)(NBARC)", "(TIR)(PLOOP)")         ~ "TN (OTHER)",
        !is.na(Domain)                                        ~ "UNDETERMINED",
        TRUE                                                  ~ NA_character_
      ),
      ## How the subclass was reached, and how much to trust it.
      subclass_support = dplyr::case_when(
        !is_nlr                                        ~ NA_character_,
        ccg10_hit & `Subclass (putative)` == "CCG10-NLR" ~ "probable",
        has_cc | has_tir | has_rpw8 | has_cbln | has_cjid ~ "domain",
        ev_cc  | ev_tir                                ~ "lineage",
        ev_cc_sup | ev_tir_sup                         ~ "supporting",
        has_ccmotif | has_rnbsd                        ~ "motif",
        TRUE                                           ~ NA_character_
      ),
      `Subclass (putative)` = dplyr::if_else(
        !is.na(subclass_support) & subclass_support == "probable",
        paste0(`Subclass (putative)`, " (probable)"),
        `Subclass (putative)`)
    )

  ## Simplified architecture string, driven by the simple_code column.
  simple <- cfg$classes %>% filter(!is.na(simple_code), simple_code != "")
  d <- d %>% mutate(Simple = Domain)
  ## Drop evidence-only classes from the printed architecture first.
  for (cl in cfg$drop_after) {
    d$Simple <- gsub(paste0("\\(", gsub("([\\W])", "\\\\\\1", cl, perl = TRUE), "\\)"),
                     "", d$Simple, perl = TRUE)
    d$Domain <- gsub(paste0("\\(", gsub("([\\W])", "\\\\\\1", cl, perl = TRUE), "\\)"),
                     "", d$Domain, perl = TRUE)
  }
  for (i in seq_len(nrow(simple))) {
    cl <- simple$domain_class[i]
    d$Simple <- gsub(paste0("\\(", gsub("([\\W])", "\\\\\\1", cl, perl = TRUE), "\\)"),
                     simple$simple_code[i], d$Simple, perl = TRUE)
  }
  ## Collapse runs of identical LRR / integrated-domain letters.
  d$Simple <- gsub("(L)\\1+", "\\1", d$Simple)
  d$Simple <- gsub("(O)\\1+", "\\1", d$Simple)

  d %>%
    filter(!is.na(Status)) %>%
    select(seqname, Status, `Subclass (putative)`, subclass_support, Domain, Simple)
}
