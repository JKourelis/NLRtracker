## 02_domains.R -- turn per-signature annotations into per-class domain intervals.
##
## Replaces the five near-identical dedup loops in the legacy script with one generic
## merge driven by the merge_group / precedence columns of module/nlr_classes.tsv.
##
## Merging runs in ordered stages, because a class can take part in more than one.
## PLOOP is merged into NB-ARC first, and only the surviving PLOOP is then compared
## with integrated domains to produce PLOOP-OTHER -- the same order as the legacy code.

## Stage definitions.
##   group     : merge_group in nlr_classes.tsv, or an explicit class vector
##   mode      : "all"      merge every member class in one pass, highest precedence wins
##               "pairwise" merge each member against `against` separately, so two
##                          specific domains can never be joined through a shared OTHER
##   relabel   : optional function(classes_present) -> class name
MERGE_STAGES <- list(
  list(name = "nbarc", group = "nbarc", mode = "all"),

  list(name = "ploop_other",
       classes = c("PLOOP", "OTHER"), mode = "all",
       relabel = function(present) {
         if (all(c("PLOOP", "OTHER") %in% present)) "PLOOP-OTHER"
         else if ("PLOOP" %in% present) "PLOOP" else "OTHER"
       }),

  list(name = "integrated", group = "id", mode = "pairwise", against = "OTHER"),

  list(name = "rpw8", group = "rpw8", mode = "all")
)

## Merge overlapping [start, end] intervals within each seqname, keeping the set of
## classes that contributed. Touching intervals (end + 1 == next start) are treated as
## separate, matching the legacy behaviour.
merge_intervals <- function(df) {
  if (nrow(df) == 0) return(df[0, c("seqname", "domain_class", "start", "end")])
  df %>%
    arrange(seqname, start, end) %>%
    group_by(seqname) %>%
    mutate(.grp = cumsum(c(TRUE, start[-1] > cummax(end)[-n()]))) %>%
    group_by(seqname, .grp) %>%
    summarise(start = min(start), end = max(end),
              .classes = list(unique(domain_class)), .groups = "drop") %>%
    select(-.grp)
}

resolve_by_precedence <- function(classes, cfg) {
  p <- cfg$prec_of[classes]
  p[is.na(p)] <- -1L
  classes[which.max(p)]
}

apply_stage <- function(domains, stage, cfg) {
  members <- if (!is.null(stage$classes)) {
    stage$classes
  } else {
    cfg$classes$domain_class[cfg$classes$merge_group == stage$group]
  }
  members <- intersect(members, unique(domains$domain_class))
  if (length(members) < 1) return(domains)

  if (identical(stage$mode, "pairwise")) {
    against <- stage$against
    specific <- setdiff(members, against)
    if (!(against %in% unique(domains$domain_class)) || length(specific) == 0) {
      return(domains)
    }
    ## Each specific class is merged against `against` on its own. The `against`
    ## rows are consumed by whichever pass claims them, so process by descending
    ## precedence and remove as we go.
    specific <- specific[order(-cfg$prec_of[specific])]
    rest <- domains %>% filter(!domain_class %in% c(specific, against))
    pool <- domains %>% filter(domain_class %in% c(specific, against))
    kept <- list()
    for (cl in specific) {
      pair <- pool %>% filter(domain_class %in% c(cl, against))
      if (!nrow(pair)) next
      merged <- merge_intervals(pair) %>%
        mutate(domain_class = purrr::map_chr(.classes, resolve_by_precedence, cfg = cfg)) %>%
        select(seqname, domain_class, start, end)
      kept[[cl]] <- merged %>% filter(domain_class == cl)
      ## anything that stayed `against` remains available to the next class
      pool <- bind_rows(
        pool %>% filter(!domain_class %in% c(cl, against)),
        merged %>% filter(domain_class == against)
      )
    }
    return(bind_rows(rest, bind_rows(kept), pool %>% filter(domain_class == against)))
  }

  rest <- domains %>% filter(!domain_class %in% members)
  merged <- domains %>%
    filter(domain_class %in% members) %>%
    merge_intervals() %>%
    mutate(domain_class = if (is.null(stage$relabel)) {
      purrr::map_chr(.classes, resolve_by_precedence, cfg = cfg)
    } else {
      purrr::map_chr(.classes, stage$relabel)
    }) %>%
    select(seqname, domain_class, start, end)
  bind_rows(rest, merged)
}

## The NB-ARC motifs, grouped into runs. A decrease in the canonical motif index
## marks a new NB-ARC, which is how tandem NB-ARCs stay separate.
nbarc_motif_intervals <- function(annotation) {
  m <- annotation %>%
    filter(source == "fimo", Name %in% names(NBARC_MOTIF_ORDER)) %>%
    mutate(idx = unname(NBARC_MOTIF_ORDER[Name])) %>%
    arrange(seqname, start)
  if (nrow(m) == 0) {
    return(tibble::tibble(seqname = character(), domain_class = character(),
                          start = integer(), end = integer()))
  }
  m %>%
    group_by(seqname) %>%
    mutate(run = cumsum(c(TRUE, diff(idx) <= 0))) %>%
    group_by(seqname, run) %>%
    summarise(start = min(start), end = max(end), .groups = "drop") %>%
    transmute(seqname, domain_class = "nbarc-motif", start, end)
}

## Motif-derived evidence classes, kept as points rather than merged regions.
motif_evidence <- function(annotation, cfg) {
  annotation %>%
    filter(source == "fimo", !is.na(domain_class)) %>%
    filter(domain_class %in% c("cc-motif", "rnbs-d", "other-motif", "linker-MHD")) %>%
    ## linker and MHD only count above a higher score
    mutate(score_num = suppressWarnings(as.numeric(score))) %>%
    filter(!(Name %in% names(MOTIF_HIGH_SCORE)) |
             score_num >= unname(MOTIF_HIGH_SCORE[Name])) %>%
    transmute(seqname, domain_class, start, end)
}

build_domains <- function(annotation, cfg) {
  ## Signature-derived intervals, dropping anything with no class assignment.
  sig_domains <- annotation %>%
    filter(source != "fimo", !is.na(domain_class)) %>%
    transmute(seqname, domain_class, start, end)

  domains <- bind_rows(
    sig_domains,
    nbarc_motif_intervals(annotation),
    motif_evidence(annotation, cfg)
  )

  if (nrow(domains) == 0) {
    stop("no NLR-relevant domains were found. Check that InterProScan ran the ",
         "expected analyses.", call. = FALSE)
  }

  ## Collapse repeats of the same class before the cross-class stages.
  domains <- domains %>%
    group_by(domain_class) %>%
    group_modify(~ merge_intervals(mutate(.x, domain_class = .y$domain_class)) %>%
                   select(seqname, start, end)) %>%
    ungroup() %>%
    select(seqname, domain_class, start, end)

  for (stage in MERGE_STAGES) domains <- apply_stage(domains, stage, cfg)

  domains %>%
    filter(!is.na(domain_class)) %>%
    arrange(seqname, start, end)
}

## Classes that only existed to extend a region are removed once merging is done.
drop_scaffold_classes <- function(domains, cfg) {
  domains %>% filter(!domain_class %in% cfg$drop_after)
}
