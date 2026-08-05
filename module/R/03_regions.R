## 03_regions.R -- promote domain classes to arbitrated, non-overlapping regions.
##
## NB-ARC is the anchor and is never trimmed. Arbitration then runs in a fixed order:
##
##   1. NB-ARC                       fixed
##   2. TIR, CC, RPW8, R1            if overlapping NB-ARC, end = NB-ARC start - 1
##   3. LRR                          if overlapping NB-ARC, start = NB-ARC end + 1
##   4. C-JID                        if overlapping LRR,    start = LRR end + CJID_LRR_OFFSET
##   5. ID                           one region per integrated-domain interval
##
## A region trimmed out of existence (end < start) is dropped and reported, never
## discarded silently: a TIR lying entirely inside an NB-ARC call is a signal worth
## seeing.
##
## The WHD signatures are deliberately NOT N-terminal regions. PF23282 and PF23559 sit
## inside the NB-ARC and are lineage markers, not domain markers -- 8 of 19 CCG10-NLRs
## carry PF23559 while none carries any true CC signature. They reach the classifier as
## subclass evidence instead. See docs/signatures.md section 3.

## Anchors are chosen per sequence: for tandem NB-ARCs each N-terminal or LRR region is
## trimmed against the nearest relevant NB-ARC rather than the first one.
trim_against <- function(regions, target_regions, anchor_region, side) {
  anchors <- regions %>% filter(region %in% anchor_region) %>%
    select(seqname, a_start = start, a_end = end)
  if (nrow(anchors) == 0) return(list(regions = regions, dropped = regions[0, ]))

  ## An anchor is never a target, so a multi-anchor set cannot trim itself.
  target_regions <- setdiff(target_regions, anchor_region)
  targets <- regions %>% filter(region %in% target_regions)
  others  <- regions %>% filter(!region %in% target_regions)
  if (nrow(targets) == 0) return(list(regions = regions, dropped = regions[0, ]))

  joined <- targets %>%
    mutate(.row = row_number()) %>%
    left_join(anchors, by = "seqname", relationship = "many-to-many")

  ## Only the target rows that actually overlap an anchor are adjusted. This is kept as its own
  ## step because `summarise()` on a zero-row input still evaluates the expression once to infer
  ## the column type, and `min()`/`max()` of nothing warns. Nothing to trim is the normal case for
  ## a sequence with no anchor overlap, so it must not produce a warning.
  overlapping <- joined %>% filter(!is.na(a_start), a_start <= end, a_end >= start)

  if (side == "before") {
    ## Truncate the target so it ends before the nearest overlapping anchor starts.
    adj <- if (nrow(overlapping) == 0) tibble::tibble(.row = integer(), new_end = integer()) else
      overlapping %>% group_by(.row) %>% summarise(new_end = min(a_start) - 1L, .groups = "drop")
    targets <- targets %>%
      mutate(.row = row_number()) %>%
      left_join(adj, by = ".row") %>%
      mutate(trimmed = !is.na(new_end),
             end = dplyr::if_else(trimmed, new_end, end)) %>%
      select(-new_end)
  } else {
    ## Truncate the target so it starts after the nearest overlapping anchor ends.
    adj <- if (nrow(overlapping) == 0) tibble::tibble(.row = integer(), new_start = integer()) else
      overlapping %>% group_by(.row) %>% summarise(new_start = max(a_end) + 1L, .groups = "drop")
    targets <- targets %>%
      mutate(.row = row_number()) %>%
      left_join(adj, by = ".row") %>%
      mutate(trimmed = !is.na(new_start),
             start = dplyr::if_else(trimmed, new_start, start)) %>%
      select(-new_start)
  }

  dropped <- targets %>% filter(end < start)
  kept    <- targets %>% filter(end >= start) %>% select(-.row, -trimmed)
  list(regions = bind_rows(others, kept),
       dropped = dropped %>% select(seqname, region, start, end))
}

## Interval difference: remove every canonical region from the intervals of
## `target_region`, splitting them where necessary. Fragments shorter than
## `min_len` are discarded as annotation noise rather than kept as slivers.
subtract_canonical <- function(regions, target_region, min_len = 10L) {
  targets <- regions %>% filter(region == target_region)
  others  <- regions %>% filter(region != target_region)
  if (nrow(targets) == 0 || nrow(others) == 0) {
    return(list(regions = regions, dropped = regions[0, ]))
  }

  mask <- others %>% select(seqname, m_start = start, m_end = end)

  kept <- targets %>%
    mutate(.row = row_number()) %>%
    group_split(.row) %>%
    purrr::map_dfr(function(one) {
      ms <- mask %>% filter(seqname == one$seqname,
                            m_start <= one$end, m_end >= one$start) %>%
        arrange(m_start)
      if (nrow(ms) == 0) return(one)
      pieces <- list()
      cursor <- one$start
      for (i in seq_len(nrow(ms))) {
        if (ms$m_start[i] > cursor) {
          pieces[[length(pieces) + 1]] <- c(cursor, ms$m_start[i] - 1L)
        }
        cursor <- max(cursor, ms$m_end[i] + 1L)
      }
      if (cursor <= one$end) pieces[[length(pieces) + 1]] <- c(cursor, one$end)
      if (!length(pieces)) return(one[0, ])
      purrr::map_dfr(pieces, function(p) {
        out <- one
        out$start <- as.integer(p[1]); out$end <- as.integer(p[2])
        out
      })
    })

  kept_ok <- kept %>% filter(end - start + 1L >= min_len) %>% select(-any_of(".row"))
  lost <- targets %>%
    anti_join(kept_ok %>% distinct(seqname), by = "seqname") %>%
    select(seqname, region, start, end)

  list(regions = bind_rows(others, kept_ok),
       dropped = lost)
}

build_regions <- function(domains, cfg, offset = CJID_LRR_OFFSET,
                          nbarc_hmm = NULL, annotation = NULL,
                          require_signature_support = TRUE) {
  regions <- domains %>%
    mutate(region = unname(cfg$region_of[domain_class])) %>%
    filter(!is.na(region), region != "-") %>%
    select(seqname, region, domain_class, start, end)

  ## Several classes can feed one region (for example every integrated domain feeds
  ## ID). Merge overlapping intervals per region so a region is a single interval
  ## wherever the underlying domains touch.
  ##
  ## `domain_class` is carried through the merge unchanged so that `.classes` collects
  ## the real contributing classes. An earlier version overwrote it with the region
  ## name before merging, which made `contributors` always equal the region name --
  ## useless for reporting, and it hid the motif-only bug fixed just below.
  regions <- regions %>%
    group_by(region) %>%
    group_modify(~ {
      merge_intervals(.x) %>%
        mutate(contributors = purrr::map_chr(.classes,
                                             ~ paste(sort(unique(.x)), collapse = ","))) %>%
        select(seqname, start, end, contributors)
    }) %>%
    ungroup() %>%
    select(seqname, region, start, end, contributors)

  ## Motif evidence may EXTEND a region but must never CREATE one.
  ##
  ## `nbarc_motif_intervals()` groups NB-ARC motifs into runs, and a single isolated
  ## motif far from the real NB-ARC forms its own run -- which then became its own
  ## NB-ARC region. Every 15 aa "NB-ARC region" in RefPlantNLR traced to exactly that:
  ## one stray RNBS-C, Kin-2 or NB-ARC-motif-12 hit, in several cases sitting inside an
  ## LRR repeat region (Pi2, Pi9, PigmR, PigmS, Piz-t, TmPm3, CD3.1, Rpp4C4).
  ##
  ## This is the VALIDATOR rule from InterPro_wrapper's region_definitions.tsv.
  ##
  ## Deliberately NOT a length filter: partial NB-ARC domains are real biology
  ## (Bph19, Pib and RGA2 carry genuine ~200 aa partial duplications), so a minimum
  ## length would discard those along with the artefacts. Evidence type is the
  ## discriminator, not size.
  motif_only_classes <- cfg$classes$domain_class[
    cfg$classes$drop_after & cfg$classes$region != "-"]
  region_motif_only <- function(contrib) {
    cls <- strsplit(contrib, ",", fixed = TRUE)[[1]]
    length(cls) > 0 && all(cls %in% motif_only_classes)
  }
  motif_only <- regions %>%
    filter(purrr::map_lgl(contributors, region_motif_only))
  regions <- regions %>%
    filter(!purrr::map_lgl(contributors, region_motif_only))

  dropped <- list()
  if (nrow(motif_only) > 0) {
    dropped[["region_from_motif_evidence_only"]] <- motif_only %>%
      select(seqname, region, start, end)
  }

  ## The HMM defines the NB-ARC boundaries where it has a hit. This happens BEFORE
  ## arbitration, so the N-terminal/LRR/C-JID trimming and the ID subtraction all work
  ## against the HMM boundaries rather than the signature ones.
  hmm_stats <- NULL
  if (!is.null(nbarc_hmm)) {
    hs <- apply_hmm_nbarc(regions, nbarc_hmm, annotation, require_signature_support)
    regions <- hs$regions
    hmm_stats <- hs
    if (!is.null(hs$rejected) && nrow(hs$rejected) > 0) {
      dropped[["nbarc_hmm_hit_rejected"]] <- hs$rejected %>%
        transmute(seqname, region = "NB-ARC", start = NA_integer_, end = NA_integer_)
    }
  }

  ## A bare P-loop overlapping a real NB-ARC is the same nucleotide-binding core with weaker
  ## evidence, and the class precedence already says NBARC beats PLOOP. Both are arbitration
  ## anchors, so neither trims the other -- drop the PLOOP instead. The HMM step can reveal
  ## this by widening the NB-ARC over a PLOOP that previously sat clear of it (seen on Pb1).
  ploop_shadowed <- regions %>% filter(region == "PLOOP") %>%
    inner_join(regions %>% filter(region == "NB-ARC") %>%
                 select(seqname, n_start = start, n_end = end),
               by = "seqname", relationship = "many-to-many") %>%
    filter(n_start <= end, n_end >= start) %>%
    distinct(seqname, region, start, end)
  if (nrow(ploop_shadowed) > 0) {
    regions <- regions %>% anti_join(ploop_shadowed, by = c("seqname", "region", "start", "end"))
    dropped[["ploop_inside_nbarc"]] <- ploop_shadowed
  }

  ## 2. N-terminal regions end before the NB-ARC.
  step <- trim_against(regions, REGION_NTERM, REGION_ANCHOR, "before")
  regions <- step$regions
  dropped[["n_terminal_inside_nbarc"]] <- step$dropped

  ## 3. LRR starts after the NB-ARC.
  step <- trim_against(regions, REGION_AFTER, REGION_ANCHOR, "after")
  regions <- step$regions
  dropped[["lrr_inside_nbarc"]] <- step$dropped

  ## 3b. Integrated domains yield to every canonical region.
  ##
  ## An ID can sit either side of the NB-ARC, so a one-sided truncation is wrong here:
  ## subtract the canonical intervals instead, which can split one ID into pieces.
  ## Observed cases are small (Pik-1's integrated HMA overlaps the NB-ARC start by a
  ## few residues; RPP2A and WRKY19 have an ID overlapping the TIR), so subtraction
  ## keeps the domain essentially intact while restoring disjointness.
  step <- subtract_canonical(regions, REGION_ID)
  regions <- step$regions
  dropped[["id_inside_canonical_region"]] <- step$dropped

  ## 4. C-JID starts after the LRR.
  if (offset != 1L) {
    ## honour a configured offset by shifting the LRR anchor end
    lrr_shift <- regions %>% filter(region == REGION_AFTER) %>%
      mutate(end = end + (offset - 1L))
    tmp <- bind_rows(regions %>% filter(region != REGION_AFTER), lrr_shift)
    step <- trim_against(tmp, REGION_TERMINAL, REGION_AFTER, "after")
    regions <- bind_rows(step$regions %>% filter(region != REGION_AFTER),
                         regions %>% filter(region == REGION_AFTER))
  } else {
    step <- trim_against(regions, REGION_TERMINAL, REGION_AFTER, "after")
    regions <- step$regions
  }
  dropped[["cjid_inside_lrr"]] <- step$dropped

  list(regions = regions %>% arrange(seqname, start, end),
       dropped  = bind_rows(dropped, .id = "rule"),
       hmm      = hmm_stats)
}

## Assert the arbitration actually produced disjoint regions. Called from the driver
## so a logic regression fails loudly instead of shipping wrong coordinates.
check_regions_disjoint <- function(regions) {
  bad <- regions %>%
    arrange(seqname, start, end) %>%
    group_by(seqname) %>%
    filter(n() > 1) %>%
    mutate(prev_end = lag(end), prev_region = lag(region)) %>%
    filter(!is.na(prev_end), start <= prev_end) %>%
    ungroup()
  bad %>% select(seqname, region, start, end, prev_region, prev_end)
}

## ------------------------------------------------------------------ HMM-defined NB-ARC

## Replace signature-derived NB-ARC intervals with the HMM envelopes.
##
## Measured on a Nicotiana proteome: the HMM is a strict superset -- 272 hits against 249 signature
## regions, 23 HMM-only, and zero signature-only. Boundaries move about 12 aa, toward greater
## completeness. It also recovers partial copies that signatures miss entirely: RRS1's second
## NB-ARC at 1068-1187 (score 36) is invisible to every member database.
##
## Multiple NB-ARCs stay separate, because hmmsearch reports one domain per copy. Bph1/9-1
## keeps a partial copy matching model positions 15-197 and a complete one matching 12-377 --
## which part of the model a copy matches is far more informative than its length.
##
## `require_signature_support` decides whether an HMM-only hit may create an NB-ARC where no
## member-database signature agrees. Of 23 HMM-only candidates in a test proteome, 6 carried CDC48
## (PF02359/PF02933) or AAA+ (PF00004/PF17862/PS00674) signatures and are not NLRs at all, so
## the default is to require corroboration.
NBARC_NLR_SUPPORT <- c("PF00931", "PF23282", "PF23559", "PF25895", "PF18052", "cd14798",
                       "G3DSA:1.20.5.4130", "PF01582", "PF13676", "SM00255", "PS50104",
                       "PF05659", "PS51153", "PF20160", "G3DSA:1.10.8.430")
NBARC_EXCLUDE_SUPPORT <- c("PF00004", "PF17862", "PS00674", "PF02933", "PF02359")

apply_hmm_nbarc <- function(regions, nbarc_hmm, annotation,
                            require_signature_support = TRUE,
                            protein_ga = NBARC_HMM_PROTEIN_GA) {
  if (is.null(nbarc_hmm) || nrow(nbarc_hmm) == 0) {
    return(list(regions = regions, replaced = 0L, added = 0L, rejected = NULL))
  }
  sig_prots <- regions %>% filter(region == "NB-ARC") %>% distinct(seqname) %>% pull()

  ## HMM-only proteins must clear the protein-level cutoff and, by default, carry a
  ## corroborating NLR signature while carrying no AAA+/CDC48 signature.
  hmm_only <- setdiff(unique(nbarc_hmm$seqname), sig_prots)
  rejected <- NULL
  if (length(hmm_only)) {
    ev <- annotation %>% filter(seqname %in% hmm_only) %>% distinct(seqname, Signature)
    best <- nbarc_hmm %>% group_by(seqname) %>%
      summarise(best_score = max(score), .groups = "drop")
    verdict <- tibble::tibble(seqname = hmm_only) %>%
      left_join(best, by = "seqname") %>%
      mutate(
        has_nlr = seqname %in% ev$seqname[ev$Signature %in% NBARC_NLR_SUPPORT],
        has_aaa = seqname %in% ev$seqname[ev$Signature %in% NBARC_EXCLUDE_SUPPORT],
        accepted = best_score >= protein_ga &
          (!require_signature_support | (has_nlr & !has_aaa)))
    rejected <- verdict %>% filter(!accepted) %>%
      mutate(reason = dplyr::case_when(
        best_score < protein_ga ~ "below protein gathering threshold",
        has_aaa & !has_nlr      ~ "AAA+/CDC48 signatures, not an NLR",
        !has_nlr                ~ "no corroborating NLR signature",
        TRUE                    ~ "rejected")) %>%
      select(seqname, best_score, reason)
    keep_only <- verdict$seqname[verdict$accepted]
    nbarc_hmm <- nbarc_hmm %>% filter(seqname %in% c(sig_prots, keep_only))
  }

  ## Adjacent HMM envelopes can overlap: Bph1/9-1's partial copy comes back as 138-388 and
  ## its complete copy as 381-768, sharing 8 residues. Both copies are real and must be kept
  ## (partial duplications are genuine biology), so the earlier one is trimmed to end just
  ## before the next begins rather than either being merged away or discarded.
  new_nbarc <- nbarc_hmm %>%
    arrange(seqname, start, end) %>%
    group_by(seqname) %>%
    mutate(next_start = dplyr::lead(start),
           end = dplyr::if_else(!is.na(next_start) & next_start <= end,
                                next_start - 1L, end)) %>%
    ungroup() %>%
    filter(end > start) %>%
    transmute(seqname, region = "NB-ARC", start, end,
              contributors = "NBARC-hmm", hmm_score = score, completeness)

  kept <- regions %>% filter(!(region == "NB-ARC" & seqname %in% new_nbarc$seqname))
  list(regions = bind_rows(kept, new_nbarc) %>% arrange(seqname, start, end),
       replaced = length(intersect(new_nbarc$seqname, sig_prots)),
       added    = length(setdiff(new_nbarc$seqname, sig_prots)),
       rejected = rejected)
}
