## 05_write.R -- all output writers.
##
## Output tree:
##   NLRtracker.tsv  Domains.tsv  summary.tsv  run_info.yaml
##   figures/      NLR.pdf, NLR-associated.pdf, by_subclass/*.pdf
##   sequences/    NLR.fasta/.lst, NLR-associated.fasta/.lst, NBARC*.fasta, <region>.fasta
##   annotations/  NLR.gff3, NLR-associated.gff3, iTOL.txt, iTOL_dedup.txt
##   raw/          interproscan/, fimo/, hmmer/, logs/

write_tsv_out <- function(x, path) {
  readr::write_tsv(x, path, na = "", escape = "none", progress = FALSE)
}

write_lines_out <- function(x, path) {
  readr::write_lines(x, path)
}

## FASTA in the same 2-line style the legacy script produced.
write_fasta <- function(df, path, id_col = "seqname", seq_col = "sequence") {
  if (nrow(df) == 0) {
    file.create(path)
    return(invisible(0L))
  }
  write_lines_out(as.vector(rbind(paste0(">", df[[id_col]]), df[[seq_col]])), path)
  invisible(nrow(df))
}

## Region-derived FASTA. Multiple copies of a region in one protein get the legacy
## "_(n)" suffix, which downstream tools already expect from NBARC.fasta.
region_fasta <- function(regions, aa, region_name) {
  r <- regions %>% filter(region == region_name)
  if (nrow(r) == 0) return(tibble::tibble(seqname = character(), sequence = character()))
  r %>%
    inner_join(aa, by = "seqname") %>%
    group_by(seqname) %>%
    arrange(seqname, start, .by_group = TRUE) %>%
    mutate(copy = row_number(), n_copy = n()) %>%
    ungroup() %>%
    mutate(label = dplyr::if_else(n_copy == 1, seqname,
                                  paste0(seqname, "_(", copy, ")")),
           sequence = stringr::str_sub(sequence, start, end)) %>%
    filter(nchar(sequence) > 0) %>%
    select(seqname = label, sequence)
}

## Identical sequences collapsed with a "//" separator, as NBARC_deduplicated.fasta did.
dedup_fasta <- function(fa) {
  if (nrow(fa) == 0) return(fa)
  fa %>%
    group_by(sequence) %>%
    arrange(seqname, .by_group = TRUE) %>%
    summarise(seqname = paste(seqname, collapse = "//"), .groups = "drop") %>%
    select(seqname, sequence)
}

## GFF3 built from the internal annotation table, so its semantics are unchanged by
## the switch to TSV input. Noisy per-repeat rows are suppressed.
write_gff3 <- function(annotation, regions, ids, cfg, path, source_label = "NLRtracker") {
  noise <- gff3_noise_signatures(cfg)
  sig_rows <- annotation %>%
    filter(seqname %in% ids, !Signature %in% noise) %>%
    transmute(seqname, source, feature = "protein_match", start, end,
              score = dplyr::coalesce(score, "."), strand, frame, attribute)

  region_rows <- regions %>%
    filter(seqname %in% ids) %>%
    transmute(seqname, source = source_label, feature = "Region", start, end,
              score = ".", strand = "+", frame = ".",
              attribute = paste0("Name=", region, ";region_name=", region,
                                 ";contributors=", contributors,
                                 ";NCBI Feature Key=Region"))

  bind_rows(region_rows, sig_rows) %>%
    arrange(seqname, start, end) %>%
    readr::write_tsv(path, col_names = FALSE, na = "", escape = "none", progress = FALSE)
}

## Domains.tsv -- the table refplantnlR::draw_nlr() consumes.
##
## Every integrated domain contributes both a REGION row named ID (one consistent
## handle for "integrated domain") and a DOMAIN row keeping its real name for labels.
build_domains_table <- function(regions, domains, annotation, classification, aa, cfg) {
  chain <- aa %>%
    filter(seqname %in% classification$seqname) %>%
    transmute(seqname, type = "CHAIN", description = "chain",
              start = 1L, end = nchar(sequence))

  region_rows <- regions %>%
    filter(seqname %in% classification$seqname) %>%
    transmute(seqname, type = "REGION", description = region, start, end)

  ## Specific integrated-domain labels, for the figure and the GFF3.
  id_classes <- cfg$classes$domain_class[cfg$classes$region == REGION_ID]
  domain_rows <- domains %>%
    filter(seqname %in% classification$seqname,
           domain_class %in% id_classes,
           !domain_class %in% cfg$drop_after) %>%
    transmute(seqname, type = "DOMAIN", description = domain_class, start, end)

  motif_rows <- annotation %>%
    filter(source == "fimo", seqname %in% classification$seqname,
           Name %in% names(NBARC_MOTIF_ORDER)) %>%
    transmute(seqname, type = "MOTIF", description = Name, start, end)

  ## Individual LRR repeats, kept as REPEAT rows so draw_nlr(label_repeats = TRUE) has
  ## something to draw. These are the per-repeat signature hits that the LRR region is
  ## assembled from, and they are the rows suppressed from the published GFF3 as noise.
  lrr_sigs <- cfg$signatures$signature[cfg$signatures$domain_class == "LRR"]
  repeat_rows <- annotation %>%
    filter(seqname %in% classification$seqname,
           Signature %in% lrr_sigs,
           ENTRY_TYPE %in% c("Repeat", "repeat") |
             stringr::str_detect(dplyr::coalesce(Name, ""), stringr::regex("repeat", TRUE))) %>%
    transmute(seqname, type = "REPEAT", description = Name, start, end)

  bind_rows(chain, region_rows, domain_rows, motif_rows, repeat_rows) %>%
    left_join(classification, by = "seqname") %>%
    left_join(aa, by = "seqname") %>%
    mutate(sequence = stringr::str_sub(sequence, start, end)) %>%
    arrange(seqname, start, end)
}

build_summary <- function(classification, regions, annotation, cfg) {
  by_class <- classification %>%
    count(Status, `Subclass (putative)`, subclass_support, name = "n") %>%
    arrange(desc(n)) %>%
    mutate(metric = "classification") %>%
    select(metric, Status, subclass = `Subclass (putative)`,
           support = subclass_support, n)

  by_region <- regions %>%
    count(region, name = "n") %>%
    transmute(metric = "region", Status = NA_character_, subclass = region,
              support = NA_character_, n) %>%
    arrange(desc(n))

  ## Per-signature hit counts, so a signature that never fires is visible.
  by_sig <- annotation %>%
    filter(Signature %in% cfg$signatures$signature) %>%
    distinct(seqname, Signature) %>%
    count(Signature, name = "n") %>%
    right_join(cfg$signatures %>% select(Signature = signature, domain_class),
               by = "Signature") %>%
    mutate(n = dplyr::coalesce(n, 0L)) %>%
    transmute(metric = "signature", Status = domain_class, subclass = Signature,
              support = NA_character_, n) %>%
    arrange(desc(n))

  bind_rows(by_class, by_region, by_sig)
}

write_run_info <- function(path, meta) {
  yaml_scalar <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x)) return('""')
    x <- as.character(x)
    ## Quote when the value contains a YAML-significant character. fixed = FALSE
    ## regex was mis-parsing the brace class, so test the characters directly.
    if (any(strsplit(x, "")[[1]] %in% c(":", "#", "[", "]", "{", "}", '"', ","))) {
      paste0('"', gsub('"', '\\\\"', x), '"')
    } else x
  }
  lines <- c("# NLRtracker run provenance")
  for (nm in names(meta)) {
    v <- meta[[nm]]
    if (length(v) > 1) {
      lines <- c(lines, paste0(nm, ":"),
                 paste0("  - ", vapply(v, yaml_scalar, character(1))))
    } else {
      lines <- c(lines, paste0(nm, ": ", yaml_scalar(v)))
    }
  }
  write_lines_out(lines, path)
}

## iTOL domain-architecture strips.
write_itol <- function(domains_tbl, template, path, dedup_map = NULL) {
  colour_of <- c(R1 = "#FBB4AE", CC = "#B3CDE3", RPW8 = "#CCEBC5", TIR = "#DECBE4",
                 PLOOP = "#FFFFCC", `NB-ARC` = "#FED9A6", LRR = "#E5D8BD",
                 `C-JID` = "#FDDAEC", MLKL = "#CCEBC5", ID = "#666666")
  shape_of  <- c(R1 = "HH", CC = "OC", RPW8 = "OC", TIR = "OC", PLOOP = "RE",
                 `NB-ARC` = "RE", LRR = "EL", `C-JID` = "HH", MLKL = "OC", ID = "HH")

  strips <- domains_tbl %>%
    filter(type == "REGION") %>%
    mutate(colour = dplyr::coalesce(unname(colour_of[description]), "#666666"),
           shape  = dplyr::coalesce(unname(shape_of[description]), "HH"),
           piece  = paste0(shape, "|", start, "|", end, "|", colour, "|", description)) %>%
    arrange(seqname, start) %>%
    group_by(seqname) %>%
    summarise(attribute = paste(piece, collapse = ","), .groups = "drop")

  lengths <- domains_tbl %>% filter(type == "CHAIN") %>% select(seqname, end, Simple)

  rows <- lengths %>% left_join(strips, by = "seqname") %>%
    filter(!is.na(attribute))

  if (!is.null(dedup_map)) {
    rows <- rows %>%
      left_join(dedup_map, by = "seqname") %>%
      mutate(group = dplyr::coalesce(group, seqname)) %>%
      group_by(group) %>%
      arrange(desc(nchar(Simple)), desc(end), .by_group = TRUE) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(seqname = group)
  }

  body <- paste0(rows$seqname, ",", rows$end, ",", rows$attribute)
  write_lines_out(c(readr::read_lines(template, progress = FALSE), body), path)
}

write_all <- function(outdir, classification, domains, regions, annotation, aa,
                      cfg, dropped, meta, itol_template) {
  ## Own the output tree here so the R stage also works when driven directly,
  ## not only through NLRtracker.sh.
  for (d in c("figures", file.path("figures", "by_subclass"), "sequences",
              "annotations", file.path("raw", "logs"))) {
    dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
  }

  ## ------------------------------------------------------------- main tables
  classification %>%
    select(seqname, Status, `Subclass (putative)`, subclass_support,
           Domain, `Domain architecture simplified` = Simple) %>%
    write_tsv_out(file.path(outdir, "NLRtracker.tsv"))

  domains_tbl <- build_domains_table(regions, domains, annotation, classification, aa, cfg)
  write_tsv_out(domains_tbl, file.path(outdir, "Domains.tsv"))
  write_tsv_out(build_summary(classification, regions, annotation, cfg),
                file.path(outdir, "summary.tsv"))

  ## Regions trimmed out of existence, reported rather than hidden.
  if (nrow(dropped) > 0) {
    write_tsv_out(dropped, file.path(outdir, "raw", "logs", "dropped_regions.tsv"))
  }

  ## ------------------------------------------------------------- id groups
  nlr <- classification %>% filter(Status %in% c("NLR", "degenerate NLR")) %>% pull(seqname)
  assoc <- classification %>%
    filter(Status %in% c("CCX", "RPW8", "TX", "MLKL", "CJID")) %>% pull(seqname)

  seqdir <- file.path(outdir, "sequences")
  anndir <- file.path(outdir, "annotations")

  write_lines_out(nlr,   file.path(seqdir, "NLR.lst"))
  write_lines_out(assoc, file.path(seqdir, "NLR-associated.lst"))
  write_fasta(aa %>% filter(seqname %in% nlr),   file.path(seqdir, "NLR.fasta"))
  write_fasta(aa %>% filter(seqname %in% assoc), file.path(seqdir, "NLR-associated.fasta"))

  ## ------------------------------------------------------------- region FASTA
  nlr_regions <- regions %>% filter(seqname %in% nlr)
  region_files <- c("NB-ARC" = "NBARC", "TIR" = "TIR", "CC" = "CC", "RPW8" = "RPW8",
                    "R1" = "R1", "LRR" = "LRR", "C-JID" = "CJID", "MLKL" = "MLKL",
                    "ID" = "ID")
  nbarc_fa <- NULL
  for (rg in names(region_files)) {
    fa <- region_fasta(nlr_regions, aa, rg)
    if (nrow(fa) == 0) next
    write_fasta(fa, file.path(seqdir, paste0(region_files[[rg]], ".fasta")))
    if (rg == "NB-ARC") nbarc_fa <- fa
  }
  dedup_map <- NULL
  if (!is.null(nbarc_fa) && nrow(nbarc_fa) > 0) {
    dd <- dedup_fasta(nbarc_fa)
    write_fasta(dd, file.path(seqdir, "NBARC_deduplicated.fasta"))
    dedup_map <- dd %>%
      transmute(group = seqname, member = seqname) %>%
      tidyr::separate_rows(member, sep = "//") %>%
      transmute(seqname = stringr::str_replace(member, "_\\([0-9]+\\)$", ""), group) %>%
      distinct(seqname, .keep_all = TRUE)
  }

  ## ------------------------------------------------------------- annotations
  write_gff3(annotation, regions, nlr,   cfg, file.path(anndir, "NLR.gff3"))
  write_gff3(annotation, regions, assoc, cfg, file.path(anndir, "NLR-associated.gff3"))

  if (file.exists(itol_template)) {
    write_itol(domains_tbl %>% filter(seqname %in% nlr), itol_template,
               file.path(anndir, "iTOL.txt"))
    write_itol(domains_tbl %>% filter(seqname %in% nlr), itol_template,
               file.path(anndir, "iTOL_dedup.txt"), dedup_map = dedup_map)
  }

  write_run_info(file.path(outdir, "run_info.yaml"), meta)

  invisible(list(domains_tbl = domains_tbl, nlr = nlr, assoc = assoc))
}
