#!/usr/bin/env Rscript
## NLRtracker.R -- driver.
##
## Called by NLRtracker.sh with named arguments. All biological logic lives in
## module/nlr_signatures.tsv, module/nlr_classes.tsv and module/nlr_subclass_rules.tsv;
## all code lives in module/R/.

## Component packages only. The tidyverse meta-package is deliberately not required:
## it adds an install dependency without adding capability.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

parse_args <- function(argv) {
  out <- list(); i <- 1
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (!startsWith(a, "--")) stop("unexpected argument: ", a, call. = FALSE)
    key <- sub("^--", "", a)
    if (i + 1 > length(argv) || startsWith(argv[[i + 1]], "--")) {
      out[[key]] <- TRUE; i <- i + 1
    } else {
      out[[key]] <- argv[[i + 1]]; i <- i + 2
    }
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
need_arg <- function(k) {
  if (is.null(args[[k]]) || !nzchar(as.character(args[[k]]))) {
    stop("missing required argument --", k, call. = FALSE)
  }
  as.character(args[[k]])
}
opt <- function(k, default = NA_character_) {
  v <- args[[k]]
  if (is.null(v) || !nzchar(as.character(v))) default else as.character(v)
}

here <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in sort(list.files(file.path(here, "R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

outdir <- need_arg("outdir")
message("[R] loading configuration")
cfg <- load_config(
  signatures     = need_arg("signatures"),
  classes        = need_arg("classes"),
  subclass_rules = need_arg("subclass-rules"),
  entry_list     = need_arg("entry-list")
)
message("[R]   ", nrow(cfg$signatures), " signatures across ",
        nrow(cfg$classes), " classes; ", nrow(cfg$rules), " multi-signature rule(s)")

seqtype <- opt("seqtype", "p")
if (!identical(seqtype, "p")) {
  stop("nucleotide input (--seqtype n) was retired; translate sequences to protein FASTA first",
       call. = FALSE)
}
message("[R] parsing InterProScan 6 TSV")
ips <- read_interproscan(need_arg("interproscan"), cfg, seqtype = seqtype)
message("[R]   ", nrow(ips), " matches over ", dplyr::n_distinct(ips$seqname), " proteins")

message("[R] parsing FIMO")
fimo <- read_fimo(need_arg("fimo"))
message("[R]   ", nrow(fimo), " motif hits above threshold")

hmmer <- read_hmmer(opt("hmmer"))
message("[R] parsing HMMER: ",
        if (is.null(hmmer)) "skipped" else paste(nrow(hmmer), "C-JID hits"))

nbarc_hmm <- read_nbarc_hmm(opt("nbarc-hmm"))
message("[R] parsing NB-ARC HMM: ",
        if (is.null(nbarc_hmm)) "skipped (no --nbarc-hmm)" else
          paste(nrow(nbarc_hmm), "domains over",
                dplyr::n_distinct(nbarc_hmm$seqname), "proteins"))

aa <- read_sequences(need_arg("fasta"))
message("[R]   ", nrow(aa), " sequences read")

annotation <- build_annotation(ips, fimo, hmmer, cfg)
unclassified <- annotation %>% filter(is.na(domain_class)) %>% distinct(Signature)
message("[R] ", nrow(annotation), " annotations; ", nrow(unclassified),
        " signature(s) not in the NLR signature table (ignored)")

message("[R] building domains")
domains <- build_domains(annotation, cfg)
message("[R]   ", nrow(domains), " domain intervals")

## Regions are built from the full domain set, because scaffold classes such as
## Winged-helix exist precisely to extend NB-ARC. They are dropped afterwards.
message("[R] arbitrating regions")
reg <- build_regions(domains, cfg, nbarc_hmm = nbarc_hmm, annotation = annotation,
                     require_signature_support =
                       !identical(opt("nbarc-hmm-unsupported"), "true"))
if (!is.null(reg$hmm)) {
  message("[R]   NB-ARC from HMM: ", reg$hmm$replaced, " boundaries replaced, ",
          reg$hmm$added, " protein(s) added",
          if (!is.null(reg$hmm$rejected)) paste0(", ", nrow(reg$hmm$rejected),
                                                 " HMM-only hit(s) rejected") else "")
}
regions <- reg$regions
overlaps <- check_regions_disjoint(regions)
if (nrow(overlaps) > 0) {
  message("[R] WARNING: ", nrow(overlaps), " residual region overlap(s) after arbitration")
  print(utils::head(as.data.frame(overlaps), 10))
} else {
  message("[R]   ", nrow(regions), " regions, all disjoint")
}
if (nrow(reg$dropped) > 0) {
  message("[R]   ", nrow(reg$dropped), " region(s) trimmed out of existence:")
  reg$dropped %>% count(rule) %>% as.data.frame() %>% print(row.names = FALSE)
}

message("[R] classifying")
domains_kept <- drop_scaffold_classes(domains, cfg)
classification <- classify(domains, annotation, cfg)
message("[R]   ", nrow(classification), " sequences classified")
classification %>% count(Status, `Subclass (putative)`) %>%
  arrange(desc(n)) %>% as.data.frame() %>% utils::head(15) %>% print(row.names = FALSE)

message("[R] writing output")
meta <- list(
  nlrtracker_version      = opt("version", "unknown"),
  interproscan_version    = opt("interproscan-version", "unknown"),
  interpro_version        = opt("interpro-version", "unknown"),
  applications            = opt("applications", "unknown"),
  interproscan_batch_size = opt("batch-size", "unknown"),
  seqtype                 = seqtype,
  input                   = opt("input-name", "unknown"),
  n_input_sequences       = nrow(aa),
  n_classified            = nrow(classification),
  n_regions               = nrow(regions),
  signature_table_rows    = nrow(cfg$signatures),
  R_version               = paste(R.version$major, R.version$minor, sep = "."),
  run_date                = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)

res <- write_all(outdir, classification, domains_kept, regions, annotation, aa,
                 cfg, reg$dropped, meta,
                 itol_template = file.path(here, "iTOL_NLR_template.txt"),
                 reference_itol = file.path(here, "reference", "NB-ARC_region.reference.itol"))

message("[R] drawing figures")
draw_all_figures(res$domains_tbl, outdir, res$nlr, res$assoc)

phylo <- NULL
if (isTRUE(args[["phylogeny"]]) || identical(opt("phylogeny"), "true")) {
  message("[R] phylogenetic classification")
  phylo <- run_phylogeny(outdir, classification, regions, aa, cfg,
                         clades_path   = file.path(here, "nlr_clades.tsv"),
                         reference_dir = file.path(here, "reference"),
                         threads       = as.integer(dplyr::coalesce(opt("threads"), "1")))
}

message("[R] writing report")
write_report(outdir, classification, res$domains_tbl, regions, annotation, cfg, meta,
             template = file.path(here, "report_template.html"), phylo = phylo)

message("[R] done")
