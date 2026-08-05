## 01_parse.R -- read InterProScan 6 TSV, FIMO GFF and hmmsearch domtblout into one
## long annotation table.
##
## InterProScan 6 TSV, 15 columns, no header. Verified against real 6.0.1 output --
## the published documentation is wrong about columns 10 and 15, see
## docs/interproscan6-notes.md.
##
##  1 seqname   2 md5          3 seq_length  4 analysis    5 signature
##  6 sig_desc  7 start        8 end         9 score      10 source (NOT "T")
## 11 date     12 ipr_entry   13 ipr_desc   14 go_terms   15 pathways ("-" if empty)

IPS6_COLS <- c("seqname", "md5", "seq_length", "analysis", "signature", "sig_desc",
               "start", "end", "score", "source", "date", "ipr_entry", "ipr_desc",
               "go_terms", "pathways")

## InterProScan strips these characters from identifiers, so the FASTA and the
## annotation tables must be normalised the same way.
normalise_id <- function(x) stringr::str_replace_all(x, "/", "")

read_interproscan <- function(path, cfg, seqtype = "p") {
  ## quote = "" for the same reason as in 00_config.R: InterPro descriptions contain bare
  ## double quotes. IPR036390 is '"Winged helix" DNA-binding domain', and with readr's
  ## default quoting the quotes are silently stripped from column 6, so SSF46785 was named
  ## 'Winged helix DNA-binding domain'. A field with an ODD number of quotes would shift
  ## the columns instead of merely losing a character.
  raw <- readr::read_tsv(path, col_names = IPS6_COLS, quote = "", col_types = readr::cols(
    .default = readr::col_character(), start = readr::col_integer(),
    end = readr::col_integer(), seq_length = readr::col_integer()),
    progress = FALSE, na = character())

  if (nrow(raw) == 0) stop("InterProScan TSV is empty: ", path, call. = FALSE)

  ## Guard against being handed an InterProScan 5 file, or a GFF3.
  if (!any(grepl("^(PF|SM|PS|PR|SSF|G3DSA|cd|PTHR|TIGR|NF)", raw$signature))) {
    stop("column 5 of ", path, " does not look like signature accessions.\n",
         "  This parser needs InterProScan 6 TSV (--formats tsv).", call. = FALSE)
  }

  out <- raw %>%
    mutate(
      seqname = normalise_id(seqname),
      ## For --nucleic, InterProScan reports <nucleotide_id>_<orf_id>.
      seqname = if (seqtype == "n") stringr::str_extract(seqname, "^(.*)(?=_orf)") %>%
                                    dplyr::coalesce(seqname) else seqname,
      ipr_entry = dplyr::na_if(ipr_entry, "-"),
      ipr_desc  = dplyr::na_if(ipr_desc, "-"),
      sig_desc  = dplyr::na_if(sig_desc, "-"),
      score     = dplyr::na_if(score, "-"),
      feature   = analysis
    ) %>%
    filter(!is.na(start), !is.na(end))

  ## Name resolution, most to least specific: signature description (column 6),
  ## InterPro entry name (column 13), a locally supplied fallback for the nine
  ## accessions that have neither, then the bare accession.
  fb <- NAME_FALLBACK[out$signature]
  out <- out %>%
    mutate(Name = dplyr::coalesce(sig_desc, ipr_desc, unname(fb), signature),
           Signature = signature,
           ENTRY_AC = dplyr::coalesce(ipr_entry, "Not integrated"))

  ## ENTRY_TYPE comes from the InterPro entry list; the TSV does not carry it.
  out <- out %>%
    left_join(cfg$entries %>% select(ENTRY_AC, ENTRY_TYPE), by = "ENTRY_AC")

  out %>%
    mutate(strand = "+", frame = ".",
           attribute = paste0("Name=", Name, ";Signature=", Signature,
                              ";InterPro=", ENTRY_AC, ";analysis=", analysis),
           attribute = stringr::str_replace_all(attribute, "[\r\n\t]+", " ")) %>%
    select(seqname, source = feature, feature, start, end, score, strand, frame,
           attribute, Name, Signature, ENTRY_AC, ENTRY_TYPE, seq_length)
}

read_fimo <- function(path) {
  raw <- readr::read_tsv(path, col_names = FALSE, comment = "#", quote = "",
                         col_types = readr::cols(.default = readr::col_character()),
                         progress = FALSE, na = character())
  if (nrow(raw) == 0) stop("FIMO output is empty: ", path, call. = FALSE)
  names(raw)[1:9] <- c("seqname", "source", "feature", "start", "end",
                       "score", "strand", "frame", "attribute")

  raw %>%
    mutate(
      motif  = stringr::str_match(attribute, "Name=([^;_]+)")[, 2],
      qvalue = suppressWarnings(as.numeric(
        stringr::str_match(attribute, "qvalue=([^;]+)")[, 2])),
      pvalue = stringr::str_match(attribute, "pvalue=([^;]+)")[, 2],
      sequence = stringr::str_match(attribute, "sequence=([^;]+)")[, 2],
      start = as.integer(start), end = as.integer(end),
      score_num = suppressWarnings(as.numeric(score))
    ) %>%
    filter(!is.na(start), !is.na(score_num),
           score_num >= FIMO_MIN_SCORE,
           is.na(qvalue) | qvalue <= FIMO_MAX_QVALUE) %>%
    mutate(
      seqname   = normalise_id(seqname),
      Signature = paste("Motif", motif),
      Name      = dplyr::coalesce(unname(MEME_MOTIF_NAMES[motif]), "OTHER"),
      source    = "fimo",
      feature   = "fimo",
      strand    = "+",
      frame     = ".",
      score     = as.character(score_num),
      ENTRY_AC  = "Not integrated",
      ENTRY_TYPE = NA_character_,
      seq_length = NA_integer_,
      attribute = paste0("Name=", Name, ";Signature=", Signature,
                         ";pvalue=", pvalue, ";qvalue=", qvalue),
      attribute = stringr::str_replace_all(attribute, "[\r\n\t]+", " ")
    ) %>%
    select(seqname, source, feature, start, end, score, strand, frame,
           attribute, Name, Signature, ENTRY_AC, ENTRY_TYPE, seq_length)
}

## hmmsearch --domtblout. Whitespace-delimited, three comment lines, 23 fields.
## Columns used: 1 target, 14 domain score, 20/21 envelope start/end.
read_hmmer <- function(path, min_score = 10) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  lines <- readr::read_lines(path, progress = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  if (!length(lines)) return(NULL)

  f <- stringr::str_split(lines, "\\s+", n = 23, simplify = TRUE)
  if (ncol(f) < 23) return(NULL)

  tibble::tibble(
    seqname = normalise_id(f[, 1]),
    score_domain = suppressWarnings(as.numeric(f[, 14])),
    start = suppressWarnings(as.integer(f[, 20])),
    end   = suppressWarnings(as.integer(f[, 21]))
  ) %>%
    filter(!is.na(start), !is.na(end), !is.na(score_domain),
           score_domain >= min_score) %>%
    mutate(
      source = "HMMER3", feature = "HMMER3", strand = "+", frame = ".",
      score = as.character(score_domain),
      Name = "C-JID", Signature = "CJID",
      ENTRY_AC = "Not integrated", ENTRY_TYPE = NA_character_,
      seq_length = NA_integer_,
      attribute = paste0("Name=C-JID;region_name=C-terminal jelly roll/Ig-like domain",
                         ";note=jelly roll/Ig-like domain [structural motif]",
                         ";NCBI Feature Key=Region")
    ) %>%
    select(seqname, source, feature, start, end, score, strand, frame,
           attribute, Name, Signature, ENTRY_AC, ENTRY_TYPE, seq_length)
}

read_sequences <- function(path) {
  aa <- Biostrings::readAAStringSet(path)
  tibble::tibble(seqname = normalise_id(names(aa)),
                 sequence = as.character(aa)) %>%
    mutate(seqname = stringr::str_extract(seqname, "^\\S+")) %>%
    distinct(seqname, .keep_all = TRUE)
}

## Combine everything into one annotation table, and attach the class each
## signature belongs to.
build_annotation <- function(ips, fimo, hmmer, cfg) {
  ann <- dplyr::bind_rows(ips, fimo, hmmer)
  ann %>%
    mutate(domain_class = unname(cfg$class_of[Signature])) %>%
    arrange(seqname, start, end)
}

## ------------------------------------------------------------------ NB-ARC HMM

## Read `hmmsearch --domtblout` from the locked NB-ARC_region.hmm.
##
## Thresholds are calibrated, not guessed (docs/nbarc-hmm.md §1):
##   domain_ga  20  -- where are the NB-ARCs inside a protein? Recovers RRS1's second,
##                     partial NB-ARC at score 36, which signatures miss entirely.
##   protein_ga 40  -- is this protein an NLR at all? Every one of the 481 RefPlantNLR
##                     proteins scores >= 40.7; AAA+ ATPases and CDC48/p97 top out at 15.8.
## One threshold cannot do both jobs, which is why there are two.
NBARC_HMM_DOMAIN_GA  <- 20
NBARC_HMM_PROTEIN_GA <- 40
NBARC_HMM_LENGTH     <- 382

read_nbarc_hmm <- function(path, domain_ga = NBARC_HMM_DOMAIN_GA) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  lines <- readr::read_lines(path, progress = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  if (!length(lines)) return(NULL)
  f <- stringr::str_split(lines, "\\s+", n = 23, simplify = TRUE)
  if (ncol(f) < 23) return(NULL)

  tibble::tibble(
    ## normalise exactly as the rest of the pipeline does. RefPlantNLR contains
    ## "Bph1/9-1", and stripping "/" here is what makes it join to the annotation table.
    seqname  = normalise_id(f[, 1]),
    score    = suppressWarnings(as.numeric(f[, 14])),
    hmm_from = suppressWarnings(as.integer(f[, 16])),
    hmm_to   = suppressWarnings(as.integer(f[, 17])),
    start    = suppressWarnings(as.integer(f[, 20])),   # envelope
    end      = suppressWarnings(as.integer(f[, 21]))
  ) %>%
    filter(!is.na(start), !is.na(end), !is.na(score), score >= domain_ga) %>%
    mutate(completeness = round((hmm_to - hmm_from + 1) / NBARC_HMM_LENGTH, 3)) %>%
    arrange(seqname, start, end)
}
