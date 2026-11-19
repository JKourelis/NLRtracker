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
  if (!identical(seqtype, "p")) {
    stop("nucleotide input (-t n / --seqtype n) was retired; translate sequences to protein FASTA first",
         call. = FALSE)
  }
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

## Evidence-based split NB-ARC HMM envelope joining.
##
## hmmsearch reports one NB-ARC as two envelopes when an insertion interrupts the
## alignment (for example Rpp4C4). Two consecutive envelopes A then B on one protein
## join when all five conditions hold:
##   1. model carries on:      B$hmm_from >= A$hmm_to - 1
##   2. model does not repeat: B$hmm_from > A$hmm_from
##   3. protein non-overlap:   B$start > A$end
##   4. gap is empty:          no unrelated domain contained in [A$end + 1, B$start - 1]
##   5. corroboration:         a FIMO NB-ARC motif run spans the gap, OR a signature-
##                             derived NBARC/NACHT domain spans the gap, OR the gap
##                             is <= NBARC_MAX_UNCORROBORATED_INSERT (30 aa backstop).
## Chains are bounded so total model span <= NBARC_HMM_LENGTH and total backstop-only
## gap <= NBARC_MAX_UNCORROBORATED_INSERT.
NBARC_MAX_UNCORROBORATED_INSERT <- 30L
NBARC_GAP_IGNORE <- c("NBARC", "NACHT", "PLOOP", "nbarc-motif",
                      "Winged-helix", "nacht-winged-helix")
NBARC_SPANNING_CLASSES <- c("NBARC", "NACHT")

nbarc_join_basis <- function(prev, row, seq_domains, seq_motifs) {
  if (is.null(prev) || !identical(prev$seqname, row$seqname)) return(NA_character_)
  if (is.na(prev$hmm_to) || is.na(row$hmm_from) ||
      is.na(prev$hmm_from) || is.na(row$hmm_to)) return(NA_character_)
  if (row$hmm_from < prev$hmm_to - 1L) return(NA_character_)
  if (row$hmm_from <= prev$hmm_from)  return(NA_character_)
  if (row$start <= prev$end)          return(NA_character_)
  if (max(prev$hmm_to, row$hmm_to) - prev$hmm_from + 1L > NBARC_HMM_LENGTH) {
    return(NA_character_)
  }

  gap_start <- prev$end + 1L
  gap_end   <- row$start - 1L
  if (!is.null(seq_domains) && nrow(seq_domains) > 0L && gap_end >= gap_start) {
    unrelated <- !(seq_domains$domain_class %in% NBARC_GAP_IGNORE)
    if (any(unrelated & seq_domains$start >= gap_start & seq_domains$end <= gap_end)) {
      return(NA_character_)
    }
  }

  if (!is.null(seq_motifs) && nrow(seq_motifs) > 0L) {
    if (any(seq_motifs$start <= prev$end & seq_motifs$end >= row$start)) {
      return("motif")
    }
  }
  if (!is.null(seq_domains) && nrow(seq_domains) > 0L) {
    if (any(seq_domains$domain_class %in% NBARC_SPANNING_CLASSES &
            seq_domains$start <= prev$end & seq_domains$end >= row$start)) {
      return("signature")
    }
  }

  gap <- gap_end - gap_start + 1L
  if (prev$uncorroborated_insert + gap <= NBARC_MAX_UNCORROBORATED_INSERT) {
    return("backstop")
  }
  NA_character_
}

join_split_nbarc <- function(rows, domains = NULL, motif_runs = NULL) {
  if (is.null(rows) || nrow(rows) == 0L) return(rows)
  rows <- rows %>%
    arrange(seqname, start, end) %>%
    mutate(join_basis = NA_character_, joined_from = 1L)
  if (nrow(rows) <= 1L) return(rows)

  counts <- rows %>% count(seqname)
  multi_seqs <- counts$seqname[counts$n > 1L]
  if (length(multi_seqs) == 0L) return(rows)

  singles <- rows %>% filter(!(seqname %in% multi_seqs))
  multis  <- rows %>% filter(seqname %in% multi_seqs)

  dom_by_seq <- if (!is.null(domains) && nrow(domains) > 0L) {
    d_sub <- domains %>% filter(seqname %in% multi_seqs)
    if (nrow(d_sub) > 0L) split(d_sub, d_sub$seqname) else list()
  } else list()

  mot_by_seq <- if (!is.null(motif_runs) && nrow(motif_runs) > 0L) {
    m_sub <- motif_runs %>% filter(seqname %in% multi_seqs)
    if (nrow(m_sub) > 0L) split(m_sub, m_sub$seqname) else list()
  } else list()

  multi_split <- split(multis, multis$seqname)
  joined_list <- lapply(names(multi_split), function(sid) {
    df <- multi_split[[sid]]
    sd <- dom_by_seq[[sid]]
    sm <- mot_by_seq[[sid]]
    acc <- vector("list", nrow(df))
    k <- 0L
    for (i in seq_len(nrow(df))) {
      row <- as.list(df[i, , drop = FALSE])
      row$uncorroborated_insert <- 0L
      if (k == 0L) {
        k <- 1L
        acc[[1L]] <- row
        next
      }
      prev <- acc[[k]]
      basis <- nbarc_join_basis(prev, row, sd, sm)
      if (is.na(basis)) {
        k <- k + 1L
        acc[[k]] <- row
      } else {
        gap <- row$start - prev$end - 1L
        prev$end <- max(prev$end, row$end)
        prev$hmm_to <- max(prev$hmm_to, row$hmm_to)
        prev$score <- prev$score + row$score
        prev$completeness <- round((prev$hmm_to - prev$hmm_from + 1L) / NBARC_HMM_LENGTH, 3)
        if (identical(basis, "backstop")) {
          prev$uncorroborated_insert <- prev$uncorroborated_insert + gap
        }
        seen <- if (is.na(prev$join_basis) || !nzchar(prev$join_basis)) character() else
          strsplit(prev$join_basis, ",", fixed = TRUE)[[1]]
        if (!(basis %in% seen)) seen <- c(seen, basis)
        prev$join_basis <- paste(seen, collapse = ",")
        prev$joined_from <- prev$joined_from + 1L
        acc[[k]] <- prev
      }
    }
    dplyr::bind_rows(lapply(acc[seq_len(k)], function(x) {
      x$uncorroborated_insert <- NULL
      tibble::as_tibble(x)
    }))
  })

  dplyr::bind_rows(singles, dplyr::bind_rows(joined_list)) %>%
    arrange(seqname, start, end)
}
