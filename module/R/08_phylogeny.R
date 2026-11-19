## 08_phylogeny.R -- classify NLRs by NB-ARC phylogeny against a RefPlantNLR scaffold.
##
## Domain architecture is a proxy for phylogeny, and it fails wherever an N-terminal domain is
## undetectable. On real proteomes that is 24-33% of NLRs, against 5% on the curated RefPlantNLR
## set. Clade membership does not depend on a domain being detectable, so it answers those cases.
##
## Alignment uses module/reference/NB-ARC_region.hmm -- built from RefPlantNLR NB-ARC regions and
## shipped frozen -- not Pfam PF00931. Two measured reasons:
##   * PF00931 is LENG 175 against our region median of 356 aa, so it covers about half of
##     NBD-HD1-WHD and discards exactly the part the region definition exists to capture.
##   * PF00931 sits on a detection cliff: the TNP outgroup scores 26.0-26.1 against its gathering
##     threshold of 23.5. Our model scores the same sequences 72-89.
##
## A frozen profile also means the alignment does not shift when the query set changes, which is
## what makes clade assignment comparable across proteomes.

PHYLO_MIN_OCCUPANCY <- 0.50   # fraction of match-state columns a tip must cover to enter the tree
PHYLO_MAX_FOREIGN   <- 0.20   # report (not reject) a clade subtree containing this many foreign anchors
PHYLO_MAX_SPAN      <- 0.90   # reject a clade whose MRCA spans this share of the reference set
PHYLO_ANCHOR_QUORUM <- 1.00   # anchors bind: smallest node holding ALL declared anchors
PHYLO_MAX_COPHEN    <- 8000   # tips above which the full patristic matrix is not built

## ------------------------------------------------------------------ helpers

read_fasta_simple <- function(path) {
  lines <- readr::read_lines(path, progress = FALSE)
  h <- grepl("^>", lines)
  ids <- sub("^>", "", lines[h])
  ids <- sub("\\s.*$", "", ids)
  grp <- cumsum(h)
  seqs <- vapply(split(lines[!h], grp[!h]), paste0, character(1), collapse = "")
  tibble::tibble(id = ids, seq = unname(seqs[as.character(seq_along(ids))]))
}

write_fasta_simple <- function(df, path) {
  readr::write_lines(as.vector(rbind(paste0(">", df$id), df$seq)), path)
}

## ------------------------------------------------------------------ alignment

## hmmalign with --outformat A2M marks insert columns in lowercase, so keeping only the
## upper-case and gap positions leaves exactly the match states. That is why A2M is used
## rather than Stockholm: esl-alimask is not installed, and this needs no parser.
## `protect` lists ids that must never be dropped. The occupancy filter exists to keep
## unplaceable QUERY fragments out of the tree, but the reference scaffold is the coordinate
## system: dropping a reference anchor silently breaks the clade it defines. This bit hard --
## the TNP outgroup aligns over only 46-48% of the match states, because an outgroup is by
## construction the most divergent thing present, so the filter removed the very sequences
## used to root the tree.
align_to_profile <- function(fasta, hmm, out_aln, protect = character()) {
  a2m <- paste0(out_aln, ".a2m")
  st <- system2("hmmalign",
                c("--amino", "--trim", "--outformat", "A2M", shQuote(hmm), shQuote(fasta)),
                stdout = a2m, stderr = FALSE)
  if (!identical(st, 0L) || !file.exists(a2m)) stop("hmmalign failed", call. = FALSE)

  fa <- read_fasta_simple(a2m)
  ## In A2M the insert columns are NOT aligned across sequences -- each row carries its own
  ## inserts, so rows have different lengths and a shared column mask is meaningless. (Using
  ## one made VeryFastTree reject the alignment: "expected 456 but have 382".) Deleting the
  ## lower-case and '.' positions per sequence leaves exactly the model's match states, so
  ## every row ends up the same width by construction.
  fa$seq <- toupper(gsub("[a-z.]", "", fa$seq))

  widths <- unique(nchar(fa$seq))
  if (length(widths) != 1) {
    stop("match-state extraction produced ragged rows: widths ",
         paste(utils::head(sort(widths), 5), collapse = ", "), call. = FALSE)
  }
  ncol_match <- widths
  occ <- vapply(strsplit(fa$seq, "", fixed = TRUE),
                function(v) mean(v != "-"), numeric(1))
  fa$occupancy <- occ
  keep_row <- occ >= PHYLO_MIN_OCCUPANCY | fa$id %in% protect
  dropped <- fa[!keep_row, c("id", "occupancy")]
  kept <- fa[keep_row, ]

  write_fasta_simple(kept[, c("id", "seq")], out_aln)
  unlink(a2m)
  list(aln = out_aln, n_cols = ncol_match, kept = kept, dropped = dropped)
}

## ------------------------------------------------------------------ safe tip names

## Newick reserves ( ) : , ; ' [ ] and whitespace, and tools sanitise them differently.
## RefPlantNLR alone contains "Prv^1", "L^1a" and the "_(n)" multi-copy suffix, which
## VeryFastTree collapsed into a duplicate name ("Non-unique name 'Prv^1_'"). Rather than
## guess which characters each tool tolerates, sequences go into the tree under opaque
## tokens and are mapped back afterwards. Arbitrary proteome identifiers then cannot break it.
make_tip_map <- function(ids) {
  tibble::tibble(id = ids, tip = sprintf("t%06d", seq_along(ids)))
}

## Quote a label for Newick when it contains a reserved character, doubling any single quote.
## A label that needs no quoting is left bare, so ordinary identifiers stay readable.
newick_label <- function(x) {
  needs <- grepl("[][(){},:;'[:space:]]", x)
  ifelse(needs, paste0("'", gsub("'", "''", x, fixed = TRUE), "'"), x)
}

## ------------------------------------------------------------------ tree

## `threads` defaults to 1, deliberately.
##
## VeryFastTree 4.0.5 is not reproducible above one thread, and on this data it also returns a
## much worse tree. Measured on the 507-tip RefPlantNLR scaffold, each thread count run twice:
##
##   threads 1     6 of 414 splits reported bad    reproducible       2.71 s
##   threads 2   401 of 414 splits reported bad    not reproducible   1.96 s
##   threads 4   378 of 414 splits reported bad    not reproducible   1.42 s
##   threads 8   378 of 414 splits reported bad    not reproducible   1.47 s
##
## The whole saving is about 1.2 seconds, against a tree that changes between runs and is
## poorly resolved. Clade assignment survived the difference, because the quorum rule tolerates
## local rearrangement, but a tree that cannot be reproduced is not one to publish a figure from.
##
## `-double-precision` is passed to VeryFastTree because 32-bit single precision can hang in
## `PairLogLk -- numerical problem!` on large alignments while remaining reproducible at
## `-threads 1`. Stderr is captured to a temporary log so failures surface diagnostics.
build_tree <- function(aln, out_tree, threads = 1) {
  prog <- if (nzchar(Sys.which("VeryFastTree"))) "VeryFastTree" else
          if (nzchar(Sys.which("FastTree"))) "FastTree" else ""
  if (!nzchar(prog)) stop("neither VeryFastTree nor FastTree found", call. = FALSE)
  args <- if (prog == "VeryFastTree") {
    c("-double-precision", "-lg", "-gamma", "-pseudo", "-quiet", "-threads", threads, shQuote(aln))
  } else {
    c("-lg", "-gamma", "-pseudo", "-quiet", shQuote(aln))
  }
  err_log <- tempfile("tree_build_", fileext = ".stderr.log")
  on.exit(unlink(err_log), add = TRUE)
  st <- system2(prog, args, stdout = out_tree, stderr = err_log)
  if (!identical(st, 0L) || !file.exists(out_tree) || file.size(out_tree) == 0) {
    err_txt <- if (file.exists(err_log)) readLines(err_log, warn = FALSE) else character()
    err_txt <- err_txt[nzchar(trimws(err_txt))]
    detail <- if (length(err_txt)) {
      tail_str <- paste(utils::tail(err_txt, 10L), collapse = " | ")
      if (nchar(tail_str) > 400L) substr(tail_str, nchar(tail_str) - 399L, nchar(tail_str)) else tail_str
    } else ""
    msg <- paste0(prog, " failed (exit ", st, ")")
    if (nzchar(detail)) msg <- paste0(msg, ": ", detail)
    stop(msg, call. = FALSE)
  }
  out_tree
}

## ------------------------------------------------------------------ clade assignment

## Tip labels are NB-ARC region ids, which carry the "_(n)" suffix when a protein has more
## than one. Strip it to recover the protein.
tip_to_protein <- function(x) sub("_\\((\\d+)\\)$", "", x)

## Match curated clade-table anchors (protein names, possibly containing "/") to tree tip
## labels (normalised identifiers, possibly carrying "_(n)" tandem-copy suffixes).
anchor_tips <- function(anchor_list, tip_labels) {
  if (length(anchor_list) == 0L || length(tip_labels) == 0L) return(character())
  tip_proteins <- tip_to_protein(tip_labels)
  want <- unique(normalise_id(trimws(anchor_list)))
  want <- want[nzchar(want)]
  tip_labels[tip_proteins %in% want]
}

load_clades <- function(path) {
  read_tsv_config(path, "clade table") %>%
    mutate(anchor_list = strsplit(dplyr::coalesce(anchors, ""), ",", fixed = TRUE))
}

## For every anchored clade: the smallest subtree holding all its anchors, the tips
## beneath it, and how many anchors of other clades it contains.
clade_subtrees <- function(tree, clades, ref_ids) {
  ntip <- length(tree$tip.label)
  ## descendant tip sets for every internal node, computed once
  parts <- ape::prop.part(tree)
  node_tips <- lapply(parts, function(idx) tree$tip.label[idx])
  node_ids <- ntip + seq_along(node_tips)

  anchor_of <- list()
  for (i in seq_len(nrow(clades))) {
    anchor_of[[clades$clade[i]]] <- anchor_tips(clades$anchor_list[[i]], tree$tip.label)
  }
  all_anchors <- unlist(anchor_of)

  out <- list()
  for (i in seq_len(nrow(clades))) {
    nm <- clades$clade[i]
    a <- anchor_of[[nm]]
    if (length(a) < 2) next
    need <- max(2, ceiling(PHYLO_ANCHOR_QUORUM * length(a)))
    hits <- vapply(node_tips, function(t) sum(t %in% a), integer(1))
    sizes <- vapply(node_tips, length, integer(1))
    ok <- which(hits >= need)
    if (!length(ok)) next
    best <- ok[which.min(sizes[ok])]
    tips <- node_tips[[best]]
    own <- sum(tips %in% a)
    foreign <- sum(tips %in% setdiff(all_anchors, a))
    n_ref <- sum(tips %in% ref_ids)
    out[[nm]] <- list(node = node_ids[best], tips = tips, n_tips = length(tips),
                      n_anchor = length(a), recovered = own / length(a),
                      foreign_frac = if (n_ref > 0) foreign / n_ref else 0)
  }
  out
}

## CC-01 and CC-02 are derived from the tree rather than hardcoded, but NOT as the two
## immediate daughters of the CC node. That was the first attempt and it failed on real data:
## NRC and MLA1 both landed in the SAME daughter, because the deepest split inside CC is not
## the NRC/MLA split. Forcing it would have invented a taxonomy the tree does not support.
##
## Instead each is the smallest subtree holding its own defining anchors -- CC-01 the NRC
## anchors (helpers and sensors together), CC-02 the MLA1 anchors -- and they are required to
## be disjoint. This states what the tree actually shows, and fails loudly when it shows
## something else.
derive_cc_split <- function(tree, sub, clades) {
  anchors_for <- function(nm) {
    i <- match(nm, clades$clade)
    if (is.na(i)) character() else anchor_tips(clades$anchor_list[[i]], tree$tip.label)
  }
  nrc <- union(anchors_for("CC-01-NRCh"), anchors_for("CC-01-NRC-sensor"))
  mla <- anchors_for("CC-02-MLA")
  if (length(nrc) < 2 || length(mla) < 2) {
    return(list(ok = FALSE, why = "too few NRC or MLA1 anchors in the tree"))
  }

  ntip <- length(tree$tip.label)
  parts <- ape::prop.part(tree)
  node_tips <- lapply(parts, function(idx) tree$tip.label[idx])
  smallest_holding <- function(a) {
    need <- max(2, ceiling(PHYLO_ANCHOR_QUORUM * length(a)))
    hits <- vapply(node_tips, function(t) sum(t %in% a), integer(1))
    ok <- which(hits >= need)
    if (!length(ok)) return(NULL)
    node_tips[[ok[which.min(vapply(node_tips[ok], length, integer(1)))]]]
  }
  t1 <- smallest_holding(nrc); t2 <- smallest_holding(mla)
  if (is.null(t1) || is.null(t2)) {
    return(list(ok = FALSE, why = "no node holds a quorum of the NRC or MLA1 anchors"))
  }
  shared <- length(intersect(t1, t2))
  if (shared > 0) {
    return(list(ok = FALSE,
                why = sprintf("the NRC and MLA1 subtrees overlap by %d tips (%s), so they are nested rather than sister",
                              shared,
                              if (all(t2 %in% t1)) "MLA1 inside NRC" else
                              if (all(t1 %in% t2)) "NRC inside MLA1" else "partial")))
  }
  list(ok = TRUE, cc01 = t1, cc02 = t2,
       note = sprintf("CC-01 %d tips from %d NRC anchors, CC-02 %d tips from %d MLA1 anchors, disjoint",
                      length(t1), length(nrc), length(t2), length(mla)))
}

## ------------------------------------------------------------------ nearest reference

nearest_reference <- function(tree, query_tips, ref_tips, k = 4) {
  if (!length(query_tips) || !length(ref_tips)) {
    return(tibble::tibble(tip = character(), nearest_ref = character(),
                          nearest_dist = numeric(), next_refs = character()))
  }
  if (length(tree$tip.label) <= PHYLO_MAX_COPHEN) {
    d <- ape::cophenetic.phylo(tree)
    d <- d[query_tips, ref_tips, drop = FALSE]
    purrr::map_dfr(seq_len(nrow(d)), function(i) {
      o <- order(d[i, ])
      nxt <- if (length(o) >= 2L && k >= 2L) {
        idx <- o[2:min(k, length(o))]
        paste(sprintf("%s(%.3f)", colnames(d)[idx], d[i, idx]), collapse = ", ")
      } else ""
      tibble::tibble(
        tip = rownames(d)[i],
        nearest_ref = colnames(d)[o[1]],
        nearest_dist = round(unname(d[i, o[1]]), 4),
        next_refs = nxt)
    })
  } else {
    ## Too many tips for a full matrix. Walk up from each query until a reference appears in
    ## the subtree, then measure only within that subtree.
    message("[R]   ", length(tree$tip.label), " tips: using subtree search instead of a ",
            "full patristic matrix")
    purrr::map_dfr(query_tips, function(q) {
      node <- which(tree$tip.label == q)
      repeat {
        parent <- tree$edge[tree$edge[, 2] == node, 1]
        if (!length(parent)) return(tibble::tibble(
          tip = q, nearest_ref = NA_character_, nearest_dist = NA_real_, next_refs = ""))
        tips <- ape::extract.clade(tree, parent)$tip.label
        hit <- intersect(tips, ref_tips)
        if (length(hit)) {
          sub <- ape::keep.tip(tree, unique(c(q, hit)))
          d <- ape::cophenetic.phylo(sub)[q, hit, drop = TRUE]
          o <- order(d)
          nxt <- if (length(o) >= 2L && k >= 2L) {
            idx <- o[2:min(k, length(o))]
            paste(sprintf("%s(%.3f)", names(d)[idx], d[idx]), collapse = ", ")
          } else ""
          return(tibble::tibble(
            tip = q, nearest_ref = names(d)[o[1]], nearest_dist = round(unname(d[o[1]]), 4),
            next_refs = nxt))
        }
        node <- parent
      }
    })
  }
}

## ------------------------------------------------------------------ categorise

## Map an architecture subclass onto the clade vocabulary so the two can be compared.
ARCH_TO_CLADE <- c(
  "CC-NLR" = "CC-NLR", "TIR-NLR" = "TIR-NLR", "CCR-NLR" = "CCR-NLR",
  "CCG10-NLR" = "CCG10-NLR", "CCG10-NLR (probable)" = "CCG10-NLR"
)

categorise <- function(arch_subclass, clade, clade_path, n_nbarc, n_discordant, basis = "") {
  amb <- is.na(arch_subclass) | arch_subclass == "UNDETERMINED" |
    grepl(" or ", arch_subclass, fixed = TRUE)
  mapped <- unname(ARCH_TO_CLADE[sub(" \\(probable\\)$", "", arch_subclass)])
  dplyr::case_when(
    basis == "outgroup"                       ~ "outgroup",
    is.na(clade) | clade == "unplaced"        ~ "unplaced",
    n_discordant > 0                          ~ "multi_nbarc_discordant",
    amb                                       ~ "phylogeny_resolves",
    !is.na(mapped) & mapped == clade          ~ "agree",
    !is.na(mapped) & purrr::map2_lgl(mapped, clade_path, ~ !is.na(.x) &&
      grepl(paste0("(^|>)", .x, "(>|$)"), .y)) ~ "phylogeny_refines",
    TRUE                                      ~ "conflict"
  )
}

## ------------------------------------------------------------------ entry point

run_phylogeny <- function(outdir, classification, regions, aa, cfg, clades_path,
                          reference_dir, threads = 1) {
  if (!nzchar(Sys.which("hmmalign"))) {
    message("[R]   hmmalign not found - skipping phylogeny"); return(NULL)
  }
  if (!nzchar(Sys.which("VeryFastTree")) && !nzchar(Sys.which("FastTree"))) {
    message("[R]   neither VeryFastTree nor FastTree found - skipping phylogeny"); return(NULL)
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    message("[R]   ape not installed - skipping phylogeny"); return(NULL)
  }
  hmm <- file.path(reference_dir, "NB-ARC_region.hmm")
  ref_fa <- file.path(reference_dir, "NB-ARC_region.reference.fasta")
  for (f in c(hmm, ref_fa)) {
    if (!file.exists(f)) {
      message("[R]   missing reference file ", f, " - skipping phylogeny"); return(NULL)
    }
  }
  pdir <- file.path(outdir, "phylogeny")
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)

  ## query NB-ARC regions, named exactly as the per-region FASTA names them
  nlr <- classification %>% filter(Status %in% c("NLR", "degenerate NLR")) %>% pull(seqname)
  q_raw <- region_fasta(regions %>% filter(seqname %in% nlr), aa, "NB-ARC")
  if (nrow(q_raw) < 2) { message("[R]   fewer than 2 query NB-ARC regions - skipping"); return(NULL) }
  ref <- read_fasta_simple(ref_fa)
  ## Keep non-colliding sequences for the FASTA, while retaining colliding query proteins
  ## via their scaffold tips after tree construction.
  q <- q_raw %>% rename(id = seqname, seq = sequence) %>% filter(!id %in% ref$id)
  combined <- file.path(pdir, "NB-ARC.input.fasta")
  write_fasta_simple(bind_rows(ref, q), combined)
  message("[R]   ", nrow(ref), " reference + ", nrow(q), " non-reference query NB-ARC regions")

  message("[R]   aligning to the NB-ARC profile")
  al <- align_to_profile(combined, hmm, file.path(pdir, "NB-ARC.aln"), protect = ref$id)
  message("[R]     ", al$n_cols, " match-state columns; ", nrow(al$kept), " tips kept, ",
          nrow(al$dropped), " below ", round(100 * PHYLO_MIN_OCCUPANCY), "% occupancy")

  ## Rename to opaque tokens for the tree, then map straight back. See make_tip_map().
  tipmap <- make_tip_map(al$kept$id)
  safe_aln <- file.path(pdir, "NB-ARC.safe.aln")
  write_fasta_simple(tibble::tibble(id = tipmap$tip, seq = al$kept$seq), safe_aln)

  message("[R]   building tree")
  tp <- build_tree(safe_aln, file.path(pdir, "NB-ARC.tree"), threads)
  raw_newick <- paste(readr::read_lines(tp, progress = FALSE), collapse = "")
  tree_tok <- ape::read.tree(tp)
  back <- setNames(tipmap$id, tipmap$tip)
  real_tips <- unname(back[tree_tok$tip.label])
  if (anyNA(real_tips)) stop("tip name mapping lost a label", call. = FALSE)
  ## Publish the tree by substituting real names into the tree builder's own output, rather
  ## than re-writing a parsed tree through ape::write.tree(), which corrupts "Bph19-1_(1)".
  for (i in seq_len(nrow(tipmap))) {
    raw_newick <- gsub(tipmap$tip[i], newick_label(tipmap$id[i]), raw_newick, fixed = TRUE)
  }
  readr::write_lines(raw_newick, tp)
  unlink(safe_aln)
  message("[R]     ", length(real_tips), " tips")

  cl <- load_clades(clades_path)
  ref_ids <- intersect(ref$id, real_tips)

  ## Root on the outgroup while tips are still opaque tokens, then write NB-ARC.rooted.tree
  ## by substituting real names back into the rooted Newick string.
  og_rows <- cl %>% filter(role == "outgroup")
  outgroup_clade <- character()
  for (i in seq_len(nrow(og_rows))) {
    tips_i <- anchor_tips(og_rows$anchor_list[[i]], real_tips)
    if (length(tips_i)) outgroup_clade[tips_i] <- og_rows$clade[i]
  }
  og <- names(outgroup_clade)
  og_tok <- tipmap$tip[match(og, tipmap$id)]
  og_tok <- og_tok[!is.na(og_tok)]
  rooted <- FALSE
  if (length(og_tok) >= 1) {
    tree_tok <- tryCatch({
      t2 <- ape::root(tree_tok, outgroup = og_tok, resolve.root = TRUE, edgelabel = TRUE)
      t2 <- ape::ladderize(t2)
      rooted <- TRUE
      t2
    }, error = function(e) {
      message("[R]   rooting failed: ", conditionMessage(e))
      tree_tok
    })
  }
  message("[R]   rooted on ", length(og), " outgroup tip(s): ", rooted)
  if (rooted) {
    if (!is.null(tree_tok$node.label) && length(tree_tok$node.label) != tree_tok$Nnode) {
      tree_tok$node.label <- NULL
    }
    rooted_nwk <- ape::write.tree(tree_tok)
    for (i in seq_len(nrow(tipmap))) {
      rooted_nwk <- gsub(tipmap$tip[i], newick_label(tipmap$id[i]), rooted_nwk, fixed = TRUE)
    }
    rooted_tp <- file.path(pdir, "NB-ARC.rooted.tree")
    readr::write_lines(rooted_nwk, rooted_tp)
    message("[R]     rooted tree written to ", basename(rooted_tp))
  }

  tree <- tree_tok
  tree$tip.label <- unname(back[tree$tip.label])
  if (anyNA(tree$tip.label)) stop("rooted tip name mapping lost a label", call. = FALSE)

  sub <- clade_subtrees(tree, cl %>% filter(role == "anchored"), ref_ids)

  n_ref_total <- max(length(ref_ids), 1)
  bad <- names(sub)[vapply(sub, function(s)
    sum(s$tips %in% ref_ids) / n_ref_total > PHYLO_MAX_SPAN, logical(1))]
  broad <- list()
  if (length(bad)) {
    message("[R]   spans >", round(100 * PHYLO_MAX_SPAN),
            "% of the reference set, kept as lineage-level fallback: ",
            paste(sort(bad), collapse = ", "))
    broad <- sub[bad]
    sub <- sub[setdiff(names(sub), bad)]
  }
  low <- names(sub)[vapply(sub, function(s)
    !is.na(s$foreign_frac) && s$foreign_frac > PHYLO_MAX_FOREIGN, logical(1))]
  if (length(low)) {
    message("[R]   note: contains anchors of other clades (expected where clades nest): ",
            paste(low, collapse = ", "))
  }

  cc_anchored <- cl %>% filter(clade %in% c("CC-01", "CC-02"), role == "anchored") %>% nrow()
  ccs <- if (cc_anchored == 2) list(ok = FALSE, why = "CC-01/CC-02 are curated, not derived")
         else derive_cc_split(tree, sub, cl)
  if (ccs$ok) {
    sub[["CC-01"]] <- list(tips = ccs$cc01, n_tips = length(ccs$cc01), n_anchor = NA,
                           recovered = NA, foreign_frac = 0)
    sub[["CC-02"]] <- list(tips = ccs$cc02, n_tips = length(ccs$cc02), n_anchor = NA,
                           recovered = NA, foreign_frac = 0)
    message("[R]   CC split derived: ", ccs$note)
  } else if (identical(ccs$why, "CC-01/CC-02 are curated, not derived")) {
    message("[R]   CC split: from curated anchors")
  } else {
    message("[R]   WARNING: CC-01/CC-02 not derivable - ", ccs$why)
  }

  parent_of <- setNames(cl$parent, cl$clade)
  path_of <- function(nm) {
    p <- nm; out <- nm
    while (!is.na(parent_of[p]) && parent_of[p] != "root" && nzchar(parent_of[p])) {
      p <- parent_of[p]; out <- paste0(p, ">", out)
    }
    out
  }
  lineage_of <- function(nm) {
    p <- nm
    while (!is.na(parent_of[p]) && parent_of[p] != "root" && nzchar(parent_of[p])) {
      p <- parent_of[p]
    }
    p
  }

  ## most specific containing clade = smallest subtree that holds the tip
  order_by_size <- names(sub)[order(vapply(sub, function(s) s$n_tips, numeric(1)))]
  specific <- function(tp_) {
    for (nm in order_by_size) if (tp_ %in% sub[[nm]]$tips) return(nm)
    NA_character_
  }

  q_tips <- intersect(q$id, tree$tip.label)
  ref_by_protein <- split(ref_ids, tip_to_protein(ref_ids))
  input_ref_tips <- ref_ids[tip_to_protein(ref_ids) %in% nlr & !(ref_ids %in% q_tips)]
  all_tips <- c(q_tips, input_ref_tips)

  informative <- ref_ids[vapply(ref_ids, function(t) !is.na(specific(t)), logical(1))]
  needs_lineage <- all_tips[!(all_tips %in% names(outgroup_clade)) &
                              vapply(all_tips, function(t) is.na(specific(t)), logical(1))]
  lineage_near_map <- list()
  if (length(needs_lineage) > 0 && length(broad) > 0 && length(informative) > 0) {
    non_ref_need <- setdiff(needs_lineage, input_ref_tips)
    if (length(non_ref_need) > 0) {
      ln_df <- nearest_reference(tree, non_ref_need, informative, k = 1)
      for (i in seq_len(nrow(ln_df))) lineage_near_map[[ln_df$tip[i]]] <- ln_df$nearest_ref[i]
    }
    ref_need <- intersect(needs_lineage, input_ref_tips)
    for (tp_ in ref_need) {
      own_tips <- ref_by_protein[[tip_to_protein(tp_)]]
      cand_ref <- setdiff(informative, own_tips)
      if (length(cand_ref) > 0) {
        ln_df <- nearest_reference(tree, tp_, cand_ref, k = 1)
        if (nrow(ln_df) > 0) lineage_near_map[[tp_]] <- ln_df$nearest_ref[1]
      }
    }
  }

  assign_tip <- function(tp_) {
    if (tp_ %in% names(outgroup_clade)) {
      return(list(clade = outgroup_clade[[tp_]], basis = "outgroup"))
    }
    sp <- specific(tp_)
    if (!is.na(sp)) {
      return(list(clade = sp, basis = "clade"))
    }
    nr <- lineage_near_map[[tp_]]
    if (!is.null(nr) && !is.na(nr)) {
      sp_ref <- specific(nr)
      if (!is.na(sp_ref)) {
        return(list(clade = lineage_of(sp_ref), basis = "lineage fallback"))
      }
    }
    list(clade = NA_character_, basis = "not placed")
  }

  assigned <- lapply(all_tips, assign_tip)
  tip_clade <- tibble::tibble(
    tip   = all_tips,
    clade = vapply(assigned, `[[`, character(1), "clade"),
    basis = vapply(assigned, `[[`, character(1), "basis")
  )

  near_q <- nearest_reference(tree, q_tips, ref_ids)
  near_ref <- if (length(input_ref_tips) > 0) {
    ref_prots <- sort(unique(tip_to_protein(input_ref_tips)))
    message("[R]   ", length(ref_prots),
            " input protein(s) match scaffold identifiers, placed against remaining scaffold")
    purrr::map_dfr(ref_prots, function(p) {
      own <- ref_by_protein[[p]]
      others <- setdiff(ref_ids, own)
      nearest_reference(tree, own, others)
    })
  } else {
    nearest_reference(tree, character(), ref_ids)
  }
  near <- dplyr::bind_rows(near_q, near_ref)

  ## collapse to one row per protein: multiple NB-ARC regions must not expand the output
  per_tip <- tip_clade %>%
    left_join(near, by = "tip") %>%
    mutate(seqname = tip_to_protein(tip)) %>%
    arrange(seqname, dplyr::coalesce(nearest_dist, Inf), tip)

  per_protein <- per_tip %>%
    group_by(seqname) %>%
    summarise(
      n_nbarc_in_tree = dplyr::n(),
      clades = paste(sort(unique(stats::na.omit(clade))), collapse = "|"),
      n_distinct_clade = dplyr::n_distinct(stats::na.omit(clade)),
      basis = dplyr::if_else(all(basis == "outgroup"), "outgroup", ""),
      nearest_ref = nearest_ref[1],
      nearest_dist = nearest_dist[1],
      next_refs = next_refs[1],
      .groups = "drop") %>%
    mutate(clade = dplyr::na_if(clades, ""),
           clade = dplyr::coalesce(clade, "unplaced"),
           clade_path = vapply(clade, function(c) if (c %in% names(parent_of)) path_of(c) else c,
                               character(1)))

  n_nbarc_total <- regions %>% filter(region == "NB-ARC", seqname %in% nlr) %>%
    count(seqname, name = "n_nbarc")
  per_protein <- per_protein %>%
    left_join(n_nbarc_total, by = "seqname") %>%
    left_join(classification %>% select(seqname, arch_subclass = `Subclass (putative)`),
              by = "seqname") %>%
    mutate(category = categorise(arch_subclass, clade, clade_path,
                                 n_nbarc, n_distinct_clade - 1L, basis = basis)) %>%
    select(seqname, clade, clade_path, arch_subclass, category, nearest_ref, nearest_dist,
           next_refs, n_nbarc, n_nbarc_in_tree, n_distinct_clade)

  write_tsv_out(per_protein, file.path(pdir, "clade_assignment.tsv"))
  write_tsv_out(per_protein %>% filter(category != "agree"),
                file.path(pdir, "discrepancies.tsv"))
  write_tsv_out(per_protein %>% filter(clade == "unplaced"), file.path(pdir, "unplaced.tsv"))
  if (nrow(al$dropped)) write_tsv_out(al$dropped, file.path(pdir, "excluded_low_occupancy.tsv"))

  summ <- tibble::tibble(clade = names(sub),
                         n_tips = vapply(sub, function(s) s$n_tips, numeric(1)),
                         n_anchor = vapply(sub, function(s) s$n_anchor %||% NA_real_, numeric(1)),
                         anchor_recovery = vapply(sub, function(s) s$recovered %||% NA_real_,
                                                  numeric(1)),
                         foreign_frac = vapply(sub, function(s) s$foreign_frac, numeric(1)))
  write_tsv_out(summ, file.path(pdir, "clade_summary.tsv"))

  message("[R]   assigned ", nrow(per_protein), " proteins")
  per_protein %>% count(category, sort = TRUE) %>% as.data.frame() %>% print(row.names = FALSE)

  per_protein
}
