## 06_draw.R -- domain-architecture figures.
##
## Ported from refplantnlR::draw_nlr() and kept in sync with it. Fixes applied here and
## mirrored in the package:
##
##   * ylim()/xlim() followed by scale_x_continuous() silently discarded the limits and
##     warned. Now one coord_cartesian() call.
##   * geom_rect(size=) is deprecated since ggplot2 3.4; use linewidth.
##   * scale_fill_manual() was a closed 9-value set, so any new region name aborted with
##     "Insufficient values in manual scale". Now a named palette with a grey fallback,
##     which matters because the region vocabulary gained ID and C-JID.
##   * max(data$end) in the ggsave sizing lacked na.rm.
##   * Height scaled linearly with protein count, so 481 NLRs produced an unusable
##     figure. Now paginated.
##   * PDF = "file.pdf" overloaded a logical as a filename; use file = and pick the
##     device from the extension.

REGION_PALETTE <- c(
  "R1"     = "#FBB4AE",
  "CC"     = "#B3CDE3",
  "RPW8"   = "#CCEBC5",
  "TIR"    = "#DECBE4",
  "PLOOP"  = "#FFFFCC",
  "NB-ARC" = "#FED9A6",
  "LRR"    = "#E5D8BD",
  "C-JID"  = "#FDDAEC",
  "MLKL"   = "#CCEBC5",
  "ID"     = "#666666"
)
REGION_FALLBACK <- "#999999"

palette_for <- function(values) {
  values <- unique(stats::na.omit(values))
  pal <- REGION_PALETTE[values]
  pal[is.na(pal)] <- REGION_FALLBACK
  names(pal) <- values
  pal
}

draw_nlr <- function(data,
                     label_regions = TRUE,
                     label_domains = TRUE,
                     label_motifs  = FALSE,
                     label_repeats = FALSE,
                     order_by = c("subclass", "name", "length"),
                     file = NULL,
                     per_page = 50,
                     width = NULL, height = NULL) {
  stopifnot(is.data.frame(data))
  required <- c("seqname", "type", "description", "start", "end")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("draw_nlr(): data is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (!any(data$type == "CHAIN")) {
    stop("draw_nlr(): data must contain at least one CHAIN row", call. = FALSE)
  }
  order_by <- match.arg(order_by)

  chains <- data %>% filter(type == "CHAIN")
  ord <- switch(
    order_by,
    name     = chains %>% arrange(desc(seqname)),
    length   = chains %>% arrange(end),
    subclass = {
      if ("Subclass (putative)" %in% names(chains)) {
        chains %>% arrange(desc(`Subclass (putative)`), desc(seqname))
      } else {
        chains %>% arrange(desc(seqname))
      }
    }
  ) %>% distinct(seqname) %>% mutate(.order = row_number())

  data <- data %>% inner_join(ord, by = "seqname")

  ## Rename the NB-ARC motifs to their canonical order number, as the package does.
  motif_number <- setNames(as.character(NBARC_MOTIF_ORDER), names(NBARC_MOTIF_ORDER))
  data <- data %>%
    mutate(description = dplyr::if_else(
      type == "MOTIF" & description %in% names(motif_number),
      unname(motif_number[description]), description))

  n_prot <- nrow(ord)
  pages <- if (is.null(per_page) || per_page <= 0) 1L else
    as.integer(ceiling(n_prot / per_page))
  max_end <- max(data$end, na.rm = TRUE)

  plots <- vector("list", pages)
  for (p in seq_len(pages)) {
    lo <- (p - 1L) * per_page + 1L
    hi <- min(p * per_page, n_prot)
    d <- data %>% filter(.order >= lo, .order <= hi) %>%
      mutate(row = .order - lo + 1L)
    n_row <- hi - lo + 1L

    regions <- d %>% filter(type == "REGION")
    pal <- palette_for(regions$description)

    g <- ggplot2::ggplot() +
      ggplot2::labs(x = "Amino acid number", y = NULL, fill = "Region") +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      ggplot2::coord_cartesian(
        xlim = c(-max_end * 0.22, max_end * 1.05),
        ylim = c(0.4, n_row + 0.6), expand = FALSE) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        legend.key.size = ggplot2::unit(0.4, "cm"),
        legend.position = "right")

    ## chain backbone
    g <- g + ggplot2::geom_rect(
      data = d %>% filter(type == "CHAIN"),
      ggplot2::aes(xmin = start, xmax = end, ymin = row - 0.12, ymax = row + 0.12),
      fill = "grey75", linewidth = 0)

    ## chain labels
    ch <- d %>% filter(type == "CHAIN")
    g <- g + ggplot2::annotate("text", x = -max_end * 0.01, y = ch$row,
                               label = ch$seqname, hjust = 1, size = 2.6)

    ## regions: draw NB-ARC last so it sits on top where boundaries touch
    if (nrow(regions)) {
      back <- regions %>% filter(!description %in% c("NB-ARC", "PLOOP"))
      front <- regions %>% filter(description %in% c("NB-ARC", "PLOOP"))
      if (nrow(back)) {
        g <- g + ggplot2::geom_rect(
          data = back,
          ggplot2::aes(xmin = start, xmax = end, ymin = row - 0.3, ymax = row + 0.3,
                       fill = description), linewidth = 0)
      }
      if (nrow(front)) {
        g <- g + ggplot2::geom_rect(
          data = front,
          ggplot2::aes(xmin = start, xmax = end, ymin = row - 0.3, ymax = row + 0.3,
                       fill = description), linewidth = 0)
      }
      g <- g + ggplot2::scale_fill_manual(values = pal, drop = FALSE)
      if (label_regions) {
        g <- g + ggplot2::geom_text(
          data = regions %>% filter(end - start > max_end * 0.02),
          ggplot2::aes(x = start + (end - start) / 2, y = row, label = description),
          colour = "white", fontface = "bold", size = 1.9)
      }
    }

    ## integrated-domain labels keep the specific domain name
    dom <- d %>% filter(type == "DOMAIN")
    if (nrow(dom) && label_domains) {
      g <- g + ggplot2::geom_label(
        data = dom,
        ggplot2::aes(x = start + (end - start) / 2, y = row, label = description),
        fill = "#666666", colour = "white", fontface = "bold",
        size = 1.8, linewidth = 0, label.padding = ggplot2::unit(0.08, "lines"))
    }

    mot <- d %>% filter(type == "MOTIF")
    if (nrow(mot) && label_motifs) {
      g <- g + ggplot2::geom_rect(
        data = mot,
        ggplot2::aes(xmin = start, xmax = end, ymin = row - 0.46, ymax = row - 0.26),
        fill = "#c1c1c1", linewidth = 0) +
        ggplot2::geom_text(
          data = mot,
          ggplot2::aes(x = start + (end - start) / 2, y = row - 0.36,
                       label = description), size = 1.5)
    }

    rep <- d %>% filter(type == "REPEAT")
    if (nrow(rep) && label_repeats) {
      g <- g + ggplot2::geom_rect(
        data = rep,
        ggplot2::aes(xmin = start, xmax = end, ymin = row - 0.46, ymax = row - 0.26),
        fill = "#c1c1c1", linewidth = 0)
    }

    if (pages > 1) {
      g <- g + ggplot2::labs(subtitle = sprintf("page %d of %d", p, pages))
    }
    plots[[p]] <- g
  }

  if (!is.null(file)) {
    w <- if (is.null(width))  min(24, max(8, max_end / 120)) else width
    h <- if (is.null(height)) min(48, max(3, per_page * 0.22)) else height
    ext <- tolower(tools::file_ext(file))
    if (ext == "pdf") {
      grDevices::pdf(file, width = w, height = h, onefile = TRUE)
      for (g in plots) print(g)
      grDevices::dev.off()
    } else {
      for (p in seq_along(plots)) {
        out <- if (length(plots) == 1) file else
          sub(paste0("\\.", ext, "$"), sprintf("_p%02d.%s", p, ext), file)
        ggplot2::ggsave(out, plots[[p]], width = w, height = h,
                        limitsize = FALSE, dpi = 150)
      }
    }
  }

  invisible(if (length(plots) == 1) plots[[1]] else plots)
}

## Figures produced automatically by the pipeline. Motifs and repeats are off, per the
## agreed default.
draw_all_figures <- function(domains_tbl, outdir, nlr_ids, assoc_ids) {
  figdir <- file.path(outdir, "figures")
  dir.create(file.path(figdir, "by_subclass"), showWarnings = FALSE, recursive = TRUE)

  safe_draw <- function(d, path, ...) {
    if (nrow(d) == 0 || !any(d$type == "CHAIN")) return(invisible(NULL))
    tryCatch(draw_nlr(d, label_motifs = FALSE, label_repeats = FALSE, file = path, ...),
             error = function(e) message("  figure failed for ", basename(path), ": ",
                                         conditionMessage(e)))
  }

  safe_draw(domains_tbl %>% filter(seqname %in% nlr_ids),
            file.path(figdir, "NLR.pdf"))
  safe_draw(domains_tbl %>% filter(seqname %in% assoc_ids),
            file.path(figdir, "NLR-associated.pdf"))

  if ("Subclass (putative)" %in% names(domains_tbl)) {
    subs <- domains_tbl %>%
      filter(seqname %in% nlr_ids, !is.na(`Subclass (putative)`)) %>%
      distinct(`Subclass (putative)`) %>% pull()
    for (s in subs) {
      slug <- gsub("[^A-Za-z0-9]+", "_", s)
      safe_draw(domains_tbl %>% filter(seqname %in% nlr_ids,
                                       `Subclass (putative)` == s),
                file.path(figdir, "by_subclass", paste0(slug, ".pdf")))
    }
  }
  invisible(NULL)
}
