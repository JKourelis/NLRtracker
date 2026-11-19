## 07_report.R -- self-contained interactive HTML report.
##
## Writes a compact JSON payload into a static HTML/JS shell that draws protein schematics
## client-side as inline SVG. Pre-rendering thousands of SVGs server-side would make the file
## unusable on a real proteome; this keeps it small and makes filtering instant.
##
## No external assets, no CDN, no fonts. One file the user can open or email.

## ------------------------------------------------------------------ payload

## Compact keys, because this is repeated once per protein:
##   i id · n length · s Status · c Subclass · u support · a architecture · m simplified
##   r regions [name,start,end] · d specific integrated domains · f motifs
build_report_payload <- function(classification, domains_tbl, regions, annotation, cfg, meta,
                                 phylo = NULL) {
  chains <- domains_tbl %>% filter(type == "CHAIN") %>%
    select(seqname, len = end) %>% distinct(seqname, .keep_all = TRUE)

  cls <- classification %>%
    left_join(chains, by = "seqname") %>%
    mutate(len = dplyr::coalesce(len, 0L))

  reg <- domains_tbl %>% filter(type == "REGION") %>%
    select(seqname, description, start, end)
  dom <- domains_tbl %>% filter(type == "DOMAIN") %>%
    select(seqname, description, start, end)
  mot <- domains_tbl %>% filter(type == "MOTIF") %>%
    select(seqname, description, start, end)

  split_by <- function(df) split(df[, c("description", "start", "end")], df$seqname)
  reg_l <- split_by(reg); dom_l <- split_by(dom); mot_l <- split_by(mot)

  as_triples <- function(d) {
    if (is.null(d) || !nrow(d)) return(list())
    lapply(seq_len(nrow(d)), function(i)
      list(d$description[i], as.integer(d$start[i]), as.integer(d$end[i])))
  }

  ph <- NULL
  if (!is.null(phylo) && nrow(phylo) > 0) {
    ph <- phylo %>% select(any_of(c("seqname", "clade", "clade_path",
                                    "nearest_ref", "nearest_dist", "category")))
  }

  proteins <- lapply(seq_len(nrow(cls)), function(k) {
    id <- cls$seqname[k]
    out <- list(
      i = id, n = as.integer(cls$len[k]),
      s = cls$Status[k], c = cls$`Subclass (putative)`[k],
      u = dplyr::coalesce(cls$subclass_support[k], ""),
      a = cls$Domain[k], m = cls$Simple[k],
      r = as_triples(reg_l[[id]]), d = as_triples(dom_l[[id]]), f = as_triples(mot_l[[id]])
    )
    if (!is.null(ph)) {
      row <- ph[ph$seqname == id, , drop = FALSE]
      if (nrow(row)) {
        out$p <- as.list(row[1, setdiff(names(row), "seqname"), drop = FALSE])
      }
    }
    out
  })

  ## region palette, taken from the class table so the HTML and the PDFs agree
  pal <- cfg$classes %>%
    filter(!is.na(region), region != "-", !is.na(colour), colour != "") %>%
    distinct(region, .keep_all = TRUE) %>%
    select(region, colour)
  palette <- setNames(as.list(pal$colour), pal$region)
  palette[["ID"]] <- palette[["ID"]] %||% "#666666"

  list(meta = meta, palette = palette, proteins = proteins,
       has_phylo = !is.null(ph))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

## ------------------------------------------------------------------ HTML shell
##
## The shell lives in module/report_template.html rather than inline here. Embedding a page
## of CSS and JavaScript in an R string means escaping every quote and backslash, which is
## unreadable and was in fact the first thing to break. Keeping it as a real .html file also
## means it can be edited with syntax highlighting and checked with `node --check`.

report_html_shell <- function(payload_json, title, template) {
  if (!file.exists(template)) {
    stop("report template not found: ", template, call. = FALSE)
  }
  html <- paste(readLines(template, warn = FALSE), collapse = "\n")
  html <- sub("__TITLE__", title, html, fixed = TRUE)
  ## payload substituted last, and with fixed = TRUE, so backslashes and $ in protein
  ## identifiers cannot be interpreted as regex replacements
  sub("__PAYLOAD__", payload_json, html, fixed = TRUE)
}

## ------------------------------------------------------------------ entry point

write_report <- function(outdir, classification, domains_tbl, regions, annotation, cfg, meta,
                         template, phylo = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    message("[R]   jsonlite not installed - skipping report.html")
    return(invisible(NULL))
  }
  payload <- build_report_payload(classification, domains_tbl, regions, annotation, cfg,
                                  meta, phylo)
  ## auto_unbox so scalars are not wrapped in arrays; digits = NA keeps numbers exact
  js <- as.character(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null",
                                      na = "null", digits = NA))
  ## A literal </script> anywhere in the payload -- a protein identifier could contain
  ## one -- would terminate the script block early and corrupt the page.
  js <- gsub("</", "<\\/", js, fixed = TRUE)
  path <- file.path(outdir, "report.html")
  writeLines(report_html_shell(js, paste0("NLRtracker \u2014 ",
                                         dplyr::coalesce(meta$input, "report")), template), path)
  message("[R]   report.html  (", round(file.size(path) / 1024), " KB, ",
          nrow(classification), " proteins)")
  invisible(path)
}
