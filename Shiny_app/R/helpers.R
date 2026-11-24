# ============================================================
# helpers.R  —  Data logic for the Limma DEG Explorer (R Shiny)
#
# Rule 3: everything that computes lives here, outside the Shiny
# server, so the same functions drive the on-screen plots, the
# downloadable report and the tests.
#
# Author:      Patricia Agudelo-Romero
# Institution: The Kids Research Institute Australia
# ============================================================

APP_VERSION <- "1.0.0"

REQUIRED_COLUMNS <- c("gene", "logFC", "adj.P.Val")
TABLE_COLUMNS    <- c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
GO_TABLE_COLUMNS <- c("Description", "GeneRatio", "BgRatio", "pvalue",
                      "p.adjust", "Count", "geneID")
COLOURS <- c(NS = "grey70", Up = "#7B2D8B", Down = "#2E8B57")
GO_ONTOLOGIES <- c(GO_BP = "BP", GO_MF = "MF", GO_CC = "CC")
GO_LABELS <- c(GO_BP = "GO Biological Process",
               GO_MF = "GO Molecular Function",
               GO_CC = "GO Cellular Component")
MIN_DEGS_FOR_ENRICHMENT <- 5


# ── Rule 4: validate at the boundary ─────────────────────────

validate_limma_table <- function(df) {
  missing <- setdiff(REQUIRED_COLUMNS, names(df))
  if (length(missing) > 0) {
    return(sprintf("Missing required column(s): %s. Required: %s.",
                   paste(missing, collapse = ", "),
                   paste(REQUIRED_COLUMNS, collapse = ", ")))
  }
  problems <- character(0)
  if (nrow(df) == 0) problems <- c(problems, "The uploaded file contains no rows.")
  for (col in c("logFC", "adj.P.Val")) {
    if (!is.numeric(df[[col]])) {
      problems <- c(problems, sprintf("Column '%s' must be numeric.", col))
    }
  }
  if (length(problems) == 0) {
    p <- df$adj.P.Val[!is.na(df$adj.P.Val)]
    if (any(p < 0 | p > 1)) {
      problems <- c(problems, "Column 'adj.P.Val' must contain values between 0 and 1.")
    }
  }
  problems
}

# Returns list(data = data.frame or NULL, problems = character())
read_limma_csv <- function(path) {
  df <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
                 error = function(e) NULL)
  if (is.null(df)) {
    return(list(data = NULL, problems = "Could not read the file as CSV."))
  }
  problems <- validate_limma_table(df)
  if (length(problems) > 0) return(list(data = NULL, problems = problems))
  df <- df[stats::complete.cases(df[, REQUIRED_COLUMNS]), , drop = FALSE]
  df$gene <- as.character(df$gene)
  list(data = df, problems = character(0))
}


# ── DEG filtering ────────────────────────────────────────────

filter_degs <- function(df, pval, lfc) {
  df[df$adj.P.Val <= pval & abs(df$logFC) >= lfc, , drop = FALSE]
}

deg_table <- function(degs) {
  out <- degs[order(degs$adj.P.Val), intersect(TABLE_COLUMNS, names(degs)), drop = FALSE]
  rownames(out) <- NULL
  out
}


# ── GO enrichment (clusterProfiler + org.Mm.eg.db, offline) ──

# Results are cached per (genes, universe, ontology), so moving the
# enrichment P-value slider re-filters instantly without recomputing.
.go_cache <- new.env(parent = emptyenv())

clear_go_cache <- function() rm(list = ls(.go_cache), envir = .go_cache)

go_cache_key <- function(genes, universe, ont) {
  rlang::hash(list(sort(unique(genes)), sort(unique(universe)), ont))
}

is_go_cached <- function(genes, universe, ont) {
  exists(go_cache_key(genes, universe, ont), envir = .go_cache, inherits = FALSE)
}

# One unfiltered enrichGO run. The background (universe) is every gene in
# the uploaded table, not the whole genome. pvalueCutoff and qvalueCutoff
# are set to 1 so the app's slider is the only significance filter.
.enrich_go <- function(genes, universe, ont) {
  db  <- org.Mm.eg.db::org.Mm.eg.db
  map <- suppressMessages(suppressWarnings(
    clusterProfiler::bitr(unique(universe), fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = db)))
  gene_ids <- unique(map$ENTREZID[map$SYMBOL %in% genes])
  if (length(gene_ids) == 0) {
    stop("none of the gene symbols could be mapped to mouse Entrez IDs")
  }
  ego <- clusterProfiler::enrichGO(
    gene = gene_ids, universe = unique(map$ENTREZID), OrgDb = db,
    ont = ont, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
    readable = TRUE)
  if (is.null(ego)) return(data.frame(p.adjust = numeric(0)))
  as.data.frame(ego)
}

# Returns list(table = data.frame or NULL, error = character or NULL)
run_go <- function(genes, universe, ont, cutoff) {
  if (length(unique(genes)) < MIN_DEGS_FOR_ENRICHMENT) {
    return(list(table = NULL, error = sprintf(
      "At least %d DEGs are needed for enrichment analysis.", MIN_DEGS_FOR_ENRICHMENT)))
  }
  key <- go_cache_key(genes, universe, ont)
  if (!exists(key, envir = .go_cache, inherits = FALSE)) {
    res <- tryCatch(.enrich_go(genes, universe, ont), error = function(e) e)
    if (inherits(res, "error")) {
      return(list(table = NULL,
                  error = paste("GO enrichment could not be run:", conditionMessage(res))))
    }
    assign(key, res, envir = .go_cache)
  }
  res <- get(key, envir = .go_cache)
  res <- res[res$p.adjust <= cutoff, , drop = FALSE]
  res <- res[order(res$p.adjust), , drop = FALSE]
  rownames(res) <- NULL
  list(table = if (nrow(res) > 0) res else NULL, error = NULL)
}

go_table <- function(df) df[, intersect(GO_TABLE_COLUMNS, names(df)), drop = FALSE]


# ── Figures (shared by the app and the report) ───────────────

volcano_plot <- function(df, pval, lfc, n_labels = 5) {
  sig <- df$adj.P.Val <= pval & abs(df$logFC) >= lfc
  df$direction <- factor(ifelse(!sig, "NS", ifelse(df$logFC > 0, "Up", "Down")),
                         levels = c("NS", "Up", "Down"))
  counts <- table(df$direction)
  labels <- sprintf("%s (%d)", names(counts), as.integer(counts))
  up  <- df[sig & df$logFC > 0, ]
  dn  <- df[sig & df$logFC < 0, ]
  top <- rbind(utils::head(up[order(up$adj.P.Val), ], n_labels),
               utils::head(dn[order(dn$adj.P.Val), ], n_labels))
  ggplot2::ggplot(df, ggplot2::aes(logFC, -log10(pmax(adj.P.Val, 1e-300)),
                                   colour = direction)) +
    ggplot2::geom_point(alpha = 0.55, size = 1.2) +
    ggplot2::scale_colour_manual(values = COLOURS, labels = labels, drop = FALSE) +
    ggplot2::geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(pval), linetype = "dashed", linewidth = 0.4) +
    ggrepel::geom_text_repel(data = top, ggplot2::aes(label = gene), size = 3,
                             colour = "black", max.overlaps = 20,
                             segment.colour = "grey50", show.legend = FALSE) +
    ggplot2::labs(title = "Volcano plot", x = expression(log[2] ~ "fold change"),
                  y = expression(-log[10] ~ "(adj. P-value)"), colour = NULL) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "top")
}

rank_barplot <- function(degs, n = 30) {
  degs$rank <- sign(degs$logFC) * -log10(pmax(degs$adj.P.Val, 1e-300))
  top <- utils::head(degs[order(-abs(degs$rank)), ], n)
  top$direction <- ifelse(top$logFC > 0, "Up", "Down")
  ggplot2::ggplot(top, ggplot2::aes(stats::reorder(gene, rank), rank, fill = direction)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = COLOURS[c("Up", "Down")]) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = sprintf("Top %d DEGs by rank", nrow(top)), x = NULL,
                  y = expression("Rank score: sign(" * log[2] * "FC)" %*%
                                   -log[10] * "(adj. P-value)"),
                  fill = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")
}

go_dotplot <- function(df, title, n = 15) {
  top <- utils::head(df[order(df$p.adjust), ], n)
  top$Description <- ifelse(nchar(top$Description) > 55,
                            paste0(substr(top$Description, 1, 54), "\u2026"),
                            top$Description)
  ggplot2::ggplot(top, ggplot2::aes(Count, stats::reorder(Description, -p.adjust),
                                    colour = p.adjust, size = Count)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low = "#B2182B", high = "#2166AC") +
    ggplot2::labs(title = title, x = "Gene count", y = NULL, colour = "adj. P-value") +
    ggplot2::theme_minimal(base_size = 12)
}

message_plot <- function(text) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = text, colour = "grey40", size = 5) +
    ggplot2::theme_void()
}


# ── Rule 6: report folder ────────────────────────────────────

pkg_version <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "not installed")
}

# Writes DEG_analysis_<stamp>.zip to `file`, containing the HTML report,
# the DEG table (CSV) and GO enrichment for the three ontologies (Excel).
build_report_zip <- function(file, raw, degs, params, go_results,
                             template, now = Sys.time()) {
  stamp   <- format(now, "%Y%m%d_%H%M%S")
  out_dir <- tempfile("deg_report_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # 1. DEG table
  utils::write.csv(deg_table(degs), file.path(out_dir, paste0("DEG_results_", stamp, ".csv")),
                   row.names = FALSE)

  # 2. GO workbook
  wb <- openxlsx::createWorkbook()
  header <- openxlsx::createStyle(fontName = "Arial", fontSize = 11, fontColour = "white",
                                  fgFill = "#2E74B5", halign = "left", textDecoration = "bold")
  for (key in names(GO_ONTOLOGIES)) {
    tab <- go_results[[key]]$table
    if (is.null(tab)) next
    tab <- go_table(tab)
    openxlsx::addWorksheet(wb, key)
    openxlsx::writeData(wb, key, tab, headerStyle = header)
    openxlsx::setColWidths(wb, key, cols = seq_along(tab), widths = "auto")
    openxlsx::freezePane(wb, key, firstRow = TRUE)
  }
  if (length(openxlsx::sheets(wb)) == 0) {
    openxlsx::addWorksheet(wb, "No_results")
    errors <- unlist(lapply(go_results, `[[`, "error"))
    msg <- c(sprintf("No enriched GO terms at adj. P-value <= %s", params$enr_pval),
             if (length(errors) > 0) paste0(names(errors), ": ", errors))
    openxlsx::writeData(wb, "No_results", data.frame(Message = msg))
  }
  openxlsx::saveWorkbook(wb, file.path(out_dir, paste0("GO_enrichment_", stamp, ".xlsx")),
                         overwrite = TRUE)

  # 3. HTML report, rendered from a temporary copy of the template so
  #    nothing is written into the app folder (read-only on servers)
  work_dir <- tempfile("rmd_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  rmd <- file.path(work_dir, "report.Rmd")
  file.copy(template, rmd)
  rmarkdown::render(
    rmd,
    output_file = paste0("DEG_report_", stamp, ".html"),
    output_dir = out_dir,
    intermediates_dir = work_dir,
    params = list(raw = raw, degs = degs, deg_pval = params$deg_pval,
                  deg_lfc = params$deg_lfc, enr_pval = params$enr_pval,
                  input_name = params$input_name, go_results = go_results,
                  generated = format(now, "%d %B %Y, %H:%M")),
    envir = new.env(parent = parent.env(environment())),
    quiet = TRUE
  )

  # 4. Bundle (files at the root of the zip)
  zip::zipr(file, list.files(out_dir, full.names = TRUE))
  invisible(file)
}
