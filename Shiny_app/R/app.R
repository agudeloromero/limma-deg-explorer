# ============================================================
# app.R  —  Limma DEG Interactive Explorer (R Shiny)
#
# Author:      Patricia Agudelo-Romero
# Email:       patricia.agudeloromero@thekids.org.au
# Institution: The Kids Research Institute Australia
#
# Rules implemented in this file:
#   Rule 2  — Framework chosen for Bioconductor ecosystem fit
#   Rule 3  — Reactive expressions separate data from display
#   Rule 4  — Input validation at the boundary (raw_data reactive)
#   Rule 5  — Every threshold is a UI control (sliders / radio)
#   Rule 6  — Analysis state captured in downloadable report folder
#   Rule 9  — App is deployed to a URL (see 02_deploy.R)
#   Rule 10 — UI surfaces decisions, not implementation details
#
# Download output (Rule 6):
#   DEG_analysis_TIMESTAMP.zip containing:
#     ├── DEG_report_TIMESTAMP.html     (figures only)
#     ├── DEG_results_TIMESTAMP.csv     (full DEG table)
#     └── GO_enrichment_TIMESTAMP.xlsx  (GO_BP / GO_MF / GO_CC sheets)
#
# Workflow: Upload CSV → Filter DEGs → Explore plots →
#           Filter enrichment → Download report folder
# ============================================================

library(shiny)
library(ggplot2)
library(ggrepel)          # Rule 2: no Python equivalent; R-specific advantage
library(DT)
library(dplyr)
library(rmarkdown)
library(clusterProfiler)  # Rule 2: Bioconductor — offline, no API call
library(org.Mm.eg.db)
library(enrichplot)
library(openxlsx)         # Rule 6: multi-sheet Excel for GO enrichment output


# ── UI ───────────────────────────────────────────────────────────────────────

ui <- fluidPage(

  titlePanel("Limma DEG Interactive Explorer"),

  tags$head(tags$style(HTML("
    body { font-family: Arial, sans-serif; }
    .section-label { font-weight: bold; color: #2C3E50; margin-top: 8px; }
    .num-box .form-group { margin-bottom: 0; }
    .num-box input { text-align: center; padding: 4px 6px; }
  "))),

  sidebarLayout(

    sidebarPanel(width = 3,

      fileInput("file", "Upload limma output (CSV)",
                accept      = ".csv",
                placeholder = "Required columns: gene, logFC, adj.P.Val"),

      hr(),

      tags$h5("DEG Filters", style = "color:#2E74B5; font-weight:bold;"),

      fluidRow(
        column(8, sliderInput("deg_pval", "Adjusted P-value \u2264",
                              min = 0.001, max = 1, value = 0.05, step = 0.001)),
        column(4, tags$div(class = "num-box", tags$br(),
                           numericInput("deg_pval_num", NULL,
                                        min = 0.001, max = 1,
                                        value = 0.05, step = 0.001)))
      ),

      fluidRow(
        column(8, sliderInput("deg_lfc", "| log\u2082FC | \u2265",
                              min = 0, max = 4, value = 1.0, step = 0.1)),
        column(4, tags$div(class = "num-box", tags$br(),
                           numericInput("deg_lfc_num", NULL,
                                        min = 0, max = 4,
                                        value = 1.0, step = 0.1)))
      ),

      hr(),

      tags$h5("Enrichment Filters", style = "color:#1E8449; font-weight:bold;"),

      selectInput("geneset", "Gene set (interactive preview):",
                  choices = c(
                    "GO Biological Process" = "GO_BP",
                    "GO Molecular Function" = "GO_MF",
                    "GO Cellular Component" = "GO_CC"
                  )),

      tags$p(style = "font-size:11px; color:#888; margin-top:-6px;",
             "All three ontologies are included in the downloaded report."),

      fluidRow(
        column(8, sliderInput("enr_pval", "Enrichment adj. P-value \u2264",
                              min = 0.001, max = 1, value = 0.05, step = 0.001)),
        column(4, tags$div(class = "num-box", tags$br(),
                           numericInput("enr_pval_num", NULL,
                                        min = 0.001, max = 1,
                                        value = 0.05, step = 0.001)))
      ),

      hr(),

      tags$h5("Report", style = "color:#2C3E50; font-weight:bold;"),
      tags$p(style = "font-size:11px; color:#888;",
             "Downloads a .zip folder containing the HTML report (figures),
              DEG CSV, and GO enrichment Excel (all three ontologies)."),
      downloadButton("dl_report", "Download Report Folder",
                     style = "width:100%; background:#2E74B5; color:white;
                              border:none; padding:8px; border-radius:4px;
                              cursor:pointer;")
    ),

    mainPanel(width = 9,
      uiOutput("validation_msg"),
      uiOutput("deg_summary"),
      uiOutput("welcome_msg"),
      uiOutput("main_tabs")
    )
  )
)


# ── Server ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Slider ↔ numeric synchronisation ────────────────────────
  observeEvent(input$deg_pval,
    updateNumericInput(session, "deg_pval_num", value = input$deg_pval),
    ignoreInit = TRUE)
  observeEvent(input$deg_pval_num, {
    req(input$deg_pval_num)
    updateSliderInput(session, "deg_pval", value = input$deg_pval_num)
  }, ignoreInit = TRUE)

  observeEvent(input$deg_lfc,
    updateNumericInput(session, "deg_lfc_num", value = input$deg_lfc),
    ignoreInit = TRUE)
  observeEvent(input$deg_lfc_num, {
    req(input$deg_lfc_num)
    updateSliderInput(session, "deg_lfc", value = input$deg_lfc_num)
  }, ignoreInit = TRUE)

  observeEvent(input$enr_pval,
    updateNumericInput(session, "enr_pval_num", value = input$enr_pval),
    ignoreInit = TRUE)
  observeEvent(input$enr_pval_num, {
    req(input$enr_pval_num)
    updateSliderInput(session, "enr_pval", value = input$enr_pval_num)
  }, ignoreInit = TRUE)


  # ==========================================================================
  # Rule 3: Data pipeline
  #   raw_data()  →  degs()  →  enr_result()  →  enr_filtered()
  # ==========================================================================

  raw_data <- reactive({
    req(input$file)
    df       <- read.csv(input$file$datapath, stringsAsFactors = FALSE)
    required <- c("gene", "logFC", "AveExpr", "P.Value", "adj.P.Val")
    missing  <- setdiff(required, names(df))
    validate(need(length(missing) == 0,
                  paste("Missing required columns:", paste(missing, collapse = ", "))))
    validate(need(nrow(df) > 0,             "Uploaded file contains no rows."))
    validate(need(is.numeric(df$logFC),     "Column 'logFC' must be numeric."))
    validate(need(is.numeric(df$adj.P.Val), "Column 'adj.P.Val' must be numeric."))
    df
  })

  degs <- reactive({
    df <- raw_data()
    df[df$adj.P.Val <= input$deg_pval & abs(df$logFC) >= input$deg_lfc, ]
  })

  # Interactive preview: runs only the selected ontology
  enr_result <- reactive({
    validate(need(nrow(degs()) >= 5,
                  "Minimum 5 DEGs required for enrichment analysis."))
    gene_symbols <- degs()$gene
    ids <- tryCatch(
      clusterProfiler::bitr(gene_symbols, fromType = "SYMBOL",
                             toType = "ENTREZID", OrgDb = org.Mm.eg.db),
      error = function(e) NULL
    )
    validate(need(!is.null(ids) && nrow(ids) > 0,
                  "Could not map gene symbols to Entrez IDs."))
    entrez <- ids$ENTREZID
    withProgress(message = "Running GO enrichment...",
                 detail  = "This may take a few seconds.", value = 0, {
      incProgress(0.4, detail = "Mapping identifiers...")
      result <- switch(input$geneset,
        "GO_BP" = enrichGO(gene = entrez, OrgDb = org.Mm.eg.db,
                           ont = "BP", pvalueCutoff = 1, readable = TRUE),
        "GO_MF" = enrichGO(gene = entrez, OrgDb = org.Mm.eg.db,
                           ont = "MF", pvalueCutoff = 1, readable = TRUE),
        "GO_CC" = enrichGO(gene = entrez, OrgDb = org.Mm.eg.db,
                           ont = "CC", pvalueCutoff = 1, readable = TRUE)
      )
      incProgress(0.6, detail = "Done.")
      result
    })
  })

  enr_filtered <- reactive({
    res <- as.data.frame(enr_result())
    res[res[["p.adjust"]] <= input$enr_pval, ]
  })

  # Helper: run one GO ontology for the download (all three)
  .run_go <- function(gene_symbols, ont, pval_cutoff) {
    ids <- tryCatch(
      clusterProfiler::bitr(gene_symbols, fromType = "SYMBOL",
                             toType = "ENTREZID", OrgDb = org.Mm.eg.db),
      error = function(e) NULL
    )
    if (is.null(ids) || nrow(ids) == 0) return(NULL)
    res <- tryCatch(
      enrichGO(gene = ids$ENTREZID, OrgDb = org.Mm.eg.db,
               ont = ont, pvalueCutoff = 1, readable = TRUE),
      error = function(e) NULL
    )
    if (is.null(res)) return(NULL)
    df <- as.data.frame(res)
    df <- df[df$p.adjust <= pval_cutoff, ]
    if (nrow(df) == 0) return(NULL)
    df
  }


  # ==========================================================================
  # Outputs
  # ==========================================================================

  output$validation_msg <- renderUI({ req(input$file); NULL })

  output$welcome_msg <- renderUI({
    if (!is.null(input$file)) return(NULL)
    tags$div(
      style = "margin-top:60px; text-align:center; color:#888;",
      tags$img(src = "https://www.r-project.org/logo/Rlogo.png",
               height = "48px", style = "opacity:0.3; margin-bottom:16px;"),
      tags$h4("Welcome to the Limma DEG Interactive Explorer",
              style = "color:#2E74B5;"),
      tags$p("To get started, upload your limma output CSV using the panel on the left."),
      tags$p(style = "font-size:13px;",
             "Required columns: ",
             tags$code("gene"), ", ", tags$code("logFC"), ", ",
             tags$code("AveExpr"), ", ", tags$code("P.Value"), ", ",
             tags$code("adj.P.Val")),
      tags$p(style = "font-size:12px; color:#aaa;",
             "This is the default output of ",
             tags$code("limma::topTable()"),
             " \u2014 no reformatting needed.")
    )
  })

  output$main_tabs <- renderUI({
    req(input$file)
    req(raw_data())
    tabsetPanel(
      tabPanel("DEG Table",        DTOutput("deg_table")),
      tabPanel("Volcano Plot",     plotOutput("volcano",    height = "460px")),
      tabPanel("Top Ranked Genes", plotOutput("barplot",    height = "460px")),
      tabPanel("GO Enrichment",
               tags$p(style = "color:#888; font-style:italic;
                               margin:12px 0 4px; font-size:13px;",
                      "Interactive preview of the selected gene set. All three
                       ontologies (GO_BP, GO_MF, GO_CC) are included in the
                       downloaded Excel file."),
               tags$p(style = "color:#888; font-style:italic;
                               margin:0 0 12px; font-size:13px;",
                      "GO enrichment may take 10\u201320 seconds to compute
                       after changing DEG filters or gene set."),
               plotOutput("enrich_dot", height = "420px"),
               hr(),
               DTOutput("enrich_table")
      )
    )
  })

  output$deg_summary <- renderUI({
    req(raw_data())
    total <- nrow(raw_data())
    n_deg <- nrow(degs())
    n_up  <- sum(degs()$logFC > 0)
    n_dn  <- sum(degs()$logFC < 0)
    tags$div(
      style = "background:#EBF3FB; padding:10px; border-radius:4px;
               margin-bottom:14px; font-size:14px;",
      tags$b(paste0(n_deg, " DEGs")), " from ", total, " genes  |  ",
      tags$span(style = "color:#7B2D8B;", paste0("\u25b2 ", n_up, " up")), "  ",
      tags$span(style = "color:#2E8B57;", paste0("\u25bc ", n_dn, " down")),
      tags$span(style = "color:#888; font-size:12px; margin-left:16px;",
                paste0("adj.P.Val \u2264 ", input$deg_pval,
                       "  |  |log\u2082FC| \u2265 ", input$deg_lfc))
    )
  })

  output$volcano <- renderPlot({
    df  <- raw_data()
    sig <- df$adj.P.Val <= input$deg_pval & abs(df$logFC) >= input$deg_lfc
    df$dir <- dplyr::case_when(!sig ~ "NS", df$logFC > 0 ~ "Up", TRUE ~ "Down")
    top_up   <- df[sig & df$logFC > 0, ] |>
      dplyr::slice_min(adj.P.Val, n = 5, with_ties = FALSE)
    top_down <- df[sig & df$logFC < 0, ] |>
      dplyr::slice_min(adj.P.Val, n = 5, with_ties = FALSE)
    top_lab  <- rbind(top_up, top_down)
    ggplot(df, aes(logFC, -log10(adj.P.Val), colour = dir)) +
      geom_point(alpha = 0.45, size = 1.3) +
      geom_point(data = df[sig, ], alpha = 0.85, size = 2) +
      scale_colour_manual(values = c(NS = "grey70",
                                     Up = "#7B2D8B", Down = "#2E8B57")) +
      geom_vline(xintercept = c(-input$deg_lfc, input$deg_lfc),
                 linetype = "dashed", colour = "black", linewidth = 0.4) +
      geom_hline(yintercept = -log10(input$deg_pval),
                 linetype = "dashed", colour = "black", linewidth = 0.4) +
      ggrepel::geom_text_repel(data = top_lab, aes(label = gene),
                               size = 3, max.overlaps = 15,
                               segment.colour = "grey50") +
      labs(title  = "Volcano Plot",
           x      = expression(log[2]~"Fold Change"),
           y      = expression(-log[10](adj.~P-value)),
           colour = NULL) +
      theme_minimal(base_size = 13) + theme(legend.position = "top")
  })

  output$barplot <- renderPlot({
    validate(need(nrow(degs()) > 0, "No DEGs at current cutoffs."))
    top <- degs() |>
      dplyr::mutate(rank      = sign(logFC) * -log10(adj.P.Val),
                    direction = ifelse(logFC > 0, "Up", "Down")) |>
      dplyr::slice_max(abs(rank), n = 30, with_ties = FALSE)
    ggplot(top, aes(x = reorder(gene, rank), y = rank, fill = direction)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c(Up = "#7B2D8B", Down = "#2E8B57")) +
      coord_flip() +
      labs(title = expression("Top 30 DEGs by rank  " ~
                                (sign(log[2]*FC) %*% -log[10](adj.P.Val))),
           x = NULL, y = "Rank score", fill = NULL) +
      theme_minimal(base_size = 13) + theme(legend.position = "top")
  })

  output$deg_table <- DT::renderDT({
    validate(need(nrow(degs()) > 0, "No DEGs at current cutoffs."))
    cols <- intersect(c("gene","logFC","AveExpr","P.Value","adj.P.Val","B"),
                      names(degs()))
    DT::datatable(
      degs()[, cols],
      options  = list(pageLength = 15,
                      order = list(list(which(cols == "adj.P.Val") - 1, "asc")),
                      scrollX = TRUE),
      rownames = FALSE
    ) |> DT::formatRound(intersect(c("logFC","AveExpr","P.Value","adj.P.Val","B"),
                                   cols), digits = 4)
  })

  output$enrich_dot <- renderPlot({
    ef <- enr_filtered()
    validate(need(nrow(ef) > 0, "No enrichment terms at current cutoff."))
    top <- ef |> dplyr::arrange(p.adjust) |> dplyr::slice_head(n = 15)
    ggplot(top, aes(x = Count, y = reorder(Description, -p.adjust),
                    colour = p.adjust, size = Count)) +
      geom_point() +
      scale_colour_gradient(low = "steelblue", high = "tomato") +
      labs(title  = paste0("Top ", nrow(top), " terms — ", input$geneset,
                           " (adj. P-value \u2264 ", input$enr_pval, ")"),
           x      = "Gene count",
           y      = NULL,
           colour = "adj. P-value") +
      theme_minimal(base_size = 12)
  })

  output$enrich_table <- DT::renderDT({
    ef <- enr_filtered()
    validate(need(nrow(ef) > 0, "No enrichment terms at current cutoff."))
    cols <- intersect(c("Description","GeneRatio","pvalue",
                        "p.adjust","Count","geneID"), names(ef))
    DT::datatable(ef[, cols],
                  options = list(pageLength = 10, scrollX = TRUE),
                  rownames = FALSE) |>
      DT::formatRound(intersect(c("pvalue","p.adjust"), cols), digits = 4)
  })


  # ==========================================================================
  # Rule 6: Download zip folder
  #   1. DEG_results_TIMESTAMP.csv
  #   2. GO_enrichment_TIMESTAMP.xlsx  (GO_BP / GO_MF / GO_CC)
  #   3. DEG_report_TIMESTAMP.html     (figures only)
  # ==========================================================================
  output$dl_report <- downloadHandler(

    filename = function()
      paste0("DEG_analysis_", format(Sys.time(), "%Y%m%d_%H%M"), ".zip"),

    content = function(file) {

      timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
      tmp_dir   <- tempfile()
      dir.create(tmp_dir)

      withProgress(message = "Preparing report folder...", value = 0, {

        # ── 1. DEG CSV ─────────────────────────────────────────
        incProgress(0.05, detail = "Writing DEG results CSV...")
        deg_cols <- intersect(c("gene","logFC","AveExpr","P.Value","adj.P.Val","B"),
                              names(degs()))
        write.csv(
          degs()[order(degs()$adj.P.Val), deg_cols],
          file.path(tmp_dir, paste0("DEG_results_", timestamp, ".csv")),
          row.names = FALSE
        )

        # ── 2. All three GO enrichments ────────────────────────
        gene_symbols <- degs()$gene

        incProgress(0.2, detail = "Running GO Biological Process...")
        enr_bp <- .run_go(gene_symbols, "BP", input$enr_pval)

        incProgress(0.2, detail = "Running GO Molecular Function...")
        enr_mf <- .run_go(gene_symbols, "MF", input$enr_pval)

        incProgress(0.2, detail = "Running GO Cellular Component...")
        enr_cc <- .run_go(gene_symbols, "CC", input$enr_pval)

        # ── 3. GO enrichment Excel ─────────────────────────────
        incProgress(0.1, detail = "Writing GO enrichment Excel...")
        go_cols <- c("Description","GeneRatio","BgRatio",
                     "pvalue","p.adjust","Count","geneID")

        header_style <- createStyle(
          fontName = "Arial", fontSize = 11,
          fontColour = "white", fgFill = "#2E74B5",
          halign = "left", textDecoration = "bold"
        )
        data_style <- createStyle(fontName = "Arial", fontSize = 10)

        wb <- createWorkbook()

        .add_go_sheet <- function(wb, sheet_name, df) {
          if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
          df <- df[order(df$p.adjust), intersect(go_cols, names(df))]
          addWorksheet(wb, sheet_name)
          writeData(wb, sheet_name, df)
          addStyle(wb, sheet_name, header_style,
                   rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
          addStyle(wb, sheet_name, data_style,
                   rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
          setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")
        }

        .add_go_sheet(wb, "GO_BP", enr_bp)
        .add_go_sheet(wb, "GO_MF", enr_mf)
        .add_go_sheet(wb, "GO_CC", enr_cc)

        if (length(sheets(wb)) == 0) {
          addWorksheet(wb, "No_results")
          writeData(wb, "No_results",
                    data.frame(Message = paste("No enriched terms found at",
                                               "adj. P-value \u2264", input$enr_pval)))
        }

        saveWorkbook(
          wb,
          file.path(tmp_dir, paste0("GO_enrichment_", timestamp, ".xlsx")),
          overwrite = TRUE
        )

        # ── 4. HTML report (figures only) ──────────────────────
        incProgress(0.15, detail = "Rendering HTML report...")
        rmarkdown::render(
          input       = "templates/report_template.Rmd",
          output_file = file.path(tmp_dir,
                                  paste0("DEG_report_", timestamp, ".html")),
          params = list(
            data     = raw_data(),
            degs     = degs(),
            deg_pval = input$deg_pval,
            deg_lfc  = input$deg_lfc,
            enr_pval = input$enr_pval,
            enr_bp   = enr_bp,
            enr_mf   = enr_mf,
            enr_cc   = enr_cc
          ),
          envir = new.env(parent = globalenv()),
          quiet = TRUE
        )

        # ── 5. Zip ─────────────────────────────────────────────
        incProgress(0.1, detail = "Compressing...")
        old_wd <- getwd()
        setwd(tmp_dir)
        zip(zipfile = file, files = list.files("."))
        setwd(old_wd)

      }) # end withProgress
    }
  )

} # end server


# ── Launch ───────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
