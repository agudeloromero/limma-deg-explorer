# ============================================================
# app.R  —  Limma DEG Explorer (R Shiny)
#
# Author:      Patricia Agudelo-Romero
# Email:       patricia.agudeloromero@thekids.org.au
# Institution: The Kids Research Institute Australia
#
# Rules illustrated in this file:
#   Rule 2  — framework chosen for Bioconductor ecosystem fit
#   Rule 3  — data logic (helpers.R) is separate from display
#   Rule 4  — uploads are validated at the boundary (read_limma_csv)
#   Rule 5  — every threshold is a UI control
#   Rule 6  — the analysis state is downloadable (build_report_zip)
#   Rule 10 — the UI describes decisions, not implementation details
#
# Run locally (from this folder):   shiny::runApp()
# ============================================================

library(shiny)
library(DT)

source("helpers.R", local = TRUE)

EXAMPLE_PATH    <- file.path("data", "GSE60450_limma.csv")
REPORT_TEMPLATE <- normalizePath(file.path("templates", "report_template.Rmd"))

threshold_input <- function(id, label, min, max, value, step) {
  fluidRow(
    column(8, sliderInput(id, label, min = min, max = max, value = value, step = step)),
    column(4, div(class = "num-box",
                  numericInput(paste0(id, "_num"), NULL, min = min, max = max,
                               value = value, step = step)))
  )
}


# ── UI ───────────────────────────────────────────────────────

ui <- fluidPage(
  title = "Limma DEG Explorer",
  tags$head(tags$style(HTML("
    body { font-family: Arial, sans-serif; }
    .num-box { margin-top: 25px; }
    .num-box input { text-align: center; padding: 4px 6px; }
    h5.section { font-weight: bold; color: #2C3E50; }
    h5.dge { color: #2E74B5; }  h5.enr { color: #1E8449; }
    .hint { font-size: 12px; color: #888; }
    .summary { background: #EBF3FB; padding: 10px; border-radius: 4px;
               margin-bottom: 14px; font-size: 14px; }
    .problem { color: #B00020; font-weight: bold; padding: 8px 0; }
    .go-bar { display: flex; gap: 16px; align-items: center; margin: 12px 0; }
  "))),
  titlePanel("Limma DEG Explorer"),

  sidebarLayout(
    sidebarPanel(width = 3,
      fileInput("file", "Upload limma results (CSV)", accept = ".csv"),
      if (file.exists(EXAMPLE_PATH))
        actionLink("load_example", "or load the example dataset (GSE60450)"),
      hr(),

      h5("DEG filters", class = "section dge"),
      threshold_input("deg_pval", "Adjusted P-value \u2264", 0.001, 1, 0.05, 0.001),
      threshold_input("deg_lfc", "|log\u2082FC| \u2265", 0, 4, 1.0, 0.1),
      hr(),

      h5("Enrichment filters", class = "section enr"),
      selectInput("geneset", "Gene set shown on screen",
                  choices = setNames(names(GO_LABELS), GO_LABELS)),
      threshold_input("enr_pval", "Enrichment adj. P-value \u2264", 0.001, 1, 0.05, 0.001),
      p(class = "hint", "All three GO ontologies are included in the download."),
      hr(),

      h5("Report", class = "section"),
      uiOutput("download_ui")
    ),

    mainPanel(width = 9,
      uiOutput("validation_msg"),
      uiOutput("deg_summary"),
      uiOutput("main_panel")
    )
  )
)


# ── Server ───────────────────────────────────────────────────

server <- function(input, output, session) {

  # Slider <-> numeric box synchronisation (Rule 5)
  for (id in c("deg_pval", "deg_lfc", "enr_pval")) local({
    slider <- id
    num    <- paste0(id, "_num")
    observeEvent(input[[slider]], {
      if (!identical(input[[num]], input[[slider]]))
        updateNumericInput(session, num, value = input[[slider]])
    }, ignoreInit = TRUE)
    observeEvent(input[[num]], {
      req(input[[num]])
      if (!identical(input[[num]], input[[slider]]))
        updateSliderInput(session, slider, value = input[[num]])
    }, ignoreInit = TRUE)
  })

  # ── Data source: uploaded file or example (latest action wins)
  source_file <- reactiveVal(NULL)          # list(path, name)
  observeEvent(input$file, {
    source_file(list(path = input$file$datapath, name = input$file$name))
  })
  observeEvent(input$load_example, {
    source_file(list(path = EXAMPLE_PATH, name = basename(EXAMPLE_PATH)))
  })

  # ── Rule 3: reactive pipeline  loaded -> raw_data -> degs -> enrichment
  loaded <- reactive({
    s <- source_file()
    if (is.null(s)) list(data = NULL, problems = character(0)) else read_limma_csv(s$path)
  })
  raw_data <- reactive(loaded()$data)
  degs <- reactive({
    req(raw_data())
    filter_degs(raw_data(), input$deg_pval, input$deg_lfc)
  })
  go_query <- reactive(list(genes = degs()$gene, universe = raw_data()$gene,
                            ont = GO_ONTOLOGIES[[input$geneset]]))

  # GO enrichment runs on request (it can take tens of seconds for
  # thousands of DEGs); results are cached, so the enrichment P-value
  # slider re-filters instantly.
  go_runs     <- reactiveVal(0)
  go_last_err <- reactiveVal(NULL)          # list(key, msg) of the last failed run
  observeEvent(input$run_go, {
    q <- go_query()
    res <- withProgress(message = "Running GO enrichment\u2026", value = 0.3, {
      run_go(q$genes, q$universe, q$ont, 1)
    })
    go_last_err(if (!is.null(res$error))
      list(key = go_cache_key(q$genes, q$universe, q$ont), msg = res$error))
    go_runs(go_runs() + 1)
  })

  enrichment <- reactive({
    go_runs()
    q <- go_query()
    if (length(unique(q$genes)) < MIN_DEGS_FOR_ENRICHMENT) {
      return(list(res = NULL, msg = sprintf(
        "At least %d DEGs are needed for enrichment.", MIN_DEGS_FOR_ENRICHMENT)))
    }
    if (!is_go_cached(q$genes, q$universe, q$ont)) {
      err <- go_last_err()
      if (!is.null(err) && identical(err$key, go_cache_key(q$genes, q$universe, q$ont))) {
        return(list(res = NULL, msg = err$msg))
      }
      return(list(res = NULL,
                  msg = "Click \u201cRun GO enrichment\u201d to analyse the current DEGs."))
    }
    res <- run_go(q$genes, q$universe, q$ont, input$enr_pval)
    if (!is.null(res$error)) return(list(res = NULL, msg = res$error))
    if (is.null(res$table)) {
      return(list(res = NULL, msg = sprintf("No enriched terms at adj. P-value \u2264 %s.",
                                            input$enr_pval)))
    }
    list(res = res, msg = NULL)
  })

  # ── Messages and summary ─────────────────────────────────
  output$validation_msg <- renderUI({
    problems <- loaded()$problems
    if (length(problems) > 0) div(class = "problem", lapply(problems, p))
  })

  output$deg_summary <- renderUI({
    req(raw_data())
    d <- degs()
    div(class = "summary",
        tags$b(sprintf("%d DEGs", nrow(d))), sprintf(" from %d genes  |  ", nrow(raw_data())),
        span(style = "color:#7B2D8B;", sprintf("\u25b2 %d up", sum(d$logFC > 0))), "  ",
        span(style = "color:#2E8B57;", sprintf("\u25bc %d down", sum(d$logFC < 0))),
        span(style = "color:#888; font-size:12px; margin-left:16px;",
             sprintf("adj.P.Val \u2264 %s  |  |log\u2082FC| \u2265 %s",
                     input$deg_pval, input$deg_lfc)))
  })

  output$main_panel <- renderUI({
    if (is.null(raw_data())) {
      if (length(loaded()$problems) > 0) return(NULL)
      return(div(style = "margin-top:60px; text-align:center; color:#666;",
        h4("Limma DEG Explorer", style = "color:#2E74B5;"),
        p("Upload a limma results table (CSV) using the panel on the left",
          if (file.exists(EXAMPLE_PATH)) ", or load the example dataset." else "."),
        p(style = "font-size:13px;", "Required columns: ",
          lapply(REQUIRED_COLUMNS, function(c) tagList(code(c), " ")),
          "\u2014 the default output of ", code("limma::topTable()"), ".")))
    }
    tabsetPanel(id = "tabs",
      tabPanel("DEG table", DTOutput("deg_table")),
      tabPanel("Volcano plot", plotOutput("volcano", height = "520px")),
      tabPanel("Top ranked genes", plotOutput("barplot", height = "620px")),
      tabPanel("GO enrichment",
        div(class = "go-bar",
            actionButton("run_go", "Run GO enrichment", class = "btn-primary"),
            span(class = "hint", "Runs clusterProfiler (offline) for the selected ontology",
                 "and the current DEGs. The enrichment P-value slider then filters",
                 "the result instantly.")),
        plotOutput("enrich_dot", height = "520px"),
        hr(),
        DTOutput("enrich_table"))
    )
  })

  # ── Plots and tables ─────────────────────────────────────
  output$deg_table <- renderDT({
    tab <- deg_table(degs())
    num <- intersect(c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"), names(tab))
    datatable(tab, rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE)) |>
      formatSignif(num, digits = 4)
  })

  output$volcano <- renderPlot({
    req(raw_data())
    volcano_plot(raw_data(), input$deg_pval, input$deg_lfc)
  })

  output$barplot <- renderPlot({
    d <- degs()
    if (nrow(d) == 0) return(message_plot("No DEGs at the current cutoffs."))
    rank_barplot(d)
  })

  output$enrich_dot <- renderPlot({
    e <- enrichment()
    if (is.null(e$res)) return(message_plot(e$msg))
    go_dotplot(e$res$table, sprintf("%s (adj. P-value \u2264 %s)",
                                    GO_LABELS[[input$geneset]], input$enr_pval))
  })

  output$enrich_table <- renderDT({
    e <- enrichment()
    if (is.null(e$res)) return(datatable(data.frame(Message = e$msg), rownames = FALSE,
                                         options = list(dom = "t")))
    datatable(go_table(e$res$table), rownames = FALSE, filter = "top",
              options = list(pageLength = 10, scrollX = TRUE)) |>
      formatSignif(c("pvalue", "p.adjust"), digits = 4)
  })

  # ── Rule 6: download the analysis folder ─────────────────
  output$download_ui <- renderUI({
    if (is.null(raw_data())) return(p(class = "hint", "Load data to enable the download."))
    tagList(
      p(class = "hint", "A .zip folder with the HTML report, the DEG table (CSV) and",
        "GO enrichment for all three ontologies (Excel)."),
      downloadButton("dl_report", "Download results folder", class = "btn-primary",
                     style = "width:100%;")
    )
  })

  download_time <- NULL
  output$dl_report <- downloadHandler(
    filename = function() {
      download_time <<- Sys.time()             # same timestamp inside the zip
      paste0("DEG_analysis_", format(download_time, "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      d <- degs()
      withProgress(message = "Preparing results folder\u2026", value = 0, {
        go_results <- list()
        for (key in names(GO_ONTOLOGIES)) {
          incProgress(0.25, detail = paste("Running", GO_LABELS[[key]]))
          go_results[[key]] <- run_go(d$gene, raw_data()$gene, GO_ONTOLOGIES[[key]],
                                      input$enr_pval)
        }
        incProgress(0.2, detail = "Rendering report")
        build_report_zip(
          file, raw = raw_data(), degs = d,
          params = list(deg_pval = input$deg_pval, deg_lfc = input$deg_lfc,
                        enr_pval = input$enr_pval, input_name = source_file()$name),
          go_results = go_results, template = REPORT_TEMPLATE,
          now = if (is.null(download_time)) Sys.time() else download_time)
      })
    },
    contentType = "application/zip"
  )
}

shinyApp(ui, server)
