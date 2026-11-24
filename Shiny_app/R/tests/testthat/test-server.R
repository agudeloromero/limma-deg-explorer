test_that("server: upload, filter, enrichment and download", {
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  a <- load_app()
  a$env$clear_go_cache()
  a$env$.enrich_go <- fake_enrich_go

  tab <- fake_limma_table()
  f <- withr::local_tempfile(fileext = ".csv")
  write.csv(tab, f, row.names = FALSE)

  withr::with_dir(APP_DIR, shiny::testServer(a$app, {
    session$setInputs(deg_pval = 0.05, deg_lfc = 1, enr_pval = 0.05, geneset = "GO_BP",
                      file = data.frame(name = "fake.csv", datapath = f))
    expect_equal(nrow(raw_data()), nrow(tab))
    expect_equal(nrow(degs()), nrow(subset(tab, adj.P.Val <= 0.05 & abs(logFC) >= 1)))

    expect_match(enrichment()$msg, "Run GO enrichment")
    session$setInputs(run_go = 1)
    expect_true(all(enrichment()$res$table$p.adjust <= 0.05))
    session$setInputs(enr_pval = 0.001)
    expect_true(all(enrichment()$res$table$p.adjust <= 0.001))

    a$env$.enrich_go <- function(...) stop("annotation unavailable")
    session$setInputs(geneset = "GO_MF", run_go = 2)
    expect_match(enrichment()$msg, "annotation unavailable")
    session$setInputs(geneset = "GO_BP")
    expect_null(enrichment()$msg)                           # cached BP result still shown

    path <- output$dl_report
    expect_true(file.exists(path))
    expect_length(zip::zip_list(path)$filename, 3)

    bad <- tempfile(fileext = ".csv")
    writeLines("gene,foo\nA,1", bad)
    session$setInputs(file = data.frame(name = "bad.csv", datapath = bad))
    expect_null(raw_data())
    expect_match(loaded()$problems, "Missing required column")
  }))
})
