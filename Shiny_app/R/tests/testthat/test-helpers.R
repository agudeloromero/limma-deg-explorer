h <- load_app()$env

test_that("valid and minimal tables pass validation", {
  tab <- fake_limma_table()
  expect_length(h$validate_limma_table(tab), 0)
  expect_length(h$validate_limma_table(tab[, c("gene", "logFC", "adj.P.Val")]), 0)
})

test_that("problems are reported clearly", {
  tab <- fake_limma_table()
  expect_match(h$validate_limma_table(tab[, -2]), "logFC")
  expect_match(h$validate_limma_table(transform(tab, logFC = "x")), "numeric")
  expect_match(h$validate_limma_table(transform(tab, adj.P.Val = 2)), "between 0 and 1")
})

test_that("read_limma_csv drops incomplete rows", {
  tab <- fake_limma_table()
  tab$logFC[1] <- NA
  f <- withr::local_tempfile(fileext = ".csv")
  write.csv(tab, f, row.names = FALSE)
  res <- h$read_limma_csv(f)
  expect_length(res$problems, 0)
  expect_equal(nrow(res$data), nrow(tab) - 1)
})

test_that("filter_degs applies both thresholds", {
  d <- h$filter_degs(fake_limma_table(), 0.05, 1)
  expect_true(all(d$adj.P.Val <= 0.05 & abs(d$logFC) >= 1))
})

test_that("run_go filters, caches and reports errors", {
  h$clear_go_cache()
  calls <- 0
  h$.enrich_go <- function(genes, universe, ont) {
    calls <<- calls + 1
    fake_enrich_go(genes, universe, ont)
  }
  tab <- fake_limma_table()
  res <- h$run_go(tab$gene[1:50], tab$gene, "BP", 0.01)
  expect_null(res$error)
  expect_true(all(res$table$p.adjust <= 0.01))
  h$run_go(tab$gene[1:50], tab$gene, "BP", 0.5)            # cache hit
  expect_equal(calls, 1)

  expect_match(h$run_go(c("A", "B"), tab$gene, "BP", 0.05)$error, "5 DEGs")

  h$clear_go_cache()
  h$.enrich_go <- function(...) stop("offline")
  expect_match(h$run_go(tab$gene[1:50], tab$gene, "BP", 0.05)$error, "offline")
})

test_that("plots build", {
  tab <- fake_limma_table()
  expect_s3_class(h$volcano_plot(tab, 0.05, 1), "ggplot")
  expect_s3_class(h$rank_barplot(h$filter_degs(tab, 0.05, 1)), "ggplot")
  expect_s3_class(h$go_dotplot(fake_enrich_go(tab$gene, NULL, "BP"), "t"), "ggplot")
})

test_that("the report folder has the HTML, CSV and Excel files", {
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  tab <- fake_limma_table()
  d <- h$filter_degs(tab, 0.05, 1)
  go <- list(GO_BP = list(table = fake_enrich_go(d$gene, NULL, "BP"), error = NULL),
             GO_MF = list(table = NULL, error = NULL),
             GO_CC = list(table = NULL, error = "GO enrichment could not be run: offline"))
  zipfile <- withr::local_tempfile(fileext = ".zip")
  h$build_report_zip(zipfile, tab, d,
                     list(deg_pval = 0.05, deg_lfc = 1, enr_pval = 0.05, input_name = "t.csv"),
                     go, template = file.path(APP_DIR, "templates", "report_template.Rmd"))
  files <- zip::zip_list(zipfile)$filename
  expect_length(files, 3)
  expect_true(all(c("html", "csv", "xlsx") %in% tools::file_ext(files)))

  out <- withr::local_tempdir()
  zip::unzip(zipfile, exdir = out)
  html_file <- list.files(out, "html$", full.names = TRUE)
  html <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
  expect_match(html, "t.csv", fixed = TRUE)
  expect_match(html, "offline", fixed = TRUE)
  expect_gte(lengths(regmatches(html, gregexpr("data:image/png", html))), 3)
  expect_equal(openxlsx::getSheetNames(list.files(out, "xlsx$", full.names = TRUE)), "GO_BP")
  csv <- read.csv(list.files(out, "csv$", full.names = TRUE))
  expect_equal(nrow(csv), nrow(d))
})

test_that("real GO enrichment finds keratin terms (needs org.Mm.eg.db)", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Mm.eg.db")
  h <- load_app()$env                      # fresh copy, without the stand-ins above
  keratins <- c("Krt5", "Krt14", "Krt17", "Krt6a", "Krt15", "Krt16",
                "Krt7", "Krt8", "Krt18", "Krt19")
  universe <- unique(c(keratins, head(AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db,
                                                          "SYMBOL"), 5000)))
  res <- h$run_go(keratins, universe, "CC", 0.05)
  expect_null(res$error)
  expect_true(any(grepl("keratin|intermediate filament", res$table$Description)))
})
