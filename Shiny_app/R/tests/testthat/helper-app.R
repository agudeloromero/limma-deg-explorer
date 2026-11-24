# Load the app (and helpers.R) into an environment, from the app folder.
APP_DIR <- normalizePath(file.path(testthat::test_path(), "..", ".."))

load_app <- function() {
  env <- new.env(parent = globalenv())
  app <- withCallingHandlers(
    withr::with_dir(APP_DIR, source("app.R", local = env)$value),
    warning = function(w) {
      if (grepl("was built under R version", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
  list(app = app, env = env)
}

fake_limma_table <- function(n = 400, seed = 1) {
  set.seed(seed)
  p <- sort(c(runif(120, 1e-12, 1e-3), runif(n - 120, 1e-3, 1)))
  data.frame(gene = paste0("Gene", seq_len(n)),
             logFC = sample(c(-1, 1), n, TRUE) * runif(n, 1.2, 4),
             AveExpr = runif(n, 0, 10), P.Value = p,
             adj.P.Val = pmin(p * 3, 1), B = rnorm(n))
}

# Offline stand-in for clusterProfiler::enrichGO()
fake_enrich_go <- function(genes, universe, ont) {
  data.frame(ID = sprintf("GO:%07d", 1:25),
             Description = paste(ont, "term", 1:25),
             GeneRatio = paste0(40:16, "/200"), BgRatio = "300/15000",
             pvalue = 10^-seq(12, 0.5, length.out = 25),
             p.adjust = 10^-seq(10, 0.2, length.out = 25),
             qvalue = 10^-seq(10, 0.2, length.out = 25),
             geneID = paste(head(genes, 4), collapse = "/"), Count = 40:16)
}
