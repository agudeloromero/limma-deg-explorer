# ============================================================
# 00_setup.R  —  Install dependencies and lock them
# Rule 8 — Lock your dependencies before you share
#
# Author:      Patricia Agudelo-Romero
# Email:       patricia.agudeloromero@thekids.org.au
# Institution: The Kids Research Institute Australia
#
# Run once from this folder (Shiny_app/R):
#   Rscript 00_setup.R
#
# 1. Installs the CRAN and Bioconductor packages used by the app,
#    the data preparation script and the tests.
# 2. Writes renv.lock, which records the exact version of every
#    package. Commit renv.lock; others restore it with
#      Rscript -e 'renv::restore()'
# ============================================================

repos <- c(CRAN = "https://cloud.r-project.org")

cran_pkgs <- c(
  "shiny", "DT",               # app framework and interactive tables
  "ggplot2", "ggrepel",        # plots and non-overlapping labels
  "rmarkdown", "knitr",        # HTML report
  "openxlsx", "zip",           # Excel output and zip folder
  "rlang",                     # cache keys for GO results
  "statmod",                   # used by limma
  "testthat", "withr",         # tests
  "renv", "rsconnect",         # dependency locking and deployment
  "BiocManager"
)
bioc_pkgs <- c(
  "limma", "edgeR",            # differential expression (01_data_prep.R)
  "clusterProfiler",           # GO enrichment (offline)
  "org.Mm.eg.db", "AnnotationDbi"   # mouse annotation
)

missing_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing_cran) > 0) install.packages(missing_cran, repos = repos)

missing_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

# renv can only lock a consistent set of Bioconductor packages
bioc_ok <- BiocManager::valid()
if (!isTRUE(bioc_ok)) {
  print(bioc_ok)
  stop("Bioconductor packages come from different releases (see above). ",
       "Update them to one release, e.g. for R 4.5:\n",
       "  Rscript -e 'BiocManager::install(version = \"3.22\", ask = FALSE)'\n",
       "then run 00_setup.R again.", call. = FALSE)
}

# Lock the versions of every package the project's R files use
renv::snapshot(type = "implicit", prompt = FALSE)

message("\nSetup complete. renv.lock written.",
        "\nNext:  Rscript 01_data_prep.R   (builds data/GSE60450_limma.csv)",
        "\nThen:  Rscript -e 'shiny::runApp()'")
