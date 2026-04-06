# ============================================================
# 00_setup.R
#
# Author:      Patricia Agudelo-Romero
# Email:       patricia.agudeloromero@thekids.org.au
# Institution: The Kids Research Institute Australia
#
# Rule 8 — Lock your dependencies before you share
# Run this ONCE before launching the app for the first time.
# ============================================================

# ── 1. CRAN packages ─────────────────────────────────────────
cran_pkgs <- c(
  "shiny",       # reactive web framework
  "ggplot2",     # plotting
  "dplyr",       # data manipulation
  "statmod",     # required by limma's eBayes()
  "DT",          # interactive tables
  "rmarkdown",   # report generation
  "knitr",       # report engine
  "ggrepel",     # non-overlapping volcano labels
  "openxlsx",   # multi-sheet Excel output for GO enrichment
  "renv"         # dependency locking (Rule 8)
)

install.packages(cran_pkgs, repos = "https://cloud.r-project.org")

# ── 2. Bioconductor packages ─────────────────────────────────
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",   # KEGG / GO enrichment
  "org.Mm.eg.db",      # mouse gene annotation (offline, no API call)
  "AnnotationDbi",     # mapIds() for Entrez ID → gene symbol conversion
  "enrichplot",        # enrichment visualisation
  "edgeR",             # count normalisation and voom weights
  "limma",             # DEG analysis (for data_prep script)
  "GEOquery"           # download GSE60450 (Rule 7 validation)
))

# ── 3. Lock environment with renv (Rule 8) ───────────────────
# Run after installing all packages above.
# renv::init()        # initialise renv in this project folder
# renv::snapshot()    # write renv.lock

# ── 4. Confirm project structure ─────────────────────────────
# Expected layout:
#   limma_app/
#   ├── 00_setup.R              <- this file
#   ├── 01_data_prep.R          <- download & prepare GSE83272
#   ├── app.R                   <- main Shiny application
#   ├── templates/
#   │   └── report_template.Rmd <- parameterised HTML report
#   ├── data/
#   │   └── GSE83272_limma.csv  <- produced by 01_data_prep.R
#   └── renv.lock               <- produced by renv::snapshot()

message("Setup complete. Run 01_data_prep.R next to prepare the test dataset.")
