# ============================================================
# 01_data_prep.R  —  Build the reference input for both apps
# Rule 7 — Validate your app against a public, citable dataset
#
# Author:      Patricia Agudelo-Romero
# Email:       patricia.agudeloromero@thekids.org.au
# Institution: The Kids Research Institute Australia
#
# Dataset: GSE60450 (Fu et al. 2015, Nat Cell Biol 17:365-375).
#   Mouse mammary gland, basal and luminal cells at three
#   developmental stages (virgin, 18.5 day pregnancy, 2 day
#   lactation), two biological replicates each (12 samples).
#
# Analysis: limma-voom (Law et al. 2014; Ritchie et al. 2015)
#   counts -> filterByExpr -> TMM -> voom -> lmFit -> eBayes
#   Design: ~ stage + celltype  (cell-type effect adjusted for stage)
#   Contrast reported: basal vs luminal (logFC > 0 = higher in basal)
#
# Output (identical copy written for each app):
#   Shiny_app/R/data/GSE60450_limma.csv
#   Shiny_app/python/data/GSE60450_limma.csv
#
# Run from this folder (Shiny_app/R):
#   Rscript 01_data_prep.R
# ============================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

# ── 1. Download the raw GEO files (skipped if already present) ──
geo_base <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450"
raw_dir  <- file.path("data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

geo_files <- c(
  soft   = "soft/GSE60450_family.soft.gz",
  counts = "suppl/GSE60450_Lactation-GenewiseCounts.txt.gz"
)
paths <- setNames(file.path(raw_dir, basename(geo_files)), names(geo_files))

# Tries R's own downloader first, then the system curl (which follows the
# computer's proxy settings). If both fail, download the two files in a
# browser and save them in data/raw/.
download_geo <- function(url, dest) {
  methods <- c("libcurl", if (nzchar(Sys.which("curl"))) "curl")
  for (m in methods) {
    ok <- tryCatch({
      utils::download.file(url, dest, mode = "wb", method = m, quiet = TRUE)
      file.exists(dest) && file.size(dest) > 0
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) return(invisible(dest))
    unlink(dest)
  }
  stop("Could not download ", url,
       "\nDownload it in a browser and save it as ", normalizePath(dest, mustWork = FALSE),
       ", then run this script again.", call. = FALSE)
}

options(timeout = max(300, getOption("timeout")))
for (k in names(geo_files)) {
  if (!file.exists(paths[[k]])) {
    message("Downloading ", basename(geo_files[[k]]), " from GEO ...")
    download_geo(file.path(geo_base, geo_files[[k]]), paths[[k]])
  }
}

# ── 2. Sample metadata from the GEO SOFT file ───────────────────
# Each sample block holds the sample name (e.g. "MCL1-LA"), the
# immunophenotype and the developmental stage.
read_soft_samples <- function(path) {
  lines  <- readLines(gzfile(path))
  starts <- grep("^\\^SAMPLE = ", lines)
  ends   <- c(starts[-1] - 1, length(lines))

  field <- function(block, key) {
    hits <- grep(paste0("^", key), block, value = TRUE)
    sub(paste0("^", key), "", hits)
  }

  rows <- lapply(seq_along(starts), function(i) {
    block <- lines[starts[i]:ends[i]]
    chars <- field(block, "!Sample_characteristics_ch1 = ")
    descr <- field(block, "!Sample_description = ")
    data.frame(
      geo_accession = sub("^\\^SAMPLE = ", "", block[1]),
      sample        = sub("^Sample name: ", "",
                          grep("^Sample name:", descr, value = TRUE)),
      immunophenotype = sub("^immunophenotype: ", "",
                            grep("^immunophenotype:", chars, value = TRUE)),
      stage_raw     = sub("^developmental stage: ", "",
                          grep("^developmental stage:", chars, value = TRUE))
    )
  })
  do.call(rbind, rows)
}

samples <- read_soft_samples(paths[["soft"]])
samples$celltype <- factor(
  ifelse(grepl("luminal", samples$immunophenotype, ignore.case = TRUE),
         "luminal", "basal"),
  levels = c("luminal", "basal")
)
samples$stage <- factor(
  c("virgin" = "virgin", "18.5 day pregnancy" = "pregnant",
    "2 day lactation" = "lactating")[samples$stage_raw],
  levels = c("virgin", "pregnant", "lactating")
)
stopifnot(!anyNA(samples$celltype), !anyNA(samples$stage))

# ── 3. Count matrix, matched to the metadata BY SAMPLE NAME ─────
# The GEO metadata lists luminal samples first, whereas the count
# file lists basal samples first, so samples must be matched by
# name, never by position.
counts <- read.delim(paths[["counts"]], row.names = 1, check.names = FALSE)
counts <- as.matrix(counts[, colnames(counts) != "Length"])

count_sample <- sub("_.*$", "", colnames(counts))      # "MCL1-DG_BC2..." -> "MCL1-DG"
samples      <- samples[match(count_sample, samples$sample), ]
stopifnot(
  "every count column must match one GEO sample" = !anyNA(samples$sample),
  identical(samples$sample, count_sample)
)
colnames(counts) <- samples$sample
rownames(samples) <- samples$sample

message("Samples (in count-matrix order):")
print(samples[, c("geo_accession", "celltype", "stage")], row.names = TRUE)

# ── 4. limma-voom ───────────────────────────────────────────────
design <- model.matrix(~ stage + celltype, data = samples)

dge  <- DGEList(counts = counts)
keep <- filterByExpr(dge, design)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge)
message(sprintf("Genes kept after filterByExpr: %d of %d", sum(keep), length(keep)))

v   <- voom(dge, design)
fit <- eBayes(lmFit(v, design))
res <- topTable(fit, coef = "celltypebasal", number = Inf,
                adjust.method = "BH", sort.by = "P")

# ── 5. Entrez ID -> gene symbol (offline, org.Mm.eg.db) ─────────
res$gene <- unname(mapIds(org.Mm.eg.db, keys = rownames(res),
                          column = "SYMBOL", keytype = "ENTREZID",
                          multiVals = "first"))
res <- res[!is.na(res$gene) & !duplicated(res$gene), ]

# ── 6. Write the reference CSV for both apps ────────────────────
cols <- c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
out  <- res[, cols]
num  <- vapply(out, is.numeric, logical(1))
out[num] <- lapply(out[num], signif, digits = 6)

targets <- c(file.path("data", "GSE60450_limma.csv"),
             file.path("..", "python", "data", "GSE60450_limma.csv"))
for (f in targets) {
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  write.csv(out, f, row.names = FALSE)
}

top <- out[1, ]
message(sprintf(
  "\nDone: %d genes written to\n  %s\nTop hit: %s (logFC = %.2f, adj.P.Val = %.2e)",
  nrow(out), paste(targets, collapse = "\n  "),
  top$gene, top$logFC, top$adj.P.Val
))
