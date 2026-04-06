# ============================================================
# 01_data_prep.R
# Rule 7 — Validate your app against a public, citable dataset
#
# Author:      Patricia Agudelo-Romero
# Email:       patricia.agudeloromero@thekids.org.au
# Institution: The Kids Research Institute Australia
#
# Downloads GSE83272 (breast cancer vs normal, GEO) and runs
# a limma differential expression analysis. Produces a CSV
# in the format the Shiny app expects.
#
# Citation:
#   GSE83272 — Breast cancer gene expression profiling.
#   GEO accession GSE83272. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83272
# ============================================================

library(GEOquery)
library(edgeR)
library(limma)

# ── 1. Download expression data from GEO ─────────────────────
message("Downloading GSE60450 from GEO (this may take a minute)...")
dir.create("data", showWarnings = FALSE)
options(timeout = 120)
gse <- getGEO("GSE60450", GSEMatrix = TRUE, destdir = "data/")[[1]]

# ── 2. Extract sample metadata ───────────────────────────────
pheno <- pData(gse)

# ── 3. Download and load raw count matrix ────────────────────
# GSE60450 stores counts as a supplementary file, not in the
# series matrix — exprs(gse) returns an empty matrix.

message("Downloading supplementary count file...")
supp       <- getGEOSuppFiles("GSE60450", makeDirectory = FALSE,
                              baseDir = "data/")
count_file <- rownames(supp)[1]
counts     <- read.delim(count_file, row.names = 1)

# ── 4. Define groups from sample metadata ────────────────────
# GSE60450 contains mouse mammary gland luminal and basal cell
# populations across developmental stages (virgin, pregnancy,
# lactation). Groups are defined from the immunophenotype column.
groups <- factor(
  ifelse(grepl("luminal", pheno$`immunophenotype:ch1`, ignore.case = TRUE),
         "luminal", "basal"),
  levels = c("luminal", "basal")
)
message("Sample groups assigned:")
print(table(groups))

# ── 5. Run limma-voom differential expression ─────────────────
# voom is required for count data — it transforms counts to
# log-CPM with associated precision weights for limma.

# Remove the gene length column before building the count matrix
counts_only <- counts[, -1]   # drop first column (Length)

dge <- edgeR::DGEList(counts = counts_only, group = groups)
dge <- edgeR::calcNormFactors(dge)

design <- model.matrix(~ groups)
v      <- voom(dge, design)
fit    <- lmFit(v, design)
fit    <- eBayes(fit)

# topTable returns all probes, BH-adjusted P-values
res <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")

# ── Map Entrez IDs → gene symbols ────────────────────────────
library(AnnotationDbi)
library(org.Mm.eg.db)

res$gene <- mapIds(
  org.Mm.eg.db,
  keys      = rownames(res),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

# Remove genes with no symbol annotation
res <- res[!is.na(res$gene), ]
#row.names(res) <- res$gene
#res$gene <- NULL

# ── 5. Validate required columns are present ─────────────────
required <- c("gene", "logFC", "AveExpr", "P.Value", "adj.P.Val", "B")
stopifnot(all(required %in% names(res)))

# ── 6. Write CSV for the Shiny app ───────────────────────────
out_path <- file.path("data", "GSE60450_limma.csv")
write.csv(res[, required], out_path, row.names = FALSE)

message(sprintf(
  "Done. Written %d rows to %s\nTop hit: %s (adj.P.Val = %.2e, logFC = %.2f)",
  nrow(res),
  out_path,
  res$gene[1],
  res$adj.P.Val[1],
  res$logFC[1]
))

# ── 7. Validation checklist (manual) ─────────────────────────
# After running this script, open the app and verify:
#   [ ] Upload GSE60450_limma.csv → summary banner shows gene count
#   [ ] Move deg_pval slider → volcano and bar update
#   [ ] Move deg_lfc slider  → volcano and bar update
#   [ ] Switch enrichment filter to "P-value" → table changes
#   [ ] Move enrichment cutoff → dot plot and table update
#   [ ] Click Download Report → HTML opens with correct values
