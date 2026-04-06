**Author:** Patricia Agudelo-Romero  
**Institution:** The Kids Research Institute Australia  
**Contact:** patricia.agudeloromero@thekids.org.au

---
# Limma DEG Interactive Explorer — R Shiny
---


> **Rule 10:** This README describes what the app helps you *decide*,
> not how the code works. Technical documentation belongs in comments.

---

## What this app does

You have a limma differential expression result. Your collaborator wants
to know which genes matter — but they can't run R. This app turns your
CSV into a point-and-click explorer they can use from any browser, with
no installation required.

**The user journey:**

1. **Upload** your limma output (CSV with gene, logFC, adj.P.Val)
2. **Filter DEGs** — drag the adjusted P-value and |log2FC| sliders
3. **Explore** — volcano plot and bar plot update live as you filter
4. **Enrich** — run KEGG, GO-BP, or GO-MF enrichment on the filtered gene list
5. **Filter enrichment** — choose P-value or adj. P-value cutoff independently
6. **Download** — get a self-contained HTML report capturing all settings and results

---

## Required input format

Your CSV must contain these columns (case-sensitive):

| Column | Description |
|---|---|
| `gene` | HGNC gene symbols (e.g. TP53, BRCA1) |
| `logFC` | log2 fold change |
| `AveExpr` | average expression across samples |
| `P.Value` | unadjusted P-value |
| `adj.P.Val` | Benjamini–Hochberg adjusted P-value |

This is the default output of `limma::topTable()`. No reformatting needed.

---

## Quick start (local)

```r
# 1. Install dependencies (once)
source("00_setup.R")

# 2. Prepare the validation dataset (optional but recommended)
source("01_data_prep.R")

# 3. Launch the app
shiny::runApp(".")
```

---

## Reproducing the paper results

The app was validated against GSE83272 (breast cancer vs. normal, GEO).
Run `01_data_prep.R` to download and process this dataset automatically.
The output CSV (`data/GSE83272_limma.csv`) is the file used in the paper.

**Default filter settings used in the paper:**
- adj.P.Val ≤ 0.05
- |log2FC| ≥ 1.0
- Gene set: KEGG (Human)
- Enrichment filter: adj. P-value ≤ 0.05

---

## Deployment

```r
source("02_deploy.R")   # follow prompts for ShinyApps.io or Docker
```

Live app: **[URL to be added after deployment]**

---

## File structure

```
limma_app/
├── 00_setup.R              — install dependencies
├── 01_data_prep.R          — download and prepare GSE83272
├── app.R                   — Shiny application (Rules 2–6, 9–10)
├── 02_deploy.R             — deployment to ShinyApps.io / Docker
├── templates/
│   └── report_template.Rmd — parameterised HTML report (Rule 6)
├── data/
│   └── GSE83272_limma.csv  — test dataset (produced by 01_data_prep.R)
└── renv.lock               — locked dependency versions (Rule 8)
```

---

## Citation

If you use this application or the accompanying paper in your work,
please cite:

> [Author names]. Ten Simple Rules for Deploying Bioinformatics Analyses
> as Interactive Web Applications, with a Comparison of R and Python Shiny.
> PLOS Computational Biology. [Year]. doi:[to be assigned]
