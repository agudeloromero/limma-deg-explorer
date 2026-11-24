# Limma DEG Explorer — R Shiny

Reference implementation of the Limma DEG Explorer. See the
[repository README](../../README.md) for the overview, input format and the
differences from the Python version.

> **Rule 10:** this app describes the decisions a user makes (thresholds,
> ontology, what to download), not how the code works. Technical notes live
> in the code comments.

## Run locally

From this folder, with R ≥ 4.4:

```bash
Rscript 00_setup.R            # once: install packages and write renv.lock
Rscript 01_data_prep.R        # once: build data/GSE60450_limma.csv from GEO
Rscript -e "shiny::runApp()"
```

Open the RStudio project `shiny_r_app.Rproj` to work on the app interactively.

## Files

| File | Purpose |
|---|---|
| `app.R` | User interface and server (reactive pipeline) |
| `helpers.R` | Validation, filtering, GO enrichment, plots and the results folder |
| `templates/report_template.Rmd` | HTML report included in the download |
| `data/GSE60450_limma.csv` | Reference input, produced by `01_data_prep.R` |
| `tests/` | testthat tests: `Rscript tests/testthat.R` |
| `.lintr` | lintr settings: `Rscript -e "lintr::lint_dir()"` |
| `00_setup.R` | Installs CRAN and Bioconductor packages, writes `renv.lock` |
| `01_data_prep.R` | limma-voom analysis of GSE60450 (writes the CSV for both apps) |
| `02_deploy.R` | Publishes the app to shinyapps.io |
| `Dockerfile` | Container image based on `rocker/shiny` and `renv.lock` |

## GO enrichment

`clusterProfiler::enrichGO()` with the offline `org.Mm.eg.db` annotation. Gene symbols
are mapped to Entrez IDs, every gene in the uploaded table is used as the background,
and P-values are adjusted with Benjamini–Hochberg. The app's enrichment slider is the
only significance filter (clusterProfiler's own cutoffs are set to 1). Results are
cached, so changing the slider re-filters instantly.
