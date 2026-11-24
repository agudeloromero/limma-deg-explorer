# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/). Both applications share one version
number, so a release always contains a matching R and Python app.

## [Unreleased]

## [1.0.0] - 2026-11-24

First release accompanying the manuscript *Golden Rules for Turning Bioinformatics
Analyses into Interactive Web Applications*.

### Added
- Downloadable results folder in both apps: self-contained HTML report, DEG table (CSV)
  and GO enrichment for three ontologies (Excel), recording every parameter used.
- Example dataset (GSE60450 limma-voom results) loadable from the app.
- GO enrichment uses every gene in the uploaded table as the background, runs on request
  and is cached, so the enrichment P-value slider re-filters instantly.
- Test suites (testthat, pytest), linting (lintr, ruff), and continuous integration.
- Release pipeline: versioned Docker images on GitHub Container Registry, GitHub
  Releases, and deployment to Hugging Face Spaces and shinyapps.io.
- Dockerfiles, `renv.lock` (R) and pinned `requirements.txt` (Python).
- `LICENSE` (MIT), `CITATION.cff` and this changelog.

### Changed
- Reference input rebuilt with limma-voom: samples matched by name, `filterByExpr`,
  TMM normalisation and the design `~ stage + celltype` (logFC > 0 = higher in basal
  cells). Both apps now use the same file.
- Only `gene`, `logFC` and `adj.P.Val` are required in an upload; other limma columns
  are shown when present. Upload problems are reported with a specific message.
- Data logic separated from the user interface (`helpers.R`; `analysis.py`, `report.py`).
- Python GO enrichment runs in a background task, so a slow Enrichr response no longer
  blocks other users; failures are reported in the app and the report.

### Fixed
- Sample labels in the R data preparation were assigned by position, which reversed the
  sign of every log fold change.
- Report download failed in both apps (missing R Markdown template; wrong template path
  in Python).
- clusterProfiler's default q-value cutoff (0.2) silently filtered GO terms in addition
  to the app's slider.

### Removed
- Python t-test approximation of limma (`01_data_prep.py`), Python setup and deployment
  helper scripts, committed raw GEO files and compiled Python files.

[Unreleased]: https://github.com/agudeloromero/limma-deg-explorer/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/agudeloromero/limma-deg-explorer/releases/tag/v1.0.0
