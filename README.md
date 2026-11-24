# limma-deg-explorer

[![CI](https://github.com/agudeloromero/limma-deg-explorer/actions/workflows/ci.yml/badge.svg)](https://github.com/agudeloromero/limma-deg-explorer/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/agudeloromero/limma-deg-explorer)](https://github.com/agudeloromero/limma-deg-explorer/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Parallel **R Shiny** and **Python Shiny** applications for exploring limma differential
expression results and Gene Ontology (GO) enrichment. Companion code for the manuscript:

> **Golden Rules for Turning Bioinformatics Analyses into Interactive Web Applications**
> Agudelo-Romero P, Caparros-Martin JA, Bates A, Sikazwe C, Blyth C.
> *Journal, year, DOI — to be added on acceptance.*

| | R Shiny (reference implementation) | Python Shiny (parallel implementation) |
|---|---|---|
| Live app | [shinyapps.io](https://agudeloromero.shinyapps.io/limma-deg-explorer) | [Hugging Face Spaces](https://huggingface.co/spaces/agudeloromero/limma-deg-explorer) |
| Source | [`Shiny_app/R`](Shiny_app/R) | [`Shiny_app/python`](Shiny_app/python) |

Companion website with implementation notes: <https://agudeloromero.github.io/limma-deg-explorer>

---

## What the apps do

1. **Load** a limma results table (CSV), or the bundled example dataset.
2. **Filter** DEGs with adjusted P-value and |log₂FC| thresholds (slider or typed value).
3. **Explore** a DEG table, a volcano plot and a ranked bar plot that update as the thresholds change.
4. **Run GO enrichment** for the current DEGs (Biological Process, Molecular Function or
   Cellular Component), using all genes in the table as the background, then filter terms
   by adjusted P-value.
5. **Download** a results folder (`.zip`) with a self-contained HTML report, the DEG table (CSV)
   and GO enrichment for all three ontologies (Excel), recording every parameter used.

## Input format

A comma-separated file with, at minimum:

| Column | Description |
|---|---|
| `gene` | Gene symbol (mouse symbols for the GO step, e.g. `Krt5`) |
| `logFC` | log₂ fold change |
| `adj.P.Val` | Benjamini–Hochberg adjusted P-value |

`AveExpr`, `t`, `P.Value` and `B` are shown when present. This is the output of
`limma::topTable()` with a `gene` column added. Problems with an upload are reported
straight away with a message naming what is wrong.

## Repository structure

```
limma-deg-explorer/
├── Shiny_app/
│   ├── R/                         R Shiny app (reference implementation)
│   │   ├── app.R                  user interface and server
│   │   ├── helpers.R              data logic, GO enrichment, plots, report
│   │   ├── templates/             R Markdown report template
│   │   ├── data/                  reference input (GSE60450_limma.csv)
│   │   ├── tests/                 testthat tests
│   │   ├── 00_setup.R             install packages, write renv.lock
│   │   ├── 01_data_prep.R         build the reference input from GEO
│   │   ├── 02_deploy.R            publish to shinyapps.io
│   │   └── Dockerfile
│   └── python/                    Python Shiny app (parallel implementation)
│       ├── app.py                 user interface and server
│       ├── analysis.py            data logic, GO enrichment, plots
│       ├── report.py              results folder and HTML report
│       ├── templates/             Jinja2 report template
│       ├── data/                  reference input (same file as R)
│       ├── tests/                 pytest tests
│       ├── requirements.txt       pinned runtime dependencies
│       ├── requirements-dev.txt   + test dependencies
│       ├── shiny_app_py_env.yml   conda environment
│       └── Dockerfile             also used by the Hugging Face Space
├── .github/
│   ├── workflows/                 CI, release, deployments, website
│   ├── scripts/check_version.sh   version consistency check
│   └── dependabot.yml             monthly dependency updates
├── *.qmd, _quarto.yml, assets/    companion website (Quarto)
├── CHANGELOG.md
├── CITATION.cff
└── LICENSE
```

## Quick start

**R** (R ≥ 4.4), from `Shiny_app/R`:

```bash
Rscript 00_setup.R             # install packages and write renv.lock (once)
Rscript 01_data_prep.R         # build data/GSE60450_limma.csv (once, ~2 min)
Rscript -e "shiny::runApp()"
```

To reproduce the exact environment instead, run `Rscript -e "renv::restore()"`.

**Python** (Python 3.11), from `Shiny_app/python`:

```bash
conda env create -f shiny_app_py_env.yml && conda activate shiny_app_py
# or: python -m venv .venv && source .venv/bin/activate && pip install -r requirements-dev.txt
shiny run app.py
```

The Python app reads the same reference CSV, which `01_data_prep.R` writes into both apps.

## Reference dataset (GSE60450)

The apps are validated on GSE60450 (Fu *et al.* 2015): mouse mammary basal and luminal
cells at three developmental stages, two replicates each. `01_data_prep.R` downloads the
counts and sample annotation from GEO, matches samples **by name**, keeps expressed genes
with `filterByExpr`, applies TMM normalisation and fits limma-voom with
`~ stage + celltype`. The reported contrast is basal vs luminal, so **logFC > 0 means
higher in basal cells** (for example *Krt5*, *Trp63*, *Acta2*) and logFC < 0 means
higher in luminal cells (for example *Csn2*, *Wap*, *Elf5*).

## Tests and code quality

```bash
(cd Shiny_app/R      && Rscript tests/testthat.R && Rscript -e 'lintr::lint_dir()')
(cd Shiny_app/python && pytest -q && ruff check . && ruff format --check .)
```

The tests cover input validation, filtering, GO enrichment (with an offline stand-in for
the enrichment service, plus a real clusterProfiler run where it is installed), the full
server logic and the downloaded results folder.

## Continuous integration and delivery

| Workflow | Runs on | What it does |
|---|---|---|
| [CI](.github/workflows/ci.yml) | every push and pull request | version check; lint (lintr, ruff) and tests for both apps |
| [Release](.github/workflows/release.yml) | a version tag (`v1.2.3`) | CI → Docker images on GitHub Container Registry → GitHub Release → deployment of both apps |
| [Deploy Python app](.github/workflows/deploy-hf-space.yml) | each release, or by hand | publishes `Shiny_app/python` to Hugging Face Spaces |
| [Deploy R app](.github/workflows/deploy-shinyapps.yml) | each release, or by hand | publishes `Shiny_app/R` to shinyapps.io with the `renv.lock` versions |
| [Website](.github/workflows/publish.yml) | every push to `main` | renders the Quarto site to GitHub Pages |

Dependabot opens monthly pull requests for GitHub Actions and Python dependencies; each
one is merged only if CI passes.

## Versioning and releases

Both apps share one [semantic version](https://semver.org/), recorded in
`Shiny_app/R/helpers.R`, `Shiny_app/python/report.py`, `CITATION.cff` and `CHANGELOG.md`
(CI fails if they disagree). Each app shows its version in the downloaded report.

To publish a release:

1. Update the version in those four files and describe the changes in `CHANGELOG.md`.
2. Commit and push.
3. Tag the commit (for example `v1.1.0`) and push the tag. In GitHub Desktop:
   **History → right-click the commit → Create Tag**, then **Push origin**.

The release workflow then checks the tag, runs CI, publishes the Docker images, creates
the GitHub Release from the changelog entry and redeploys both apps. With the repository
enabled on [Zenodo](https://zenodo.org/account/settings/github/), every release is
archived with its own DOI.

Released images:

```bash
docker run --rm -p 7860:7860 ghcr.io/agudeloromero/limma-deg-explorer-py:1.0.0
docker run --rm -p 3838:3838 ghcr.io/agudeloromero/limma-deg-explorer-r:1.0.0
```

## R and Python: what differs

| | R | Python |
|---|---|---|
| Plots | ggplot2 + ggrepel | matplotlib + adjustText |
| GO enrichment | clusterProfiler + org.Mm.eg.db, **offline** | gseapy → Enrichr web service, **internet required** |
| Background gene set | all genes in the uploaded table | all genes in the uploaded table |
| Report | R Markdown | Jinja2 template, figures embedded as PNG |
| Dependency lock | `renv.lock` | `requirements.txt` (pinned) |
| Hosting | shinyapps.io | Hugging Face Spaces (Docker) |

GO term lists differ between the two because the annotation sources differ
(org.Mm.eg.db versus Enrichr's GO 2023 libraries). The R app is the reference
implementation for the analysis described in the manuscript.

## Deployment

Releases deploy both apps automatically once these repository secrets exist
(**Settings → Secrets and variables → Actions**); without them the deployment step is
skipped with a warning:

| Secret | Where to get it |
|---|---|
| `HF_TOKEN` | huggingface.co → Settings → Access Tokens (write) |
| `SHINYAPPS_ACCOUNT`, `SHINYAPPS_TOKEN`, `SHINYAPPS_SECRET` | shinyapps.io → Account → Tokens |

To deploy by hand instead: `Rscript 02_deploy.R` from `Shiny_app/R` (with the same three
values in `~/.Renviron`), or run a deploy workflow from the Actions tab. Each app folder
also has a `Dockerfile`.

## Troubleshooting

- **The upload is rejected:** the message names the missing or invalid column; column
  names are case-sensitive.
- **Python GO enrichment fails:** the Enrichr web service must be reachable from the
  machine running the app. The rest of the app, and the download, still work.
- **R GO enrichment fails:** check that `clusterProfiler` and `org.Mm.eg.db` are installed
  (`Rscript 00_setup.R`) and that the gene symbols are mouse symbols.

## References

1. Agudelo-Romero P, Caparros-Martin JA, Bates A, Sikazwe C, Blyth C. Golden Rules for Turning Bioinformatics Analyses into Interactive Web Applications. *Journal, year, DOI to be added.*
2. Chang W, Cheng J, Allaire J, *et al.* shiny: Web Application Framework for R. https://CRAN.R-project.org/package=shiny
3. Posit Software, PBC. shiny: Web Application Framework for Python. https://pypi.org/project/shiny/
4. Ritchie ME, Phipson B, Wu D, *et al.* limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res.* 2015;43:e47. doi:10.1093/nar/gkv007
5. Law CW, Chen Y, Shi W, Smyth GK. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biol.* 2014;15:R29. doi:10.1186/gb-2014-15-2-r29
6. Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics.* 2010;26:139–140. doi:10.1093/bioinformatics/btp616
7. Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS.* 2012;16:284–287. doi:10.1089/omi.2011.0118
8. Fang Z, Liu X, Peltz G. GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics.* 2023;39:btac757. doi:10.1093/bioinformatics/btac757
9. Kuleshov MV, Jones MR, Rouillard AD, *et al.* Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. *Nucleic Acids Res.* 2016;44:W90–W97. doi:10.1093/nar/gkw377
10. Fu NY, Rios AC, Pal B, *et al.* EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival. *Nat Cell Biol.* 2015;17:365–375. doi:10.1038/ncb3117
11. Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. *Nucleic Acids Res.* 2002;30:207–210. doi:10.1093/nar/30.1.207

## Licence and contact

Code released under the [MIT License](LICENSE). Citation metadata: [`CITATION.cff`](CITATION.cff).

Patricia Agudelo-Romero — patricia.agudeloromero@thekids.org.au
The Kids Research Institute Australia, Perth, Western Australia
