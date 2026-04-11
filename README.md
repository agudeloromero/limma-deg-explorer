# limma-deg-explorer

Parallel R Shiny [2] and Python Shiny [3] applications for interactive differential gene expression and Gene Ontology enrichment analysis. Companion code for the manuscript:

> **Ten Simple Rules for Turning Bioinformatics Analyses into Interactive Web Applications** [1]  
> Agudelo-Romero P, Caparros-Martin JA, Bates A, Sikazwe C, Blyth C.  
> [journal, year, doi] — to be updated on acceptance

---

## Live applications

| Implementation | URL |
|---|---|
| R Shiny [2] | [https://agudeloromero.shinyapps.io/limma-deg-explorer](https://agudeloromero.shinyapps.io/limma-deg-explorer) |
| Python Shiny [3] | [https://huggingface.co/spaces/agudeloromero/limma-deg-explorer](https://huggingface.co/spaces/agudeloromero/limma-deg-explorer) |

Supporting website with expanded implementation notes and reproducibility instructions:  
[https://agudeloromero.github.io/limma-deg-explorer](https://agudeloromero.github.io/limma-deg-explorer)

---

## Repository structure

```
limma-deg-explorer/
├── Shiny_app/
│   ├── R/                  # R Shiny application (reference implementation)
│   │   ├── app.R
│   │   ├── renv.lock       # Locked R environment
│   │   └── 00_setup.R
│   │   └── 01_data_prep.R
│   └── Python/             # Python Shiny application (parallel implementation)
│       ├── app.py
│       ├── shiny_app_py_env.yml   # Conda environment definition
│       ├── requirements.txt       # pip dependencies
│       └── 00_setup.py
│       └── 01_data_prep.py
├── assets/                 # Figures and supporting assets
├── index.qmd               # Companion website source (Quarto)
├── methods.qmd
├── reproducibility.qmd
├── supplementary-code-examples.qmd
└── availability.qmd
```

---

## What the applications do

Both applications accept a CSV file produced by a limma [4] differential expression analysis and provide:

- Interactive adjustment of adjusted P-value and log2 fold-change thresholds
- Volcano plot and bar plot updated reactively to threshold changes
- Gene Ontology enrichment analysis with interactive filter controls
- Live summary banner showing total genes, DEGs, and up/down-regulated subsets
- Downloadable HTML report capturing the current analysis state

The R implementation uses `limma` [4], `clusterProfiler` [5], and R Markdown for reporting. The Python implementation uses `PyDESeq2` (a Python reimplementation of DESeq2) [6,7] and `gseapy` [8] for enrichment, and Jinja2 for reporting. Both are validated against the public dataset GSE60450 [9].

---

## Input format

Both applications accept a **comma-separated CSV file** with at minimum the following columns:

| Column | Description |
|---|---|
| `gene` | Gene identifier |
| `logFC` | Log2 fold-change |
| `adj.P.Val` | Adjusted P-value |

Missing or misnamed columns are caught at upload with an informative error message.

---

## Restore the R environment

The R implementation uses `renv` for dependency locking.

```r
# From the R/ directory
Rscript -e "renv::restore()"
```

---

## Restore the Python environment

```bash
# Create and activate the conda environment
conda env create -f Python/shiny_app_py_env.yml
conda activate shiny_app_py

# Install pip dependencies
pip install -r Python/requirements.txt
```

---

## Prepare the validation dataset

The manuscript uses GSE60450 [9], a publicly available limma tutorial dataset from GEO [10]. Data preparation completes in under five minutes on a personal computer.

**R workflow:**

```bash
Rscript R/00_setup.R
Rscript R/01_data_prep.R
```

**Python workflow:**

```bash
python Python/00_setup.py
conda activate shiny_app_py
python Python/01_data_prep.py
```

---

## Run the applications locally

**R Shiny [2]:**

```r
Rscript -e "shiny::runApp('R')"
```

**Python Shiny [3]:**

```bash
shiny run Python/app.py
```

---

## Troubleshooting

**Uploaded file is rejected**  
Check that the CSV contains the required columns (`gene`, `logFC`, `adj.P.Val`) with exact spelling. The application validates columns at upload and returns an informative message if any are missing.

**Python enrichment does not run**  
The Python implementation uses the Enrichr API via `gseapy` [8], which requires an internet connection. Check network access before running.

**Environment restores but application still fails**  
Recreate the environment from the lock files rather than reusing a modified environment. For R, run `renv::restore()` in a fresh session. For Python, delete and recreate the conda environment.

**Results differ between R and Python**  
The two implementations illustrate the same application-design principles but use different analytical libraries (`clusterProfiler` [5] vs `gseapy` [8]). The R implementation is the reference implementation for the genomics workflow described in the manuscript.

---

## References

1. Agudelo-Romero P, Caparros-Martin JA, Bates A, Sikazwe C, Blyth C. Ten Simple Rules for Turning Bioinformatics Analyses into Interactive Web Applications. [journal, year, doi] — to be updated on acceptance.

2. Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Xie Y, Allen J, McPherson J, Dipert A, Borges B. shiny: Web Application Framework for R. R package version 1.8.0. 2023. https://CRAN.R-project.org/package=shiny

3. Posit Software, PBC. shiny: Web Application Framework for Python. Python package version 1.0.0. 2024. https://pypi.org/project/shiny/

4. Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43:e47. doi:10.1093/nar/gkv007

5. Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012;16:284–287. doi:10.1089/omi.2011.0118

6. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15. doi:10.1186/s13059-014-0550-8

7. Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. Bioinformatics. 2023;39. doi:10.1093/bioinformatics/btad547

8. Fang Z, Liu X, Peltz G. GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. Bioinformatics. 2023;39. doi:10.1093/bioinformatics/btac757

9. Fu NY, Rios AC, Pal B, Soetanto R, Lun ATL, Liu K, et al. EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival. Nat Cell Biol. 2015;17:365–375. doi:10.1038/ncb3117

10. Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. Nucleic Acids Res. 2002;30:207–210. doi:10.1093/nar/30.1.207

---

## License

This repository is made available under the [MIT License](LICENSE).

---

## Contact

Patricia Agudelo-Romero — patricia.agudeloromero@thekids.org.au  
The Kids Research Institute Australia, Perth, Western Australia
