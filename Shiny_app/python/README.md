# Limma DEG Explorer — Python Shiny

Parallel implementation of the Limma DEG Explorer in Python Shiny. See the
[repository README](../../README.md) for the overview, input format and the
differences from the R version.

## Run locally

From this folder, with Python 3.11:

```bash
conda env create -f shiny_app_py_env.yml && conda activate shiny_app_py
# or: python -m venv .venv && source .venv/bin/activate && pip install -r requirements-dev.txt
shiny run app.py --reload
```

The example dataset (`data/GSE60450_limma.csv`) is produced by
`../R/01_data_prep.R`, so both apps are validated on the same limma-voom results.

## Files

| File | Purpose |
|---|---|
| `app.py` | User interface and server (reactive pipeline) |
| `analysis.py` | Validation, filtering, GO enrichment and plots |
| `report.py` | Builds the results folder (HTML report, CSV, Excel) |
| `templates/report_template.html` | Jinja2 template for the HTML report |
| `tests/test_app.py` | pytest tests: `pytest -q` |
| `pyproject.toml` | pytest and ruff settings (`ruff check .`, `ruff format .`) |
| `requirements.txt` | Pinned runtime dependencies (used by Docker / Hugging Face) |
| `requirements-dev.txt` | Runtime + test dependencies |
| `shiny_app_py_env.yml` | Conda environment |
| `Dockerfile` | Container image; also used by the Hugging Face Space |

## GO enrichment

gseapy queries the Enrichr web service (GO 2023 libraries), so the server needs
internet access. Every gene in the uploaded table is sent as the background, and
symbols are matched case-insensitively. Because it is a network call, enrichment runs
in a background task when you click **Run GO enrichment**; the rest of the app stays
responsive, and changing the enrichment slider re-filters the cached result instantly.
If Enrichr cannot be reached, the app and the report say so instead of showing an
empty result.
