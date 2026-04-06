# ============================================================
# 00_setup.py
#
# Author:      Patricia Agudelo-Romero
# Email:       patricia.agudeloromero@thekids.org.au
# Institution: The Kids Research Institute Australia
#
# Rule 8 — Lock your dependencies before you share
# Run this ONCE before launching the app for the first time.
#
# Usage:
#   python 00_setup.py
#
# This script:
#   1. Creates (or updates) the conda environment from
#      shiny_app_py_env.yml
#   2. Installs all required pip packages into that environment
#   3. Locks the environment by exporting requirements.txt
#      and updating the yml file
#
# Difference from R (Rule 2):
#   R uses Bioconductor packages (offline, no API call).
#   Python uses gseapy, which queries the Enrichr web API
#   and therefore requires an internet connection for
#   GO enrichment. All other steps work offline.
# ============================================================

import subprocess
import sys
import os
import importlib

# ── Configuration ─────────────────────────────────────────────
ENV_NAME = "shiny_app_py"
YML_FILE = "shiny_app_py_env.yml"

# ── Required packages ────────────────────────────────────────
# (pip_name, import_name, purpose)
PACKAGES = [
    ("shiny",        "shiny",       "reactive web application framework"),
    ("pandas",       "pandas",      "data manipulation (R: dplyr)"),
    ("numpy",        "numpy",       "numerical operations"),
    ("matplotlib",   "matplotlib",  "plotting (R: ggplot2)"),
    ("adjustText",   "adjustText",  "non-overlapping labels (R: ggrepel — approximate)"),
    ("gseapy",       "gseapy",      "GO enrichment via Enrichr API (R: clusterProfiler — INTERNET REQUIRED)"),
    ("jinja2",       "jinja2",      "HTML report templating (R: rmarkdown)"),
    ("openpyxl",     "openpyxl",    "multi-sheet Excel output (R: openxlsx)"),
    ("GEOparse",     "GEOparse",    "download GEO datasets (R: GEOquery)"),
    ("scipy",        "scipy",       "statistical testing for data_prep (R: limma-voom — approximate)"),
    ("statsmodels",  "statsmodels", "multiple testing correction (R: limma BH adjustment)"),
    ("mygene",       "mygene",      "Entrez ID → gene symbol mapping (R: org.Mm.eg.db)"),
]


def run_cmd(cmd, description=""):
    """Run a shell command and return True if successful."""
    print(f"  → {description or ' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ✘  Error:\n{result.stderr.strip()}")
        return False
    return True


def env_exists():
    """Check if the conda environment already exists."""
    result = subprocess.run(
        ["conda", "env", "list"],
        capture_output=True, text=True
    )
    return ENV_NAME in result.stdout


def get_env_python():
    """Get the Python executable path inside the conda environment."""
    result = subprocess.run(
        ["conda", "run", "-n", ENV_NAME, "which", "python"],
        capture_output=True, text=True
    )
    return result.stdout.strip()


# ── Step 1: Create or update conda environment ────────────────
print("=" * 60)
print(f"  Setting up conda environment: {ENV_NAME}")
print("=" * 60)

if not os.path.exists(YML_FILE):
    print(f"\n✘  {YML_FILE} not found in current directory.")
    print(f"   Make sure you are running this script from:")
    print(f"   /mnt/scratch/Projects/App/shiny_python_app/")
    sys.exit(1)

if env_exists():
    print(f"\n  Environment '{ENV_NAME}' already exists — updating...")
    ok = run_cmd(
        ["conda", "env", "update", "-n", ENV_NAME,
         "--file", YML_FILE, "--prune"],
        "Updating environment from yml"
    )
else:
    print(f"\n  Environment '{ENV_NAME}' not found — creating...")
    ok = run_cmd(
        ["conda", "env", "create", "-f", YML_FILE],
        "Creating environment from yml"
    )

if not ok:
    print("\n✘  Failed to create/update conda environment. Exiting.")
    sys.exit(1)

print(f"  ✔  Conda environment '{ENV_NAME}' is ready.\n")


# ── Step 2: Install pip packages into the environment ─────────
print("Checking and installing pip packages...\n")
env_python = get_env_python()

if not env_python:
    print("✘  Could not locate Python in the conda environment.")
    sys.exit(1)

print(f"  Using Python: {env_python}\n")

for pip_name, import_name, purpose in PACKAGES:
    # Check if already installed inside the environment
    check = subprocess.run(
        ["conda", "run", "-n", ENV_NAME,
         "python", "-c", f"import {import_name}"],
        capture_output=True
    )
    if check.returncode == 0:
        print(f"  ✔  {pip_name:<14}  already installed  ({purpose})")
    else:
        print(f"  ↓  {pip_name:<14}  installing...      ({purpose})")
        ok = run_cmd(
            ["conda", "run", "-n", ENV_NAME,
             "pip", "install", "--quiet", pip_name],
            f"pip install {pip_name}"
        )
        if ok:
            print(f"  ✔  {pip_name:<14}  done")
        else:
            print(f"  ✘  {pip_name:<14}  FAILED — check your internet connection")


# ── Step 3: Lock the environment (Rule 8) ─────────────────────
print("\nLocking environment → requirements.txt")
result = subprocess.run(
    ["conda", "run", "-n", ENV_NAME, "pip", "freeze"],
    capture_output=True, text=True
)
with open("requirements.txt", "w") as f:
    f.write(result.stdout)
print("  ✔  requirements.txt written")

print(f"\nUpdating {YML_FILE} with current environment state...")
run_cmd(
    ["conda", "env", "export", "-n", ENV_NAME,
     "--file", YML_FILE],
    f"Exporting environment to {YML_FILE}"
)
print(f"  ✔  {YML_FILE} updated")


# ── Confirm project structure ─────────────────────────────────
print(f"""
✔  Setup complete.

   Conda environment : {ENV_NAME}
   Python            : {env_python}
   Locked files      : requirements.txt  |  {YML_FILE}

   To activate the environment manually:
     conda activate {ENV_NAME}

   Expected project layout:
   shiny_python_app/
   ├── 00_setup.py                  ← this file
   ├── 01_data_prep.py              ← download & prepare GSE60450
   ├── app.py                       ← main Python Shiny application
   ├── 02_deploy.py                 ← deployment options
   ├── shiny_app_py_env.yml         ← conda environment definition
   └── report_template.html         ← Jinja2 HTML report template

   IMPORTANT: GO enrichment requires an internet connection.
   The Enrichr API (used by gseapy) must be reachable.
   This differs from the R implementation, which uses
   org.Mm.eg.db locally (see Rule 2 in the paper).

   Next step: run  python 01_data_prep.py
""")
