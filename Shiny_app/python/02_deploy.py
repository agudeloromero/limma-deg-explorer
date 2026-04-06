# ============================================================
# 02_deploy.py
# Rule 9 — Deploy to a URL, not an inbox
#
# Author:      Patricia Agudelo-Romero
# Email:       patricia.agudeloromero@thekids.org.au
# Institution: The Kids Research Institute Australia
#
# Three deployment options in order of recommended use:
#   A. Hugging Face Spaces — free, permanent, Git-versioned,
#                            citable URL (RECOMMENDED)
#   B. Posit Connect Cloud — rsconnect-python, simplest managed option
#   C. Docker              — maximum control, any cloud provider
#
# Difference from R (Rule 9):
#   R → ShinyApps.io is the natural target.
#   Python → Hugging Face Spaces is the best free option:
#            it is Git-versioned, permanently citable, and
#            free with no usage limits for public apps.
# ============================================================

import os
import subprocess
import sys

# ── Option A: Hugging Face Spaces (recommended) ──────────────
#
# 1. Create a free account at https://huggingface.co
# 2. Go to Spaces → New Space → Docker (not Gradio or Streamlit)
# 3. Set Space name, e.g.  limma-deg-explorer
# 4. Run the commands below from the project root
#
# The Dockerfile below is written automatically by this script.

HF_USERNAME   = "YOUR_HF_USERNAME"       # ← replace
HF_SPACE_NAME = "limma-deg-explorer"     # ← rename if desired
HF_SPACE_URL  = f"https://huggingface.co/spaces/{HF_USERNAME}/{HF_SPACE_NAME}"

dockerfile_hf = f"""\
FROM python:3.11-slim

WORKDIR /app

# System dependencies
RUN apt-get update && apt-get install -y \\
    libxml2-dev libcurl4-openssl-dev \\
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application files
COPY . .

# Expose Shiny default port
EXPOSE 7860

CMD ["shiny", "run", "app.py", "--host", "0.0.0.0", "--port", "7860"]
"""

# ── Option B: Posit Connect Cloud ─────────────────────────────
#
# pip install rsconnect-python
# rsconnect deploy shiny . --name {HF_SPACE_NAME}

posit_deploy_cmd = [
    "rsconnect", "deploy", "shiny", ".",
    "--name", HF_SPACE_NAME,
]

# ── Option C: Docker (any cloud) ─────────────────────────────
#
# Uses python:3.11-slim — much smaller than rocker/shiny (~200 MB vs ~800 MB)
# This is a key advantage of the Python implementation (Rule 2).

dockerfile_docker = """\
FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y \\
    libxml2-dev libcurl4-openssl-dev \\
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8000

CMD ["shiny", "run", "app.py", "--host", "0.0.0.0", "--port", "8000"]
"""


def write_dockerfile(content, path="Dockerfile"):
    with open(path, "w") as f:
        f.write(content)
    print(f"Dockerfile written to {path}")


def deploy_huggingface():
    """Deploy to Hugging Face Spaces via git push."""
    print("\n── Option A: Hugging Face Spaces ─────────────────────")
    write_dockerfile(dockerfile_hf)

    remote_url = f"https://huggingface.co/spaces/{HF_USERNAME}/{HF_SPACE_NAME}"
    cmds = [
        ["git", "init"],
        ["git", "add", "."],
        ["git", "commit", "-m", "Initial deployment"],
        ["git", "remote", "add", "hf", remote_url],
        ["git", "push", "hf", "main"],
    ]
    for cmd in cmds:
        print(f"  Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  Error: {result.stderr}")
            return
    print(f"\n  ✔  Deployed to: {HF_SPACE_URL}")


def deploy_posit():
    """Deploy to Posit Connect Cloud."""
    print("\n── Option B: Posit Connect Cloud ─────────────────────")
    print("  Running: " + " ".join(posit_deploy_cmd))
    subprocess.run(posit_deploy_cmd)


def deploy_docker():
    """Build and run Docker container locally."""
    print("\n── Option C: Docker ───────────────────────────────────")
    write_dockerfile(dockerfile_docker)
    subprocess.run(["docker", "build", "-t", "limma-shiny-py", "."])
    print("  To run locally:")
    print("    docker run -p 8000:8000 limma-shiny-py")
    print("    Open: http://localhost:8000")


# ── Run ───────────────────────────────────────────────────────
if __name__ == "__main__":
    print("""
Limma DEG Explorer — Python Shiny Deployment
=============================================

Choose a deployment option:
  A  Hugging Face Spaces (free, permanent, citable URL)
  B  Posit Connect Cloud (managed, simplest)
  C  Docker (local test or any cloud)
""")
    choice = input("Enter A, B, or C: ").strip().upper()

    if choice == "A":
        deploy_huggingface()
    elif choice == "B":
        deploy_posit()
    elif choice == "C":
        deploy_docker()
    else:
        print("Invalid choice. Please enter A, B, or C.")

# ── Post-deployment checklist ─────────────────────────────────
# After deploying, verify the live URL:
#   [ ] Upload GSE60450_limma.csv → summary banner appears
#   [ ] DEG sliders respond without timeout
#   [ ] GO Enrichment tab loads (requires internet from server)
#   [ ] Download Results Folder → ZIP with HTML, CSV, Excel
#   [ ] Record the URL in the paper's Data Availability statement
