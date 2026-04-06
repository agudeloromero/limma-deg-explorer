# ============================================================
# 02_deploy.R
# Rule 9 — Deploy to a URL, not an inbox
#
# Author:      Patricia Agudelo-Romero
# Email:       patricia.agudeloromero@thekids.org.au
# Institution: The Kids Research Institute Australia
#
# Three deployment options in order of recommended use:
#   A. ShinyApps.io  — easiest, free tier, instant citable URL
#   B. Posit Connect — institutional, behind authentication
#   C. Docker        — maximum control, any cloud provider
# ============================================================


# ── Option A: ShinyApps.io (recommended for open research) ───
#
# Prerequisites:
#   1. Create a free account at https://www.shinyapps.io
#   2. Go to Account → Tokens → Show token
#   3. Fill in your credentials below

library(rsconnect)

rsconnect::setAccountInfo(
  name   = "YOUR_SHINYAPPS_USERNAME",   # ← replace
  token  = "YOUR_TOKEN",                # ← replace
  secret = "YOUR_SECRET"                # ← replace
)

# Deploy the app — produces a live URL:
# https://YOUR_USERNAME.shinyapps.io/limma_deg_explorer
rsconnect::deployApp(
  appDir  = ".",                        # current project folder
  appName = "limma_deg_explorer",
  launch.browser = FALSE
)

# ── Option B: Posit Connect (institutional) ───────────────────
#
# rsconnect::deployApp(
#   appDir = ".",
#   server = "your-posit-connect-hostname.com"
# )

# ── Option C: Docker (cloud or HPC) ──────────────────────────
#
# The Dockerfile below produces a self-contained image (~800 MB).
# Build and run locally first, then push to your cloud provider.
#
# docker build -t limma-shiny-r .
# docker run -p 3838:3838 limma-shiny-r
# Open: http://localhost:3838
#
# Dockerfile content (save as 'Dockerfile' in project root):
dockerfile_content <- '
FROM rocker/shiny:4.3.0

# System dependencies for Bioconductor packages
RUN apt-get update && apt-get install -y \\
    libcurl4-openssl-dev \\
    libssl-dev \\
    libxml2-dev \\
    libfontconfig1-dev \\
    && rm -rf /var/lib/apt/lists/*

# Install renv for reproducibility (Rule 8)
RUN R -e "install.packages(\'renv\', repos = \'https://cloud.r-project.org\')"

# Copy renv lockfile and restore packages before copying app code
# This layer is cached unless renv.lock changes — fast rebuilds.
WORKDIR /srv/shiny-server/limma_app
COPY renv.lock .
RUN R -e "renv::restore()"

# Copy application files
COPY . .

# Expose the default Shiny port
EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
'

writeLines(trimws(dockerfile_content), "Dockerfile")
message("Dockerfile written to project root.")

# ── Post-deployment checklist ─────────────────────────────────
# After deploying, verify the live URL:
#   [ ] Upload GSE83272_limma.csv → summary banner appears
#   [ ] All sliders are responsive (no timeout on enrichment)
#   [ ] Download Report produces a valid self-contained HTML
#   [ ] App URL is accessible without login (for ShinyApps free tier)
#   [ ] Record the URL in the paper's Data Availability statement
