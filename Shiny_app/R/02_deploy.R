# ============================================================
# 02_deploy.R  —  Publish the app to shinyapps.io
# Rule 9 — Deploy to a URL, not an inbox
#
# Author:      Patricia Agudelo-Romero
# Email:       patricia.agudeloromero@thekids.org.au
# Institution: The Kids Research Institute Australia
#
# Credentials are read from environment variables so they never end
# up in the repository. Get them at shinyapps.io -> Account -> Tokens,
# then add them to ~/.Renviron (not to this file):
#   SHINYAPPS_ACCOUNT=agudeloromero
#   SHINYAPPS_TOKEN=...
#   SHINYAPPS_SECRET=...
#
# Run from this folder (Shiny_app/R):
#   Rscript 02_deploy.R
#
# A Docker alternative is provided in ./Dockerfile.
# ============================================================

library(rsconnect)

account <- Sys.getenv("SHINYAPPS_ACCOUNT")
token   <- Sys.getenv("SHINYAPPS_TOKEN")
secret  <- Sys.getenv("SHINYAPPS_SECRET")
if (!nzchar(account) || !nzchar(token) || !nzchar(secret)) {
  stop("Set SHINYAPPS_ACCOUNT, SHINYAPPS_TOKEN and SHINYAPPS_SECRET in ~/.Renviron first.")
}
setAccountInfo(name = account, token = token, secret = secret)

# Only the files the running app needs are uploaded
app_files <- c(
  "app.R", "helpers.R",
  "templates/report_template.Rmd",
  "data/GSE60450_limma.csv"
)
stopifnot(all(file.exists(app_files)))

deployApp(
  appDir   = ".",
  appFiles = app_files,
  appName  = "limma-deg-explorer",
  account  = account,
  launch.browser = FALSE
)

# After deploying, check the live URL:
#   [ ] "load the example dataset" -> summary banner appears
#   [ ] sliders and numeric boxes stay in sync
#   [ ] Run GO enrichment -> dot plot and table
#   [ ] Download results folder -> zip with HTML, CSV and Excel
