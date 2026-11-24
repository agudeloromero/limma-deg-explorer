# Run from Shiny_app/R:   Rscript tests/testthat.R
library(testthat)
test_dir(file.path("tests", "testthat"), stop_on_failure = TRUE)
