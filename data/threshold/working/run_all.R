# ======================================================================
# Script:  run_all.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-13
# ======================================================================
#
# PURPOSE
# ---
# Runs the full pipeline 01 -> 06 in a clean order and knits the report.
# Open this file in RStudio and Source it, or run each script separately.
#
# INPUTS
# ---
#   MYOII_ACC_OUTCOMES_MASTER_V1.rds (see README.md, "Kde má být master")
#
# ANALYSES
# ---
#   01 QC -> 02 PRIMARY (bootstrap) -> 03 SENSITIVITY -> 04 DESCRIPTIVE
#   -> 05 FIGURES -> 06 VALIDATION -> report
#
# OUTPUTS
# ---
#   everything under working/output/ and working/reports/*.html
#
# ======================================================================

# Optional overrides (uncomment before running):
# Sys.setenv(MYOII_MASTER_RDS = "C:/path/to/MYOII_ACC_OUTCOMES_MASTER_V1.rds")
# Sys.setenv(MYOII_B = 50)                         # quick test run (default 2000)
# Sys.setenv(MYOII_NSIM = 20, MYOII_BSIM = 100)    # quick power simulation (default 500 / 1000)

if (!exists("WORKING_DIR") || !file.exists(file.path(WORKING_DIR, "R", "00_config.R"))) WORKING_DIR <- local({
  # Locate working/: env MYOII_WORKING_DIR -> folder of the file being source()d ->
  # active RStudio document -> getwd(), getwd()/working, getwd()/data/threshold/working.
  ofile <- NULL
  for (i in rev(seq_len(sys.nframe()))) { f <- sys.frame(i)$ofile; if (!is.null(f)) { ofile <- f; break } }
  ed <- if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "") else ""
  cand <- c(Sys.getenv("MYOII_WORKING_DIR"), if (!is.null(ofile)) dirname(normalizePath(ofile)),
            if (nzchar(ed)) dirname(ed), getwd(), file.path(getwd(), "working"),
            file.path(getwd(), "data", "threshold", "working"))
  cand <- cand[nzchar(cand)]
  hit <- cand[file.exists(file.path(cand, "R", "00_config.R"))]
  if (!length(hit)) stop("Cannot locate the 'working' folder; set Sys.setenv(MYOII_WORKING_DIR = '<path>/working').")
  normalizePath(hit[1], winslash = "/")
})
source(file.path(WORKING_DIR, "R", "00_config.R"))

steps <- c("01_QC_data_dictionary_BOTH.R", "02_PRIMARY_bootstrap_BOTH.R", "03_SENSITIVITY_BOTH.R",
           "04_DESCRIPTIVE_heterogeneity_BINS.R", "05_FIGURES_TABLES_BOTH.R", "06_VALIDATION_BOTH.R")
timing <- list()
for (s in steps) {
  message("\n==================== ", s, " ====================")
  t0 <- Sys.time()
  source(file.path(WORKING_DIR, s), local = new.env())
  timing[[s]] <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
}
# 05 is re-run after 06 is not needed: the report reads validation_objects directly.

check_packages(PKG_REPORT)
rmarkdown::render(file.path(WORKING_DIR, "reports", "261013_thresholds_report.Rmd"),
                  output_format = "html_document", envir = new.env())
print(tibble::enframe(unlist(timing), name = "script", value = "minutes"))
