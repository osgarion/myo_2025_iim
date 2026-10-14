# ======================================================================
# Script:  01_QC_data_dictionary_BOTH.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-13
# ======================================================================
#
# PURPOSE
# ---
# SCOPE: BOTH (cumulative thresholds + exclusive bins).
# Machine-generated data dictionary, locked data-integrity checks
# (PASS/WARN/FAIL), cohort flow and descriptive spectrum profiles.
# Stops if any FAIL is found — later scripts must not run on a failed input.
#
# INPUTS
# ---
#   MYOII_ACC_OUTCOMES_MASTER_V1.rds — via R/00_config.R
#
# ANALYSES
# ---
#   1. Data dictionary of all master components
#   2. QC checks D1–D19
#   3. Cohort flow and visit patterns
#   4. Descriptive profiles (median, P10–P90) and coverage by visit
#
# OUTPUTS
# ---
#   RData:   output/RData/qc_objects.RData
#   Tables:  output/tables/YYMMDD_01_QC_data_dictionary.xlsx
#   Figures: output/figures/YYMMDD_01_descriptive_profile_{bins,cumulative}.png/.tiff
#   Logs:    output/logs/YYMMDD_01_QC_data_dictionary_BOTH_sessionInfo.txt
#
# ======================================================================

# **********************************************************************
# 1. Setup ----
# **********************************************************************
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
source(file.path(WORKING_DIR, "R", "01_helpers_data.R"))

acc <- load_master()

# **********************************************************************
# 2. Data dictionary ----
# **********************************************************************
describe_df <- function(df, component) {
  purrr::imap_dfr(df, \(x, nm) tibble::tibble(
    component = component, variable = nm, class = class(x)[1], n_rows = length(x),
    n_missing = sum(is.na(x)), n_unique = dplyr::n_distinct(x),
    min = if (is.numeric(x) && any(is.finite(x))) signif(min(x, na.rm = TRUE), 5) else NA_real_,
    max = if (is.numeric(x) && any(is.finite(x))) signif(max(x, na.rm = TRUE), 5) else NA_real_
  ))
}
data_dictionary <- purrr::map_dfr(
  c("record_level", "threshold_long", "bin_long", "acc_outcome_long", "cohort"),
  \(nm) describe_df(acc[[nm]], nm)
)
units_key <- tibble::tribble(
  ~variable,          ~unit,                ~meaning,
  "minutes_above",    "min/day",            "Time above cumulative ENMO threshold T (overlapping exposure)",
  "threshold_mg",     "mg",                 "Cumulative ENMO threshold",
  "minutes_bin",      "min/day",            "Time within exclusive bin [lower_mg, upper_mg)",
  "minutes_per_10mg", "min/day per 10 mg",  "minutes_bin normalised to 10 mg width (density, for coverage)",
  "lower_mg",         "mg",                 "Lower bin edge",
  "upper_mg",         "mg",                 "Upper bin edge",
  "bin_width_mg",     "mg",                 "Bin width (10/20/25 mg)"
)
data_dictionary <- data_dictionary |> left_join(units_key, by = "variable")

# **********************************************************************
# 3. QC checks ----
# **********************************************************************
qc_table <- qc_master(acc)
print(qc_table, n = Inf, width = 200)
if (any(qc_table$status == "FAIL")) stop("Data QC FAIL — see qc_table; downstream analyses must not run.")

# **********************************************************************
# 4. Cohort flow ----
# **********************************************************************
rec <- record_info(acc)
flow <- tibble::tibble(
  step = c("Records in master", "QC-eligible records", "Main-visit records (M0-M18)",
           "Patients in main analysis", "Exploratory records (M30-M60)"),
  n = c(nrow(rec), sum(rec$eligible %in% TRUE),
        sum(rec$eligible %in% TRUE & rec$visit_code %in% MAIN_VISITS),
        n_distinct(rec$patient_id[rec$eligible %in% TRUE & rec$visit_code %in% MAIN_VISITS]),
        sum(rec$eligible %in% TRUE & rec$visit_code %in% EXPLORATORY_VISITS))
)
visit_patterns <- rec |>
  filter(eligible %in% TRUE, visit_code %in% MAIN_VISITS) |>
  group_by(patient_id) |>
  summarise(pattern = paste(MAIN_VISITS[MAIN_VISITS %in% visit_code], collapse = "-"), .groups = "drop") |>
  count(pattern, name = "n_patients", sort = TRUE)
pairs_n <- tibble::tibble(visit = FOLLOWUP_VISITS,
  n_pairs = vapply(FOLLOWUP_VISITS, \(v) sum(tapply(rec$visit_code, rec$patient_id, \(z) all(c("M0", v) %in% z))), numeric(1)))
subtype_n <- rec |> filter(visit_code %in% MAIN_VISITS) |> distinct(patient_id, subtype_code, subtype) |>
  count(subtype_code, subtype, name = "n_patients")

# **********************************************************************
# 5. Descriptive profiles + coverage ----
# **********************************************************************
spectra <- list(bins = prepare_spectrum(acc, "bins"), cumulative = prepare_spectrum(acc, "cumulative"))

describe_profile <- function(d) {
  d |> group_by(kind, visit_code, descriptor, intensity) |>
    summarise(n = n(), p10 = quantile(value, .10), p25 = quantile(value, .25), median = median(value),
              p75 = quantile(value, .75), p90 = quantile(value, .90), mean = mean(value),
              pct_zero = 100 * mean(value == 0), .groups = "drop")
}
describe_coverage <- function(d) {
  tidyr::crossing(d |> select(kind, visit_code, descriptor, intensity, density), limit = COVERAGE_LIMITS) |>
    group_by(kind, visit_code, descriptor, intensity, limit) |>
    summarise(n = n(), pct_above = 100 * mean(density > limit), .groups = "drop")
}
profiles <- purrr::map_dfr(spectra, describe_profile)
coverage <- purrr::map_dfr(spectra, describe_coverage)

plot_profile <- function(k) {
  dd <- profiles |> filter(kind == k)
  ylab <- if (k == "bins") "Minutes/day within bin" else "Minutes/day above threshold"
  ggplot(dd, aes(intensity, median, linetype = visit_code, group = visit_code)) +
    geom_ribbon(aes(ymin = p10, ymax = p90), fill = "grey50", alpha = .08, colour = NA) +
    geom_line(linewidth = .8) + geom_point(aes(shape = visit_code), size = 1.6) +
    scale_visit_linetype() + ggplot2::scale_shape_manual(values = VISIT_SHAPES) +
    scale_y_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500)) +
    x_axis_mg() +
    labs(title = paste("Descriptive spectrum profile:", toupper(k)),
         subtitle = "Median (line) and P10-P90 (band) per visit; y on log(1+x) scale",
         x = if (k == "bins") "Exclusive bin midpoint (mg)" else "Cumulative threshold (mg)",
         y = ylab, linetype = "Visit", shape = "Visit") +
    theme_myo()
}
fig_profile <- purrr::map(c(bins = "bins", cumulative = "cumulative"), plot_profile)
purrr::iwalk(fig_profile, \(p, k) save_figure(p, paste0("01_descriptive_profile_", k)))

# **********************************************************************
# 6. Save ----
# **********************************************************************
save_tables(list(qc = qc_table, flow = flow, visit_patterns = visit_patterns, pairs = pairs_n,
                 subtypes = subtype_n, profiles = profiles, coverage = coverage,
                 data_dictionary = data_dictionary), "01_QC_data_dictionary")
qc_objects <- list(qc_table = qc_table, flow = flow, visit_patterns = visit_patterns, pairs_n = pairs_n,
                   subtype_n = subtype_n, profiles = profiles, coverage = coverage,
                   data_dictionary = data_dictionary, fig_profile = fig_profile,
                   input = list(path = INPUT_RDS, md5 = INPUT_MD5, master = acc_meta_string()))
save_public("qc_objects", file = file.path(PATHS$rdata, "qc_objects.RData"))
write_run_meta("01_QC_data_dictionary_BOTH.R")
message("01 QC finished: ", sum(qc_table$status == "WARN"), " WARN, 0 FAIL")
