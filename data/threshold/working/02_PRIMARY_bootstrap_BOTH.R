# ======================================================================
# Script:  02_PRIMARY_bootstrap_BOTH.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-13
# ======================================================================
#
# PURPOSE
# ---
# SCOPE: BOTH (exclusive bins = primary; cumulative thresholds = complementary).
# Primary analysis of the locked plan v2: estimands E1 (absolute, MMRM), E2 (ratio of
# marginal means, MMRM), E3 (clr shape change vs smooth null), E4 (dz), with ONE
# patient-level cluster bootstrap shared by all estimands, bins, thresholds and visits.
# Produces pointwise and simultaneous (max-t) intervals, inner/outer confidence sets
# and the pre-specified localisation verdict.
#
# INPUTS
# ---
#   MYOII_ACC_OUTCOMES_MASTER_V1.rds — via R/00_config.R
#   output/RData/qc_objects.RData    — must exist (run 01 first; 01 stops on any FAIL)
#
# ANALYSES
# ---
#   1. Observed estimates (all estimands)
#   2. Cluster bootstrap, B = B_BOOT (default 2000), checkpointed in chunks
#   3. Pointwise percentile CI, bootstrap SE, max-t band per estimand x block x visit
#   4. Inner / outer confidence sets on delta grids
#   5. Localisation decision (E3 residuals vs natural-spline null)
#
# OUTPUTS
# ---
#   RData:  output/RData/primary_objects.RData
#           output/RData/boot_chunks/*.rds (checkpoints; safe to delete after a full run)
#   Tables: output/tables/YYMMDD_02_primary_results.xlsx
#   Logs:   output/logs/YYMMDD_02_PRIMARY_bootstrap_BOTH_sessionInfo.txt
#
# RUNTIME
# ---
#   Windows, 31 workers: ~5.5 min per chunk of 100 replicates -> B = 2000 ~ 1.5-2 h.
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
for (f in c("01_helpers_data.R", "02_helpers_bootstrap.R", "03_helpers_estimands.R"))
  source(file.path(WORKING_DIR, "R", f))

qc_file <- file.path(PATHS$rdata, "qc_objects.RData")
if (!file.exists(qc_file)) stop("Run 01_QC_data_dictionary_BOTH.R first.")

acc <- load_master()
dl  <- build_analysis_long(acc)

# **********************************************************************
# 2. Observed estimates ----
# **********************************************************************
t_obs <- system.time(observed <- estimate_all(dl))
message("Observed estimates: ", nrow(observed), " rows in ", round(t_obs[3]), " s")

# **********************************************************************
# 3. Cluster bootstrap ----
# **********************************************************************
boot <- run_cluster_bootstrap(dl, estimate_all, B = B_BOOT, seed = SEED, workers = N_WORKERS,
                              checkpoint_dir = file.path(PATHS$rdata, "boot_chunks"), tag = "primary")
boot_settings <- attr(boot, "settings")
print(boot_settings)

# Fit-method log across replicates (fallbacks and failures per descriptor).
fit_log <- boot |>
  filter(!failed, estimand %in% c("E1", "E2")) |>
  count(estimand, block, descriptor, method, name = "n_replicate_visit") |>
  tidyr::pivot_wider(names_from = method, values_from = n_replicate_visit, values_fill = 0)

# **********************************************************************
# 4. Intervals and simultaneous bands ----
# **********************************************************************
keys <- c("estimand", "block", "visit", "descriptor")
summary_all <- summarise_bootstrap(observed, boot, keys = keys, band_by = c("estimand", "block", "visit")) |>
  left_join(dl |> distinct(block, descriptor, lower_mg, upper_mg, intensity) |>
              bind_rows(tibble::tibble(block = "bins", descriptor = ">=400", lower_mg = 400, upper_mg = Inf, intensity = NA_real_)) |>
              distinct(block, descriptor, .keep_all = TRUE),
            by = c("block", "descriptor")) |>
  mutate(visit = factor(visit, levels = FOLLOWUP_VISITS))
# E3 parts are defined on the composition grid; attach their bounds explicitly.
e3_bounds <- e3_composition(dl) |> distinct(descriptor, lower_mg, upper_mg)
summary_all <- summary_all |>
  rows_update(summary_all |> filter(block == "composition") |> select(-lower_mg, -upper_mg) |>
                left_join(e3_bounds, by = "descriptor"), by = keys, unmatched = "ignore")

# Ratio-scale presentation for E2 (exp of log-ratio and its bounds).
summary_all <- summary_all |>
  mutate(across(c(estimate, ci_low, ci_high, band_low, band_high),
                ~ ifelse(grepl("^E2", estimand), exp(.x), NA_real_), .names = "{.col}_ratio"))

# **********************************************************************
# 5. Confidence sets (where does change lie?) ----
# **********************************************************************
describe_set <- function(flag, labels) {
  flag[is.na(flag)] <- FALSE
  if (!any(flag)) return("none")
  r <- rle(flag); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  paste0(sub("-.*", "", labels[starts[r$values]]), "-", sub(".*-", "", labels[ends[r$values]]), collapse = "; ")
}
delta_grid <- list(
  E1_bins       = c(0, 0.25, 0.5, 1),       # min/day per 10 mg
  E1_cumulative = c(0, 1, 5),               # min/day above T
  E2_bins       = log(c(1, 1.1, 1.25)),     # ratio of marginal means
  E2_cumulative = log(c(1, 1.1, 1.25))
)
conf_sets <- purrr::imap_dfr(delta_grid, \(dg, nm) {
  est <- sub("_.*", "", nm); blk <- sub("^[^_]*_", "", nm)
  lev <- if (blk == "bins") c(descriptor_levels(dl, "bins"), ">=400") else descriptor_levels(dl, "cumulative")
  confidence_sets(summary_all |> filter(estimand == est, block == blk), dg) |>
    mutate(descriptor = factor(descriptor, levels = lev)) |> arrange(visit, delta, descriptor) |>
    group_by(estimand, block, visit, delta) |>
    summarise(inner_set = describe_set(inner, as.character(descriptor)),
              outer_set = describe_set(outer, as.character(descriptor)), .groups = "drop")
}) |> mutate(delta_display = ifelse(estimand == "E2", round(exp(delta), 3), delta))

# **********************************************************************
# 6. Localisation decision (plan v2 §3.4) ----
# **********************************************************************
e3_levels <- e3_bounds |> filter(is.finite(upper_mg), lower_mg > 0) |> arrange(lower_mg) |> pull(descriptor)
decision <- decide_localisation(
  e3_resid  = summary_all |> filter(estimand == "E3_residual"),
  e3_slope  = summary_all |> filter(estimand == "E3_slope"),
  bin_levels = e3_levels
)
print(decision, width = 200)

# **********************************************************************
# 7. Save ----
# **********************************************************************
save_tables(list(summary = summary_all |> select(-any_of("se_model")) |> left_join(observed |> select(all_of(keys), se_model), by = keys),
                 confidence_sets = conf_sets, decision = decision, fit_log = fit_log,
                 bootstrap_settings = boot_settings), "02_primary_results")
primary_objects <- list(observed = observed, summary = summary_all, conf_sets = conf_sets,
                        decision = decision, fit_log = fit_log, boot_settings = boot_settings,
                        delta_grid = delta_grid, e3_levels = e3_levels,
                        input = list(path = INPUT_RDS, md5 = INPUT_MD5, master = acc_meta_string()))
save_public("primary_objects", file = file.path(PATHS$rdata, "primary_objects.RData"))
saveRDS(boot, file.path(PATHS$rdata, "primary_boot_replicates.rds"))
write_run_meta("02_PRIMARY_bootstrap_BOTH.R", list(B = B_BOOT, seed = SEED, workers = N_WORKERS,
                                                   failed_replicates = boot_settings$n_failed_replicates))
message("02 PRIMARY finished.")
