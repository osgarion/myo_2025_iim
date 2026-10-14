# ======================================================================
# Script:  04_DESCRIPTIVE_heterogeneity_BINS.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-13
# ======================================================================
#
# PURPOSE
# ---
# SCOPE: BINS (exclusive bins; the ">=300"/"<20" parts only where stated).
# Additional pre-specified questions of plan v2:
#   C3 (priority 1) — is change a stable shift or driven by a minority?
#   C1              — does the SHAPE of change differ between M3 and M6? (within patients)
#   C5              — is change consistent with a smooth tilt summarised by IG / mean ENMO?
#   Subtypes (DM / IMNM / ASyS; PM and UCTD/SSc excluded) and M30–M60 — descriptive only.
# No clinical variables are modelled.
#
# INPUTS
# ---
#   MYOII_ACC_OUTCOMES_MASTER_V1.rds — via R/00_config.R
#
# ANALYSES
# ---
#   1. C3: median change, % patients with change > 0, share of total gain from top-10 % patients
#   2. C1: paired difference of clr-change profiles M6 vs M3 (patients with M0, M3, M6),
#          cluster bootstrap with max-t band
#   3. C5: patient-level clr slope vs change in intensity gradient and mean ENMO
#   4. Subtype (pre-specified summaries: T100, tilt, dIG; full bin profiles) and exploratory visits
#
# OUTPUTS
# ---
#   RData:  output/RData/descriptive_objects.RData
#   Tables: output/tables/YYMMDD_04_descriptive_heterogeneity.xlsx
#   Logs:   output/logs/YYMMDD_04_DESCRIPTIVE_heterogeneity_BINS_sessionInfo.txt
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

acc <- load_master()
dl  <- build_analysis_long(acc)
bins <- dl[dl$block == "bins", , drop = FALSE]
bin_levels <- descriptor_levels(dl, "bins")

# **********************************************************************
# 2. C3: heterogeneity / concentration of change ----
# **********************************************************************
pairs <- make_pairs(bins |> mutate(descriptor = factor(descriptor, levels = bin_levels))) |>
  mutate(delta10 = (value - baseline) * 10 / width_mg)
top_share <- function(x, q = 0.10) {
  gain <- sum(x); if (gain <= 0) return(NA_real_)
  k <- max(1L, ceiling(q * length(x)))
  sum(sort(x, decreasing = TRUE)[seq_len(k)]) / gain
}
c3 <- pairs |> group_by(visit, descriptor, intensity) |>
  summarise(n_pairs = n(), mean_delta10 = mean(delta10), median_delta10 = median(delta10),
            q25_delta10 = quantile(delta10, .25), q75_delta10 = quantile(delta10, .75),
            pct_positive = 100 * mean(delta10 > 0), top10_share_of_gain = top_share(delta10),
            .groups = "drop")

# **********************************************************************
# 3. C1: shape of change, M6 vs M3 (same patients) ----
# **********************************************************************
c1_stat <- function(d) {
  comp <- e3_composition(d)
  w <- comp |> select(patient_id, visit_code, descriptor, lower_mg, x_log_mg, clr) |>
    tidyr::pivot_wider(names_from = visit_code, values_from = clr) |>
    filter(!is.na(M0), !is.na(M3), !is.na(M6))
  w |> group_by(descriptor, lower_mg) |>
    summarise(estimate = mean((M6 - M0) - (M3 - M0)), n_patients = n(), .groups = "drop") |>
    mutate(contrast = "M6 minus M3 (clr change)")
}
c1_obs  <- c1_stat(dl)
c1_boot <- run_cluster_bootstrap(dl[dl$block %in% c("bins", "cumulative"), ], c1_stat, B = B_BOOT, seed = SEED,
                                 workers = N_WORKERS, checkpoint_dir = file.path(PATHS$rdata, "boot_chunks"),
                                 tag = "c1", chunk_size = 500L)
c1 <- summarise_bootstrap(c1_obs, c1_boot, keys = c("contrast", "descriptor"), band_by = "contrast") |>
  arrange(lower_mg)

# **********************************************************************
# 4. C5: is a smooth tilt (IG-like) enough? ----
# **********************************************************************
rl <- acc$record_level |> tibble::as_tibble() |>
  transmute(patient_id = as.character(patient_id), visit_code = as.character(visit_code),
            ig_gradient = .data[["part2_summary_AD_ig_gradient_ENMO_0-24hr"]],
            enmo_mean = .data[["part2_summary_ENMO_fullRecordingMean"]])
comp <- e3_composition(dl)
pslope <- comp |> filter(is.finite(x_log_mg)) |>
  select(patient_id, visit_code, descriptor, x_log_mg, clr) |>
  tidyr::pivot_wider(names_from = visit_code, values_from = clr) |>
  tidyr::pivot_longer(any_of(FOLLOWUP_VISITS), names_to = "visit", values_to = "clr_fu") |>
  filter(!is.na(M0), !is.na(clr_fu)) |>
  group_by(patient_id, visit) |>
  summarise(clr_slope = unname(coef(lm(I(clr_fu - M0) ~ x_log_mg))[2]), .groups = "drop")
dscal <- rl |> filter(visit_code %in% MAIN_VISITS) |>
  tidyr::pivot_wider(names_from = visit_code, values_from = c(ig_gradient, enmo_mean)) |>
  tidyr::pivot_longer(-patient_id, names_to = c(".value", "visit"), names_pattern = "(ig_gradient|enmo_mean)_(.*)") |>
  group_by(patient_id) |>
  mutate(d_ig = ig_gradient - ig_gradient[visit == "M0"][1], d_enmo = enmo_mean - enmo_mean[visit == "M0"][1]) |>
  ungroup() |> filter(visit != "M0")
c5 <- pslope |> inner_join(dscal, by = c("patient_id", "visit")) |> filter(is.finite(d_ig)) |>
  group_by(visit = factor(visit, levels = FOLLOWUP_VISITS)) |>
  summarise(n = n(), mean_d_ig = mean(d_ig), mean_d_enmo = mean(d_enmo), mean_clr_slope = mean(clr_slope),
            spearman_slope_vs_d_ig = cor(clr_slope, d_ig, method = "spearman"),
            spearman_slope_vs_d_enmo = cor(clr_slope, d_enmo, method = "spearman"), .groups = "drop")

# **********************************************************************
# 5. Subtypes and exploratory visits (descriptive) ----
# **********************************************************************
subtype_profile <- pairs |> filter(!is.na(subtype)) |>
  group_by(subtype, visit, descriptor, intensity) |>
  summarise(n_pairs = n(), mean_delta10 = mean(delta10), median_delta10 = median(delta10),
            ratio_of_means = mean(value) / mean(baseline), .groups = "drop")
# Pre-specified per-subtype summaries (round-3 review, item 11): T100 change (min/day),
# patient-level clr tilt and change in intensity gradient — complete pairs, descriptive.
subj_subtype <- dl |> distinct(patient_id, subtype) |> filter(!is.na(subtype))
t100 <- dl |> filter(block == "cumulative", descriptor == "T100") |>
  select(patient_id, visit_code, value) |>
  tidyr::pivot_wider(names_from = visit_code, values_from = value) |>
  tidyr::pivot_longer(any_of(FOLLOWUP_VISITS), names_to = "visit", values_to = "fu") |>
  filter(!is.na(M0), !is.na(fu)) |> transmute(patient_id, visit, d_t100 = fu - M0)
subtype_summary <- t100 |>
  full_join(pslope, by = c("patient_id", "visit")) |>
  full_join(dscal |> select(patient_id, visit, d_ig), by = c("patient_id", "visit")) |>
  inner_join(subj_subtype, by = "patient_id") |>
  mutate(visit = factor(visit, levels = FOLLOWUP_VISITS)) |>
  group_by(subtype, visit) |>
  summarise(n_pairs = sum(!is.na(d_t100)),
            t100_mean = mean(d_t100, na.rm = TRUE), t100_median = median(d_t100, na.rm = TRUE),
            t100_q25 = quantile(d_t100, .25, na.rm = TRUE), t100_q75 = quantile(d_t100, .75, na.rm = TRUE),
            t100_pct_gain = 100 * mean(d_t100 > 0, na.rm = TRUE),
            tilt_mean = mean(clr_slope, na.rm = TRUE), tilt_median = median(clr_slope, na.rm = TRUE),
            d_ig_mean = mean(d_ig, na.rm = TRUE), d_ig_median = median(d_ig, na.rm = TRUE), .groups = "drop") |>
  arrange(visit, subtype)

expl <- prepare_spectrum(acc, "bins", visits = c(MAIN_VISITS, EXPLORATORY_VISITS)) |>
  filter(visit_code %in% EXPLORATORY_VISITS | visit_code == "M0") |>
  group_by(visit_code, descriptor, intensity) |>
  summarise(n = n(), median_per10 = median(density), p25 = quantile(density, .25), p75 = quantile(density, .75),
            .groups = "drop")

# **********************************************************************
# 6. Save ----
# **********************************************************************
save_tables(list(C3_heterogeneity = c3, C1_shape_M6_vs_M3 = c1, C5_tilt_vs_IG = c5,
                 subtype_profiles = subtype_profile, subtype_summary = subtype_summary,
                 exploratory_visits = expl), "04_descriptive_heterogeneity")
descriptive_objects <- list(c3 = c3, c1 = c1, c5 = c5, subtype_profile = subtype_profile,
                            subtype_summary = subtype_summary, expl = expl,
                            c1_boot_settings = attr(c1_boot, "settings"))
save_public("descriptive_objects", file = file.path(PATHS$rdata, "descriptive_objects.RData"))
write_run_meta("04_DESCRIPTIVE_heterogeneity_BINS.R", list(B = B_BOOT, seed = SEED))
message("04 DESCRIPTIVE finished.")
