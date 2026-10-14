# ======================================================================
# Script:  03_SENSITIVITY_BOTH.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-13
# ======================================================================
#
# PURPOSE
# ---
# SCOPE: BOTH. Pre-specified sensitivity analyses S1–S10 of plan v2 and the
# invariance check required by the localisation rule (§3.4, condition 2):
# E3 bumps must persist across zero replacement {epoch, 0.01, 0.5 min/day},
# >= 5 valid days (S7) and exclusion of the UCTD/SSc patient (S8).
#
# INPUTS
# ---
#   MYOII_ACC_OUTCOMES_MASTER_V1.rds — via R/00_config.R
#   output/RData/primary_objects.RData — from 02_PRIMARY_bootstrap_BOTH.R
#
# ANALYSES
# ---
#   1. S1 complete pairs; S2 random-intercept LMM; S3 Tweedie; S4 log(y + c)
#   2. S5 season (factor, harmonics) + corrected anchor test
#   3. S6 time since M0 (natural spline); S7 >= 5 valid days; S8 without UCTD/SSc (E1/E2)
#   4. S9 tipping point M18; S10 M3 change by later M18 availability
#   5. E3 invariance bootstrap (same draws/seed as primary) and final verdict
#
# OUTPUTS
# ---
#   RData:  output/RData/sensitivity_objects.RData
#   Tables: output/tables/YYMMDD_03_sensitivity_results.xlsx
#   Logs:   output/logs/YYMMDD_03_SENSITIVITY_BOTH_sessionInfo.txt
#
# RUNTIME
# ---
#   model sensitivities ~5-10 min (single core); E3 invariance bootstrap is closed-form
#   (~1 s per replicate per worker).
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
for (f in c("01_helpers_data.R", "02_helpers_bootstrap.R", "03_helpers_estimands.R", "04_helpers_sensitivity.R"))
  source(file.path(WORKING_DIR, "R", f))

prim_file <- file.path(PATHS$rdata, "primary_objects.RData")
if (!file.exists(prim_file)) stop("Run 02_PRIMARY_bootstrap_BOTH.R first.")
load(prim_file)

acc <- load_master()
dl  <- build_analysis_long(acc)
subset_dl <- function(dl, keep) { out <- dl[keep, , drop = FALSE]; attr(out, "levels") <- attr(dl, "levels"); out }
dl_s7 <- subset_dl(dl, dl$n_valid_days >= 5)
dl_s8 <- subset_dl(dl, !(dl$subtype_code %in% "UCTD/SSc?"))

# **********************************************************************
# 2. Model-based sensitivities (observed data, Wald CI) ----
# **********************************************************************
message("S1-S2 ...")
s1 <- bind_rows(sens_complete_pairs(dl, "bins"), sens_complete_pairs(dl, "cumulative"))
s2 <- bind_rows(est_E1_E2(dl, "bins", "ri"), est_E1_E2(dl, "cumulative", "ri"))
message("S3 Tweedie (slow) ...")
s3 <- bind_rows(est_E2_tweedie(dl, "composition"), est_E2_tweedie(dl, "cumulative")) |>
  mutate(ci_low = estimate - 1.96 * se_model, ci_high = estimate + 1.96 * se_model)
message("S4 log(y + c) ...")
s4 <- bind_rows(sens_log_c(dl, "bins"), sens_log_c(dl, "cumulative"))
message("S5 season ...")
s5 <- bind_rows(sens_season(dl, "bins", "factor"), sens_season(dl, "bins", "harmonic"),
                sens_season(dl, "cumulative", "factor"), sens_season(dl, "cumulative", "harmonic"))
s5_anchor <- bind_rows(season_anchor_test(dl, "cumulative", "T100"), season_anchor_test(dl, "bins", "100-110"))
message("S6 time spline ...")
s6 <- bind_rows(sens_time_spline(dl, "bins"), sens_time_spline(dl, "cumulative"))
message("S7-S8 subsets ...")
s7 <- bind_rows(est_E1_E2(dl_s7, "bins"), est_E1_E2(dl_s7, "cumulative")) |> mutate(sensitivity = "S7 >= 5 valid days")
s8 <- bind_rows(est_E1_E2(dl_s8, "bins"), est_E1_E2(dl_s8, "cumulative")) |> mutate(sensitivity = "S8 without UCTD/SSc")
message("S9-S10 ...")
s9  <- bind_rows(sens_tipping_m18(dl, "bins"), sens_tipping_m18(dl, "cumulative"))
s10 <- bind_rows(sens_dropout_m3(dl, "bins"), sens_dropout_m3(dl, "cumulative"))

# Primary vs sensitivity comparison table for E1/E2 (same descriptor x visit).
prim <- primary_objects$summary |> filter(estimand %in% c("E1", "E2")) |>
  transmute(estimand, block, visit = as.character(visit), descriptor, primary = estimate, primary_ci_low = ci_low, primary_ci_high = ci_high)
compare_E <- bind_rows(
  s1 |> transmute(block, visit, descriptor, E1 = E1_pairs, E2 = E2_pairs, sensitivity = "S1 complete pairs") |>
    tidyr::pivot_longer(c(E1, E2), names_to = "estimand", values_to = "estimate"),
  s2 |> transmute(estimand, block, visit, descriptor, estimate, sensitivity = "S2 random intercept"),
  s7 |> transmute(estimand, block, visit, descriptor, estimate, sensitivity),
  s8 |> transmute(estimand, block, visit, descriptor, estimate, sensitivity)
) |> inner_join(prim, by = c("estimand", "block", "visit", "descriptor")) |>
  mutate(same_sign = sign(estimate) == sign(primary),
         within_primary_ci = estimate >= primary_ci_low & estimate <= primary_ci_high)
compare_summary <- compare_E |> group_by(estimand, block, sensitivity) |>
  summarise(n = n(), pct_same_sign = 100 * mean(same_sign, na.rm = TRUE),
            pct_within_primary_ci = 100 * mean(within_primary_ci, na.rm = TRUE), .groups = "drop")

# **********************************************************************
# 3. E3 invariance bootstrap (localisation rule, condition 2) ----
# **********************************************************************
E3_VARIANTS <- list(primary_epoch = list(zero = "epoch", sub = "all"),
                    zero_0.01     = list(zero = 0.01,    sub = "all"),
                    zero_0.5      = list(zero = 0.5,     sub = "all"),
                    S7_validdays5 = list(zero = "epoch", sub = "S7"),
                    S8_no_UCTD    = list(zero = "epoch", sub = "S8"))
e3_variants_stat <- function(d) {
  purrr::imap_dfr(E3_VARIANTS, \(v, nm) {
    dd <- switch(v$sub, all = d, S7 = d[d$n_valid_days >= 5, , drop = FALSE],
                 S8 = d[!(d$subtype_code %in% "UCTD/SSc?"), , drop = FALSE])
    est_E3(dd, zero = v$zero) |> filter(estimand %in% c("E3_residual", "E3_slope")) |> mutate(variant = nm)
  })
}
dl_e3 <- dl[dl$block %in% c("bins", "cumulative"), , drop = FALSE]
e3_obs  <- e3_variants_stat(dl_e3)
e3_boot <- run_cluster_bootstrap(dl_e3, e3_variants_stat, B = B_BOOT, seed = SEED, workers = N_WORKERS,
                                 checkpoint_dir = file.path(PATHS$rdata, "boot_chunks"), tag = "e3inv", chunk_size = 500L)
e3_summ <- summarise_bootstrap(e3_obs, e3_boot, keys = c("variant", "estimand", "visit", "descriptor"),
                               band_by = c("variant", "estimand", "visit"))
e3_decisions <- purrr::map_dfr(names(E3_VARIANTS), \(nm) {
  s <- e3_summ |> filter(variant == nm)
  decide_localisation(s |> filter(estimand == "E3_residual"), s |> filter(estimand == "E3_slope"),
                      primary_objects$e3_levels) |> mutate(variant = nm)
})
# Final verdict: a bump is claimed only if the primary bump is reproduced (identical
# label) in every variant; otherwise the descriptive outcome is reported.
final_verdict <- e3_decisions |>
  group_by(visit) |>
  summarise(primary_bumps = paste(c(positive_bumps[variant == "primary_epoch"], negative_bumps[variant == "primary_epoch"]) |>
                                    (\(z) z[nzchar(z)])(), collapse = "; "),
            invariant = all(positive_bumps == positive_bumps[variant == "primary_epoch"] &
                              negative_bumps == negative_bumps[variant == "primary_epoch"]),
            slope_positive_all = all(slope_ci_low > 0),
            .groups = "drop") |>
  mutate(final = case_when(
    nzchar(primary_bumps) & invariant ~ paste("localised departure:", primary_bumps),
    nzchar(primary_bumps) & !invariant ~ "candidate bump not invariant -> no evidence of a localised departure",
    slope_positive_all ~ "no evidence of a localised departure; smooth redistribution towards higher intensities",
    TRUE ~ "no evidence of a localised departure; shape change not identifiable"))
print(final_verdict, width = 200)

# **********************************************************************
# 4. Save ----
# **********************************************************************
save_tables(list(compare_summary = compare_summary, compare_detail = compare_E, S1_pairs = s1, S2_ri = s2,
                 S3_tweedie = s3, S4_log_c = s4, S5_season = s5, S5_anchor = s5_anchor, S6_time = s6,
                 S7_validdays = s7, S8_no_uctd = s8, S9_tipping = s9, S10_dropout = s10,
                 E3_invariance = e3_decisions, final_verdict = final_verdict,
                 E3_variant_bands = e3_summ), "03_sensitivity_results")
sensitivity_objects <- list(s1 = s1, s2 = s2, s3 = s3, s4 = s4, s5 = s5, s5_anchor = s5_anchor, s6 = s6,
                            s7 = s7, s8 = s8, s9 = s9, s10 = s10, compare_E = compare_E,
                            compare_summary = compare_summary, e3_summ = e3_summ,
                            e3_decisions = e3_decisions, final_verdict = final_verdict,
                            e3_boot_settings = attr(e3_boot, "settings"))
save_public("sensitivity_objects", file = file.path(PATHS$rdata, "sensitivity_objects.RData"))
write_run_meta("03_SENSITIVITY_BOTH.R", list(B = B_BOOT, seed = SEED, workers = N_WORKERS))
message("03 SENSITIVITY finished.")
