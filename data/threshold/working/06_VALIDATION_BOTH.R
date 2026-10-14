# ======================================================================
# Script:  06_VALIDATION_BOTH.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-13
# ======================================================================
#
# PURPOSE
# ---
# SCOPE: BOTH. Pre-defined validation plan (phase 5): data/unit tests, model checks,
# reproducibility checks, a plasmode power simulation for the E3 localisation rule,
# and the issue log. Produces the PASS/WARN/FAIL table used in the report.
#
# INPUTS
# ---
#   MYOII_ACC_OUTCOMES_MASTER_V1.rds; output/RData/primary_objects.RData
#   (sensitivity_objects.RData optional)
#
# ANALYSES
# ---
#   V1 data integrity (qc_master: no FAIL)
#   V2 identity: cumulative paired change = sum of bin paired changes above T (+ >=400)
#   V3 no partial splitting of bins; composition sums to 1440 min
#   V4 reproducibility: estimate_all() twice identical; bootstrap 1 vs N workers identical
#   V5 model checks: MMRM convergence across replicates, residual skewness/kurtosis (RI fit)
#   V6 influence: leave-one-patient-out max |change| in E1 (bins)
#   V7 true-null calibration (type-I) and power of the E3 bump rule (2-6 bins, +0.1/+0.2 log)
#   V8 issue log
#
# OUTPUTS
# ---
#   RData:  output/RData/validation_objects.RData
#   Tables: output/tables/YYMMDD_06_validation.xlsx
#
# RUNTIME
# ---
#   V6 ~ 71 x 23 MMRM fits (parallel, ~5 min); V7 matrix engine, 8 scenarios x 3 visits x N_SIM x B_SIM (~minutes).
#   Lower N_SIM / B_SIM below for a quick check.
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
load(file.path(PATHS$rdata, "primary_objects.RData"))
sens_file <- file.path(PATHS$rdata, "sensitivity_objects.RData")
if (file.exists(sens_file)) load(sens_file)

N_SIM <- as.integer(Sys.getenv("MYOII_NSIM", unset = "500"))   # plasmode datasets per scenario
B_SIM <- as.integer(Sys.getenv("MYOII_BSIM", unset = "1000"))   # bootstrap replicates per dataset
RUN_LOO <- TRUE

acc <- load_master()
dl  <- build_analysis_long(acc)
val <- list()
vrow <- function(id, check, status, detail) tibble::tibble(id = id, check = check, status = status, detail = detail)

# **********************************************************************
# 2. V1–V3 data and unit tests ----
# **********************************************************************
qc <- qc_master(acc)
val$V1 <- vrow("V1", "Data integrity checks (D1–D19)", if (any(qc$status == "FAIL")) "FAIL" else "PASS",
               sprintf("%d PASS, %d WARN, %d FAIL", sum(qc$status == "PASS"), sum(qc$status == "WARN"), sum(qc$status == "FAIL")))

s1b <- sens_complete_pairs(dl, "bins"); s1c <- sens_complete_pairs(dl, "cumulative")
widths <- dl |> filter(block == "bins") |> distinct(descriptor, lower_mg, width_mg)
pairs_t400 <- s1c |> filter(descriptor == "T400") |> select(visit, t400 = E1_pairs)
recon <- tidyr::crossing(T = THRESHOLD_GRID_MG, visit = FOLLOWUP_VISITS) |>
  purrr::pmap_dfr(\(T, visit) {
    b <- s1b |> filter(visit == !!visit) |> inner_join(widths, by = "descriptor") |> filter(lower_mg >= T)
    tibble::tibble(T = T, visit = visit, sum_bins = sum(b$E1_pairs * b$width_mg / 10) + pairs_t400$t400[pairs_t400$visit == visit],
                   cumulative = s1c$E1_pairs[s1c$descriptor == paste0("T", T) & s1c$visit == visit])
  })
err_v2 <- max(abs(recon$sum_bins - recon$cumulative))
val$V2 <- vrow("V2", "Paired change: cumulative T = sum of bins >= T + (>=400)", if (err_v2 < 1e-8) "PASS" else "FAIL",
               sprintf("max |error| = %.2e min/day (identity holds for complete pairs; MMRM per descriptor is not additive by construction)", err_v2))

comp <- e3_composition(dl)
tot <- comp |> group_by(patient_id, visit_code) |> summarise(s = sum(value), r = sum(value_r), .groups = "drop")
val$V3 <- vrow("V3", "Composition closes to 1440 min (raw and after zero replacement); bins = native GGIR edges",
               if (max(abs(tot$s - 1440)) < 1e-6 && max(abs(tot$r - 1440)) < 1e-6) "PASS" else "FAIL",
               sprintf("max |sum-1440| raw %.1e, replaced %.1e; %d parts", max(abs(tot$s - 1440)), max(abs(tot$r - 1440)), n_distinct(comp$descriptor)))

# **********************************************************************
# 3. V4 reproducibility ----
# **********************************************************************
e_a <- est_E3(dl); e_b <- est_E3(dl)
stat_fast <- function(d) est_E3(d) |> filter(estimand == "E3_residual")
b1 <- run_cluster_bootstrap(dl[dl$block != "composition", ], stat_fast, B = 24, workers = 1)
bN <- run_cluster_bootstrap(dl[dl$block != "composition", ], stat_fast, B = 24, workers = min(4L, N_WORKERS))
same_est <- isTRUE(all.equal(e_a$estimate, e_b$estimate, tolerance = 1e-12))
same_boot <- isTRUE(all.equal(b1$estimate, bN$estimate, tolerance = 1e-10))
val$V4 <- vrow("V4", "Reproducibility: same input = same estimates; bootstrap 1 vs N workers identical",
               if (same_est && same_boot) "PASS" else "FAIL",
               sprintf("estimates identical: %s; 1 vs %d workers identical: %s; seed %d", same_est, min(4L, N_WORKERS), same_boot, SEED))

# **********************************************************************
# 4. V5 model checks ----
# **********************************************************************
fl <- primary_objects$fit_log
col_sum <- function(df, nm) if (nm %in% names(df)) sum(df[[nm]], na.rm = TRUE) else 0
n_fallback <- col_sum(fl, "LMM_ri"); n_failed <- col_sum(fl, "failed")
val$V5a <- vrow("V5a", "MMRM convergence across bootstrap replicates",
                if (n_failed > 0) "WARN" else if (n_fallback > 0) "WARN" else "PASS",
                sprintf("replicate x visit fits: random-intercept fallback %d, failed %d (see fit_log)", n_fallback, n_failed))
resid_diag <- dl |> filter(block == "bins") |> split(~descriptor) |>
  purrr::map_dfr(\(x) {
    f <- glmmTMB::glmmTMB(density ~ visit_code + (1 | patient_id), data = x, REML = TRUE)
    r <- residuals(f, type = "pearson")
    tibble::tibble(descriptor = x$descriptor[1], skewness = mean(((r - mean(r)) / sd(r))^3),
                   excess_kurtosis = mean(((r - mean(r)) / sd(r))^4) - 3,
                   heterosc_rho = suppressWarnings(cor(abs(r), fitted(f), method = "spearman")))
  })
val$V5b <- vrow("V5b", "Residual distribution (RI LMM, per bin)", if (any(abs(resid_diag$skewness) > 2)) "WARN" else "PASS",
                sprintf("max |skewness| %.2f, max excess kurtosis %.1f, max Spearman(|r|, fitted) %.2f -> uncertainty from cluster bootstrap, not Wald",
                        max(abs(resid_diag$skewness)), max(resid_diag$excess_kurtosis), max(resid_diag$heterosc_rho, na.rm = TRUE)))

# **********************************************************************
# 5. V6 leave-one-patient-out influence (E1 bins) ----
# **********************************************************************
# Parallel plan for V6/V7. Set explicitly (on.exit() at top level of a sourced script
# fires immediately and would silently fall back to sequential execution).
if (N_WORKERS > 1L) future::plan(future::multisession, workers = N_WORKERS) else future::plan(future::sequential)
message("[06] parallel plan: ", class(future::plan())[1], " with ", future::nbrOfWorkers(), " workers")

loo <- NULL
if (RUN_LOO) {
  message("[06] V6 leave-one-patient-out: ", n_distinct(dl$patient_id), " refits x 23 bins ...")
  pts <- sort(unique(dl$patient_id))
  e1_full <- primary_objects$observed |> filter(estimand == "E1", block == "bins") |> select(visit, descriptor, full = estimate, se_full = se_model)
  one_loo <- function(p) {
    d <- dl[dl$patient_id != p & dl$block == "bins", , drop = FALSE]; attr(d, "levels") <- attr(dl, "levels")
    est_E1_E2(d, "bins") |> filter(estimand == "E1") |> transmute(visit, descriptor, estimate, left_out = p)
  }
  loo <- furrr::future_map_dfr(pts, one_loo, .options = furrr::furrr_options(seed = NULL, packages = c("dplyr", "tidyr", "purrr", "tibble", "glmmTMB"))) |>
    inner_join(e1_full, by = c("visit", "descriptor")) |>
    group_by(visit, descriptor) |>
    summarise(max_abs_change = max(abs(estimate - full)), change_in_se = max(abs(estimate - full)) / first(se_full),
              .groups = "drop")
  # No patient identifiers are stored (outputs are shared publicly).
  # Influence is scaled by the model SE of the full-data estimate (relative change is unstable near 0).
  val$V6 <- vrow("V6", "Leave-one-patient-out influence on E1 (in SE units)", if (any(loo$change_in_se > 1)) "WARN" else "PASS",
                 sprintf("max shift %.2f SE (median %.2f SE); bins/visits with shift > 1 SE: %d",
                         max(loo$change_in_se), median(loo$change_in_se), sum(loo$change_in_se > 1)))
}

# **********************************************************************
# 6. V7 plasmode calibration and power of the E3 bump rule ----
# **********************************************************************
# Round-3 review (items 5–6): TRUE-NULL plasmode — each patient's paired clr change has
# the observed mean residual removed, so the population residual profile is exactly zero
# while the smooth tilt and the between-patient covariance are kept. Bumps of +0.1 / +0.2
# log are added to 2, 3, 4 or 6 adjacent bins. Detection requires a flagged positive bump
# that OVERLAPS the injected bins. Uses the matrix engine (R/05_helpers_e3fast.R), which
# resamples the complete pairs of each visit (V7a checks equality with est_E3()).
source(file.path(WORKING_DIR, "R", "05_helpers_e3fast.R"))
comp_e3 <- e3_composition(dl)
D_list <- purrr::set_names(purrr::map(FOLLOWUP_VISITS, \(v) e3_pair_matrix(comp_e3, v)), FOLLOWUP_VISITS)
labs_e3 <- attr(D_list[[1]], "labels")
R_e3 <- e3_projector(attr(D_list[[1]], "x"))

# V7a: matrix engine reproduces est_E3() residuals exactly.
chk <- primary_objects$observed |> filter(estimand == "E3_residual")
eq_err <- max(purrr::map_dbl(FOLLOWUP_VISITS, \(v) {
  r <- drop(R_e3 %*% colMeans(D_list[[v]]))
  o <- chk |> filter(visit == v) |> (\(x) x$estimate[match(labs_e3, x$descriptor)])()
  max(abs(r - o))
}))
val$V7a <- vrow("V7a", "Matrix E3 engine = est_E3() (observed residuals)", if (eq_err < 1e-10) "PASS" else "FAIL",
                sprintf("max |difference| = %.1e", eq_err))

bins_of <- function(lo, hi) which(as.numeric(sub("-.*", "", labs_e3)) >= lo & as.numeric(sub(".*-", "", labs_e3)) <= hi)
scen <- tibble::tribble(
  ~scenario,               ~lo, ~hi, ~size,
  "null (true null)",        NA,  NA,  0,
  "2 bins 120-160, +0.2",   120, 160, 0.2,
  "2 bins 120-160, +0.1",   120, 160, 0.1,
  "3 bins 180-240, +0.2",   180, 240, 0.2,
  "4 bins 80-120, +0.2",     80, 120, 0.2,
  "4 bins 80-120, +0.1",     80, 120, 0.1,
  "6 bins 70-140, +0.2",     70, 140, 0.2,
  "6 bins 70-140, +0.1",     70, 140, 0.1
)
scen$hit <- purrr::map2(scen$lo, scen$hi, \(a, b) if (is.na(a)) integer(0) else bins_of(a, b))
scen$n_bins <- lengths(scen$hit)
jobs <- tidyr::crossing(scen |> select(scenario, size, hit, n_bins), visit = FOLLOWUP_VISITS)
message(sprintf("[06] V7 plasmode: %d scenario x visit jobs, N_SIM = %d, B_SIM = %d — %s",
                nrow(jobs), N_SIM, B_SIM, format(Sys.time(), "%H:%M")))
power <- furrr::future_pmap_dfr(jobs, \(scenario, size, hit, n_bins, visit) {
  simulate_e3(D_list[[visit]], R_e3, hit = hit, size = size, n_sim = N_SIM, b_sim = B_SIM,
              seed = SEED + 1000L * match(visit, FOLLOWUP_VISITS)) |>
    mutate(scenario = scenario, visit = visit, n_bins = n_bins)
}, .options = furrr::furrr_options(seed = NULL, packages = c("dplyr", "tidyr", "purrr", "tibble", "splines")))
mc_ci <- function(k, n) { ci <- stats::binom.test(k, n)$conf.int; sprintf("%.0f-%.0f", 100 * ci[1], 100 * ci[2]) }
power_summary <- power |>
  group_by(scenario, visit, n_bins) |>
  summarise(n_sim = n(), n_pairs = NA_integer_,
            k_any = sum(any_bump), k_det = sum(detected), .groups = "drop") |>
  mutate(n_pairs = purrr::map_int(visit, \(v) nrow(D_list[[v]])),
         pct_any_bump = 100 * k_any / n_sim, any_bump_mc95 = purrr::map2_chr(k_any, n_sim, mc_ci),
         pct_detected = ifelse(n_bins > 0, 100 * k_det / n_sim, NA_real_),
         detected_mc95 = ifelse(n_bins > 0, purrr::map2_chr(k_det, n_sim, mc_ci), NA_character_),
         visit = factor(visit, levels = FOLLOWUP_VISITS),
         scenario = factor(scenario, levels = scen$scenario)) |>
  arrange(scenario, visit) |> select(-k_any, -k_det)
null_fp <- power_summary |> filter(scenario == "null (true null)")
val$V7 <- vrow("V7", sprintf("Power of E3 bump rule, true-null plasmode (N_SIM = %d, B_SIM = %d)", N_SIM, B_SIM), "INFO",
               paste(sprintf("%s %s: %.0f%%", power_summary$scenario, power_summary$visit, power_summary$pct_detected)[power_summary$n_bins > 0], collapse = "; "))
val$V7b <- vrow("V7b", "Type-I rate of the bump rule under a true null (nominal <= 5%)",
                if (any(null_fp$pct_any_bump > 5 & !grepl("^[0-4]-", null_fp$any_bump_mc95))) "WARN" else "PASS",
                paste(sprintf("%s: %.1f%% (MC 95%% CI %s)", null_fp$visit, null_fp$pct_any_bump, null_fp$any_bump_mc95), collapse = "; "))

future::plan(future::sequential)

# **********************************************************************
# 7. V8 issue log ----
# **********************************************************************
issue_log <- tibble::tribble(
  ~id,  ~observation, ~severity, ~cause, ~test, ~fix, ~impact,
  "K1", "Earlier claim: high-intensity M3 peak is an SD artefact", "critical (scientific)", "scale-dependent location of maxima", "relative (E2), dz (E4), clr (E3) profiles", "claim withdrawn; estimands chosen a priori", "interpretation H1 revised",
  "K2", "Revised R00 hard-codes R:/ path, ignores MYOII_MASTER_RDS", "critical (reproducibility)", "regression in revised helper", "code review", "new config with lookup order + env override", "old R05/R06 not portable",
  "K3", "season_anchor_test sheet empty in A04/A07", "high", "anova.merMod column is Df, not `Chi Df`", "xlsx inspection + code", "fixed in season_anchor_test()", "seasonality test now available (S5)",
  "V1", "Bootstrap results depend on worker count (R03/R06)", "high", "clusterSetRNGStream + parLapply", "1 vs N workers", "draws pre-generated in main process", "results portable across machines",
  "V2", "R05 bootstrap per visit x bin, not per patient across spectrum", "high", "group_modify resampling", "code review", "one draw for all bins/visits", "simultaneous bands now valid",
  "MO2", "Absolute peak at 120-140 mg is a bin-width artefact", "major", "bins 10/20/25 mg", "per-10 mg normalisation", "E1 per 10 mg", "absolute peak 90-100 mg",
  "MO4", "Zeros concentrated at M0 in 375-400 mg (14%)", "major", "sparse tail", "zero tables by visit", ">=300 mg amalgamated for E3; replacement grid", "E3 tail stable",
  "MO5", "Trajectory-dependent dropout to M18", "major", "observational follow-up", "S10", "MMRM (MAR) primary; tipping point S9", "M18 interpreted with caveat",
  "R1", "Linear trend null produces U-shaped residuals", "major", "S-shaped change profile", "observed residual signs", "natural-spline null, bump rule", "no false edge 'intervals'",
  "MO9", "L/R wrist labels differ within 10 patients", "info", "filename labels (device always non-dominant per user)", "filename parsing", "none (documented)", "none expected",
  "R3-5", "'Null' plasmode resampled observed data (not a true null); small N_SIM/B_SIM", "major", "simulation design", "true-null plasmode, N_SIM 500 x B_SIM 1000, MC CI", "V7 rewritten (matrix engine)", "type-I rate now estimable",
  "R3-6", "Spline df = 3 null may absorb broad bumps; detection matched edge strings", "major", "smooth null resolution", "4- and 6-bin scenarios; overlap-based detection", "V7 scenarios extended", "E3 conclusion limited to the tested resolution",
  "R3-7", "Tipping point E/p is a mechanical identity", "major", "definition on point estimate", "algebra", "report delta at which the 95% lower bound reaches 0", "M18 robustness re-stated",
  "R3-11", "IMNM single-bin highlight was post hoc", "major", "selection from 23 x 3 x 3 cells", "-", "pre-specified per-subtype summaries (T100, tilt, dIG) + full profiles", "subtype text reframed"
)

# **********************************************************************
# 8. Save ----
# **********************************************************************
validation_table <- bind_rows(val)
print(validation_table, n = Inf, width = 200)
save_tables(list(validation = validation_table, power = power_summary, residual_diagnostics = resid_diag,
                 loo_influence = if (!is.null(loo)) loo else tibble::tibble(),
                 V2_reconstruction = recon, issue_log = issue_log), "06_validation")
validation_objects <- list(validation_table = validation_table, power_summary = power_summary,
                           resid_diag = resid_diag, loo = loo, recon = recon, issue_log = issue_log,
                           settings = list(N_SIM = N_SIM, B_SIM = B_SIM))
save_public("validation_objects", file = file.path(PATHS$rdata, "validation_objects.RData"))
write_run_meta("06_VALIDATION_BOTH.R", list(N_SIM = N_SIM, B_SIM = B_SIM))
message("06 VALIDATION finished.")
