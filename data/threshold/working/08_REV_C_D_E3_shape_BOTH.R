# ======================================================================
# Script:  08_REV_C_D_E3_shape_BOTH.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-16
# ======================================================================
#
# PURPOSE
# ---
# SCOPE: BOTH (E3 composition uses bins + cumulative T20 / T300 / T400).
# A07 revision, blocks C and D (docs/261016_A07_critical_review.md).
#   C — verifies the E3 construction (closure, zeros, invariance to bin width and to the
#       clr centre) and tests whether the localisation verdict depends on the smooth-null
#       specification (df 1-5, width-weighted fit) or on the bin grid (uniform 20 mg,
#       bounded bins up to 400 mg). Plasmode for the choice of df, including a
#       mis-specified null (true smooth df = 5, tested with df = 2-5).
#   D — formal paired comparison of the shape of change between M3, M6 and M18
#       (bin-wise, linear tilt, global 3-df smooth difference), Holm-adjusted, plus the
#       common sample with all four visits and a timing / season description.
# The primary E3 estimates (02) are not changed; F3 is re-drawn with revised explanations.
#
# INPUTS
# ---
#   MYOII_ACC_OUTCOMES_MASTER_V1.rds — via R/00_config.R
#   output/RData/primary_objects.RData, descriptive_objects.RData — originals
#
# ANALYSES
# ---
#   1. C: construction checks and invariance (numerical)
#   2. C: specification / grid variants, matrix cluster bootstrap B = B_BOOT, bump rule
#   3. C: plasmode type-I / power by df (N_SIM x B_SIM, env MYOII_NSIM / MYOII_BSIM)
#   4. D: shape contrasts M3-M6, M3-M18, M6-M18 (+ common sample), timing and season
#   5. Figures (F3v2, variants, df simulation, shape contrasts)
#
# OUTPUTS
# ---
#   RData:   output/RData/rev_CD_objects.RData
#   Tables:  output/tables/YYMMDD_08_REV_C_D_E3_shape.xlsx
#   Figures: output/figures/YYMMDD_08_C1_F3v2_E3_shape.png/.tiff
#            output/figures/YYMMDD_08_C2_E3_variants.png/.tiff
#            output/figures/YYMMDD_08_C3_df_simulation.png/.tiff
#            output/figures/YYMMDD_08_D1_shape_contrasts.png/.tiff
#   Logs:    output/logs/YYMMDD_08_REV_C_D_E3_shape_BOTH_sessionInfo.txt
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
for (f in c("01_helpers_data.R", "02_helpers_bootstrap.R", "03_helpers_estimands.R",
            "05_helpers_e3fast.R", "06_helpers_revision.R"))
  source(file.path(WORKING_DIR, "R", f))
for (o in c("primary_objects", "descriptive_objects")) {
  f <- file.path(PATHS$rdata, paste0(o, ".RData"))
  if (!file.exists(f)) stop("Missing ", f, " — run 02 and 04 first.")
  load(f)
}
N_SIM <- as.integer(Sys.getenv("MYOII_NSIM", unset = "500"))
B_SIM <- as.integer(Sys.getenv("MYOII_BSIM", unset = "1000"))

acc <- load_master()
dl  <- build_analysis_long(acc)
S   <- primary_objects$summary

# **********************************************************************
# 2. C — construction checks and invariance ----
# **********************************************************************
comp <- e3_composition(dl)
D_list <- purrr::set_names(purrr::map(FOLLOWUP_VISITS, \(v) e3_pair_matrix(comp, v)), FOLLOWUP_VISITS)
x_e3 <- attr(D_list[[1]], "x"); lab_e3 <- attr(D_list[[1]], "labels")
w_e3 <- comp |> filter(is.finite(x_log_mg)) |> distinct(descriptor, lower_mg, upper_mg) |>
  arrange(lower_mg) |> mutate(width_mg = upper_mg - lower_mg, log_width = log(upper_mg / lower_mg))
R3 <- e3_projector_w(x_e3, 3)

closure <- comp |> group_by(patient_id, visit_code) |>
  summarise(sum_clr = sum(clr), sum_raw = sum(value), sum_repl = sum(value_r), n_parts = n(), n_zero = sum(is0), .groups = "drop")
zeros_by_part <- comp |> group_by(visit_code, descriptor, lower_mg) |>
  summarise(n = n(), n_zero_replaced = sum(is0), pct_zero = 100 * mean(is0), .groups = "drop") |> arrange(lower_mg, visit_code)

# Invariance: (i) clr centre, (ii) per-10-mg normalisation, (iii) multiplicative-replacement
# rescaling are constants per recording or per bin and cancel in the residual from a spline
# that contains an intercept. Check numerically against the saved primary residuals.
lr_mat <- function(v, transform) {
  b <- comp |> filter(is.finite(x_log_mg)) |> left_join(w_e3 |> select(descriptor, width_mg), by = "descriptor") |>
    mutate(y = transform(value_r, width_mg)) |> select(patient_id, visit_code, descriptor, y) |>
    tidyr::pivot_wider(names_from = descriptor, values_from = y)
  m0 <- b[b$visit_code == "M0", c("patient_id", lab_e3)]; mv <- b[b$visit_code == v, c("patient_id", lab_e3)]
  j <- intersect(m0$patient_id, mv$patient_id)
  as.matrix(mv[match(j, mv$patient_id), lab_e3]) - as.matrix(m0[match(j, m0$patient_id), lab_e3])
}
prim_res <- S |> filter(estimand == "E3_residual")
invariance <- purrr::map_dfr(FOLLOWUP_VISITS, \(v) {
  o <- prim_res |> filter(visit == v); o <- o$estimate[match(lab_e3, o$descriptor)]
  r_clr  <- drop(R3 %*% colMeans(D_list[[v]]))
  r_log  <- drop(R3 %*% colMeans(lr_mat(v, \(x, w) log(x))))
  r_10mg <- drop(R3 %*% colMeans(lr_mat(v, \(x, w) log(x * 10 / w))))
  tibble::tibble(visit = v, n_pairs = nrow(D_list[[v]]),
                 max_diff_clr_vs_saved = max(abs(r_clr - o)),
                 max_diff_log_minutes_vs_clr = max(abs(r_log - r_clr)),
                 max_diff_per10mg_vs_clr = max(abs(r_10mg - r_clr)),
                 clr_centre_shift = mean(colMeans(lr_mat(v, \(x, w) log(x))) - colMeans(D_list[[v]])))
})
same_basis <- all(vapply(D_list, \(D) isTRUE(all.equal(attr(D, "x"), x_e3)), logical(1)))

# **********************************************************************
# 3. C — specification and grid variants ----
# **********************************************************************
GRID_NATIVE <- c(seq(20, 120, 10), seq(140, 300, 20))
GRID_20MG   <- seq(20, 300, 20)
GRID_400    <- c(seq(20, 120, 10), seq(140, 300, 20), seq(325, 400, 25))
VARIANTS <- tibble::tribble(
  ~variant,              ~grid,        ~df, ~weighted, ~class,
  "primary (native, df 3)",   "native", 3, FALSE, "primary",
  "linear null (df 1)",       "native", 1, FALSE, "smooth-null specification",
  "df 2",                     "native", 2, FALSE, "smooth-null specification",
  "df 4",                     "native", 4, FALSE, "smooth-null specification",
  "df 5",                     "native", 5, FALSE, "smooth-null specification",
  "width-weighted (df 3)",    "native", 3, TRUE,  "smooth-null specification",
  "uniform 20-mg grid (df 3)", "g20",   3, FALSE, "resolution change (estimand changes)",
  "bounded bins to 400 mg (df 3)", "g400", 3, FALSE, "composition range change (estimand changes)"
)
grids <- list(native = GRID_NATIVE, g20 = GRID_20MG, g400 = GRID_400)
D_grid <- purrr::imap(grids, \(g, nm) {
  cg <- if (nm == "native") comp else e3_comp_grouped(dl, g)
  purrr::set_names(purrr::map(FOLLOWUP_VISITS, \(v) e3_pair_matrix(cg, v)), FOLLOWUP_VISITS)
})
stopifnot(isTRUE(all.equal(D_grid$native$M3, D_list$M3)))
# Native grid rebuilt through e3_comp_grouped() must equal e3_composition() exactly.
grouped_native_err <- max(abs(e3_pair_matrix(e3_comp_grouped(dl, GRID_NATIVE), "M6") - D_list$M6))

var_runs <- tidyr::crossing(VARIANTS, visit = FOLLOWUP_VISITS)
var_res <- purrr::pmap(var_runs, \(variant, grid, df, weighted, class, visit) {
  D <- D_grid[[grid]][[visit]]; x <- attr(D, "x")
  w <- if (weighted) { l <- attr(D, "labels"); log(as.numeric(sub(".*-", "", l)) / as.numeric(sub("-.*", "", l))) } else NULL
  R <- e3_projector_w(x, df, w)
  dec <- e3_variant_decision(D, R, B = B_BOOT, seed = SEED + match(visit, FOLLOWUP_VISITS))
  list(variant = variant, class = class, visit = visit, grid = grid, df = df, weighted = weighted, dec = dec)
})
var_decisions <- purrr::map_dfr(var_res, \(r) tibble::tibble(
  variant = r$variant, class = r$class, visit = r$visit, n_pairs = r$dec$n_pairs, n_bins = nrow(r$dec$table),
  positive_bumps = r$dec$positive_bumps, negative_bumps = r$dec$negative_bumps,
  max_abs_residual = max(abs(r$dec$table$estimate)),
  n_bins_band_excl0 = sum(r$dec$table$band_low > 0 | r$dec$table$band_high < 0),
  verdict = ifelse(nzchar(r$dec$positive_bumps) | nzchar(r$dec$negative_bumps), "bump flagged", "no localised departure"))) |>
  mutate(variant = factor(variant, levels = VARIANTS$variant), visit = factor(visit, levels = FOLLOWUP_VISITS)) |>
  arrange(variant, visit)
var_profiles <- purrr::map_dfr(var_res, \(r) r$dec$table |> mutate(variant = r$variant, class = r$class, visit = r$visit)) |>
  mutate(variant = factor(variant, levels = VARIANTS$variant), visit = factor(visit, levels = FOLLOWUP_VISITS),
         intensity = exp(x_log_mg))

# The matrix engine resamples complete pairs of each visit; the primary bootstrap resampled
# all 71 patients. Compare the primary-variant bands with the saved primary bands.
engine_vs_primary <- var_profiles |> filter(variant == "primary (native, df 3)") |>
  inner_join(prim_res |> transmute(visit = factor(as.character(visit), levels = FOLLOWUP_VISITS), descriptor,
                                   saved_est = estimate, saved_low = band_low, saved_high = band_high),
             by = c("visit", "descriptor")) |>
  group_by(visit) |>
  summarise(max_diff_estimate = max(abs(estimate - saved_est)),
            median_halfwidth_engine = median((band_high - band_low) / 2),
            median_halfwidth_saved = median((saved_high - saved_low) / 2), .groups = "drop")

# **********************************************************************
# 4. C — plasmode for the df choice ----
# **********************************************************************
bins_of <- function(lo, hi) which(as.numeric(sub("-.*", "", lab_e3)) >= lo & as.numeric(sub(".*-", "", lab_e3)) <= hi)
scen <- tibble::tribble(
  ~scenario,              ~lo, ~hi, ~size,
  "null",                  NA,  NA,  0,
  "2 bins 120-160, +0.2", 120, 160, 0.2,
  "2 bins 120-160, +0.1", 120, 160, 0.1,
  "6 bins 70-140, +0.2",   70, 140, 0.2,
  "6 bins 70-140, +0.1",   70, 140, 0.1)
scen$hit <- purrr::map2(scen$lo, scen$hi, \(a, b) if (is.na(a)) integer(0) else bins_of(a, b))
jobs <- bind_rows(
  tidyr::crossing(scen |> select(scenario, size, hit), df_true = 3L, df_test = 2:5, visit = FOLLOWUP_VISITS),
  tidyr::crossing(scen |> filter(scenario == "null") |> select(scenario, size, hit), df_true = 5L, df_test = 2:5, visit = FOLLOWUP_VISITS)
)
message(sprintf("[08] df plasmode: %d jobs, N_SIM = %d, B_SIM = %d — %s", nrow(jobs), N_SIM, B_SIM, format(Sys.time(), "%H:%M")))
if (N_WORKERS > 1L) future::plan(future::multisession, workers = N_WORKERS) else future::plan(future::sequential)
df_sim <- furrr::future_pmap_dfr(jobs, \(scenario, size, hit, df_true, df_test, visit) {
  simulate_e3_df(D_list[[visit]], x_e3, df_true, df_test, hit = hit, size = size, n_sim = N_SIM, b_sim = B_SIM,
                 seed = SEED + 1000L * match(visit, FOLLOWUP_VISITS)) |>
    mutate(scenario = scenario, df_true = df_true, df_test = df_test, visit = visit)
}, .options = furrr::furrr_options(seed = NULL, packages = c("dplyr", "tidyr", "purrr", "tibble", "splines")))
future::plan(future::sequential)
mc_ci <- function(k, n) { ci <- stats::binom.test(k, n)$conf.int; sprintf("%.0f-%.0f", 100 * ci[1], 100 * ci[2]) }
df_sim_summary <- df_sim |>
  group_by(scenario, df_true, df_test, visit) |>
  summarise(n_sim = n(), k_any = sum(any_bump), k_det = sum(detected), .groups = "drop") |>
  mutate(pct_any_bump = 100 * k_any / n_sim, any_bump_mc95 = purrr::map2_chr(k_any, n_sim, mc_ci),
         pct_detected = ifelse(scenario == "null", NA_real_, 100 * k_det / n_sim),
         visit = factor(visit, levels = FOLLOWUP_VISITS), scenario = factor(scenario, levels = scen$scenario)) |>
  arrange(df_true, scenario, df_test, visit) |> select(-k_any, -k_det)

# **********************************************************************
# 5. D — shape contrasts between visits ----
# **********************************************************************
rec <- record_info(acc) |> filter(eligible %in% TRUE, visit_code %in% MAIN_VISITS)
all4 <- rec |> group_by(patient_id) |> filter(all(MAIN_VISITS %in% visit_code)) |> distinct(patient_id) |> pull()
PAIRS <- list(c("M3", "M6"), c("M3", "M18"), c("M6", "M18"))
shape_runs <- tidyr::crossing(pair = seq_along(PAIRS), sample = c("all with M0 + both visits", "common sample (all four visits)"))
shape_res <- purrr::pmap(shape_runs, \(pair, sample) {
  a <- PAIRS[[pair]][1]; b <- PAIRS[[pair]][2]
  r <- shape_contrast(D_list[[a]], D_list[[b]], B = B_BOOT, seed = SEED + 10L * pair,
                      keep = if (grepl("common", sample)) all4 else NULL)
  r$summary <- r$summary |> mutate(contrast = paste(b, "minus", a), sample = sample)
  r$bins <- r$bins |> mutate(contrast = paste(b, "minus", a), sample = sample)
  r
})
shape_summary <- purrr::map_dfr(shape_res, "summary") |>
  group_by(sample) |>
  mutate(p_wald_boot_holm = p.adjust(p_wald_boot, "holm"), slope_p_holm = p.adjust(slope_p_boot, "holm")) |>
  ungroup() |> relocate(contrast, sample)
shape_bins <- purrr::map_dfr(shape_res, "bins") |> mutate(intensity = exp(x_log_mg))

# The bin-wise M6 - M3 contrast must equal the original C1 estimates (bounded bins).
c1_check <- shape_bins |> filter(contrast == "M6 minus M3", grepl("^all", sample)) |>
  inner_join(descriptive_objects$c1 |> transmute(descriptor, c1 = estimate), by = "descriptor") |>
  summarise(max_diff = max(abs(estimate - c1)), n = n())

# Visit-specific smooth curves (spline fitted to each visit's mean clr change; same basis).
smooth_curves <- purrr::map_dfr(FOLLOWUP_VISITS, \(v) {
  m <- colMeans(D_list[[v]])
  tibble::tibble(visit = v, descriptor = lab_e3, intensity = exp(x_e3), observed = m, smooth = m - drop(R3 %*% m),
                 n_pairs = nrow(D_list[[v]]))
}) |> mutate(visit = factor(visit, levels = FOLLOWUP_VISITS))

# Timing and season (descriptive; season is tied to visit timing by design).
timing <- purrr::map_dfr(PAIRS, \(p) {
  w <- rec |> filter(visit_code %in% c("M0", p)) |> select(patient_id, visit_code, season, months_from_m0) |>
    tidyr::pivot_wider(names_from = visit_code, values_from = c(season, months_from_m0)) |>
    filter(!is.na(.data[["season_M0"]]), !is.na(.data[[paste0("season_", p[1])]]), !is.na(.data[[paste0("season_", p[2])]]))
  sa <- w[[paste0("season_", p[1])]]; sb <- w[[paste0("season_", p[2])]]
  ta <- w[[paste0("months_from_m0_", p[1])]]; tb <- w[[paste0("months_from_m0_", p[2])]]
  tibble::tibble(contrast = paste(p[2], "minus", p[1]), n_patients = nrow(w),
                 months_a_median = median(ta), months_a_range = sprintf("%.1f-%.1f", min(ta), max(ta)),
                 months_b_median = median(tb), months_b_range = sprintf("%.1f-%.1f", min(tb), max(tb)),
                 gap_months_median = median(tb - ta), n_gap_below_1_month = sum(tb - ta < 1),
                 pct_same_season_a_b = 100 * mean(sa == sb),
                 pct_same_season_M0_a = 100 * mean(sa == w$season_M0), pct_same_season_M0_b = 100 * mean(sb == w$season_M0))
})
season_by_visit <- rec |> count(visit_code, season) |> group_by(visit_code) |> mutate(pct = 100 * n / sum(n)) |> ungroup()

# **********************************************************************
# 6. Figures ----
# **********************************************************************
lt_follow <- VISIT_LINETYPES[FOLLOWUP_VISITS]
e3c <- S |> filter(estimand == "E3_clr", is.finite(upper_mg), lower_mg > 0) |> mutate(intensity = sqrt(lower_mg * upper_mg))
e3r <- S |> filter(estimand == "E3_residual") |> mutate(intensity = sqrt(lower_mg * upper_mg))
f3a <- ggplot(e3c, aes(intensity, estimate, linetype = visit, group = visit)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, fill = "grey92", alpha = .5) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = .4) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "grey40", alpha = .14, colour = NA) +
  geom_line(data = smooth_curves, aes(intensity, smooth), colour = "grey30", linewidth = .7, linetype = "solid") +
  geom_point(size = 1.6) + geom_line(linewidth = .5) + facet_wrap(~visit, ncol = 1) +
  scale_linetype_manual(values = lt_follow) + x_axis_mg(c(20, 300)) +
  labs(x = "Bin (geometric-mean intensity, mg)", y = "Mean paired clr change vs M0",
       title = "A. Smooth redistribution",
       subtitle = "points = observed; dark curve = visit-specific smooth trend (ns, df 3)") +
  theme_myo() + theme(legend.position = "none", plot.subtitle = element_text(size = 11))
f3b <- ggplot(e3r, aes(intensity, estimate, group = visit)) +
  geom_ribbon(aes(ymin = band_low, ymax = band_high), fill = "grey40", alpha = .18, colour = NA) +
  geom_hline(yintercept = 0, colour = "grey20", linewidth = .5) +
  geom_line(aes(linetype = visit), linewidth = .8) + geom_point(size = 1.4) + facet_wrap(~visit, ncol = 1) +
  scale_linetype_manual(values = lt_follow) + x_axis_mg(c(20, 300)) +
  labs(x = "Bin (geometric-mean intensity, mg)", y = "Residual from smooth trend",
       title = "B. Localised departures",
       subtitle = "shaded = simultaneous 95% band (max-t); a departure needs the band off 0") +
  theme_myo() + theme(legend.position = "none", plot.subtitle = element_text(size = 11))
fig_C1 <- patchwork::wrap_plots(f3a, f3b, ncol = 2) +
  patchwork::plot_annotation(caption = paste(
    "A: clr change > 0 (grey area) = the bin takes a larger share of the 24-h day relative to the geometric mean of all 21 parts; < 0 = smaller share.",
    "B: residual > 0 = more gain than the smooth trend predicts at that intensity; < 0 = less.",
    sprintf("Localised departure = >= 2 adjacent bins with the band excluding 0, bounded by bins that do not, away from the grid edges, invariant to zero handling / valid days / UCTD exclusion. Cluster bootstrap B = %d.",
            primary_objects$boot_settings$B), sep = "\n"))

fig_C2 <- ggplot(var_profiles |> filter(class != "primary"), aes(intensity, estimate)) +
  geom_ribbon(data = var_profiles |> filter(class == "primary") |> select(-variant),
              aes(ymin = band_low, ymax = band_high), fill = "grey40", alpha = .18) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(data = var_profiles |> filter(class == "primary") |> select(-variant), linewidth = 1.1, colour = "black") +
  geom_line(aes(colour = variant), linewidth = .6) +
  facet_wrap(~visit, ncol = 1) + x_axis_mg(c(20, 400)) +
  scale_colour_grey(start = 0.15, end = 0.7) +
  labs(x = "Bin (geometric-mean intensity, mg)", y = "Residual from smooth null",
       colour = "Variant", title = "E3 localised-departure residuals under alternative specifications",
       subtitle = "black + shaded band = primary (native grid, df 3); thin lines = variants") +
  theme_myo() + theme(legend.position = "right")

d_c3 <- df_sim_summary |>
  mutate(value = ifelse(scenario == "null", pct_any_bump, pct_detected),
         panel = ifelse(scenario == "null", paste0("Any bump flagged, true smooth df = ", df_true), paste("Detected:", scenario))) |>
  filter(!(scenario != "null" & df_true != 3))
fig_C3 <- ggplot(d_c3, aes(df_test, value, linetype = visit, group = visit)) +
  geom_hline(data = tibble::tibble(panel = unique(d_c3$panel[grepl("^Any", d_c3$panel)]), y = 5), aes(yintercept = y),
             colour = "grey50") +
  geom_line(linewidth = .7) + geom_point(size = 1.6) + facet_wrap(~panel, ncol = 3) +
  scale_linetype_manual(values = lt_follow) + scale_x_continuous(breaks = 2:5) +
  labs(x = "df of the smooth null used in the test", y = "% of simulated datasets", linetype = "Visit",
       title = "Choice of df: false alarms and power of the bump rule",
       subtitle = sprintf("True-null plasmode; %d datasets x %d bootstrap replicates per cell; grey line = 5 %%", N_SIM, B_SIM)) +
  theme_myo() + theme(legend.position = "bottom", strip.text = element_text(size = 9))

CONTRAST_LEVELS <- c("M6 minus M3", "M18 minus M3", "M18 minus M6")
d_d1 <- shape_bins |> filter(grepl("^all", sample)) |> mutate(contrast = factor(contrast, levels = CONTRAST_LEVELS))
fig_D1a <- ggplot(smooth_curves, aes(intensity, smooth, linetype = visit, group = visit)) +
  geom_hline(yintercept = 0, colour = "grey40") + geom_line(linewidth = .9) +
  geom_point(aes(y = observed, shape = visit), size = 1.4) +
  scale_linetype_manual(values = lt_follow) + scale_shape_manual(values = VISIT_SHAPES[FOLLOWUP_VISITS]) +
  x_axis_mg(c(20, 300)) +
  labs(x = "Bin (geometric-mean intensity, mg)", y = "Mean clr change vs M0", linetype = "Visit", shape = "Visit",
       title = "A. Visit-specific smooth trends (descriptive)",
       subtitle = sprintf("pairs: %s", paste(FOLLOWUP_VISITS, purrr::map_int(D_list, nrow), sep = " n=", collapse = ", "))) +
  theme_myo() + theme(legend.position = "bottom", plot.subtitle = element_text(size = 11))
fig_D1b <- ggplot(d_d1, aes(intensity, estimate)) +
  geom_ribbon(aes(ymin = band_low, ymax = band_high), fill = "grey40", alpha = .15) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "grey40", alpha = .20) +
  geom_hline(yintercept = 0, colour = "grey30") + geom_line(linewidth = .8) + geom_point(size = 1.2) +
  facet_wrap(~contrast, ncol = 1) + x_axis_mg(c(20, 300)) +
  geom_text(data = shape_summary |> filter(grepl("^all", sample)) |> mutate(contrast = factor(contrast, levels = CONTRAST_LEVELS)),
            aes(x = 25, y = Inf, label = sprintf("n = %d; tilt diff %.2f (%.2f to %.2f); smooth 3-df p = %.2f (Holm %.2f)",
                                                 n_patients, slope_diff, slope_ci_low, slope_ci_high, p_wald_boot, p_wald_boot_holm)),
            hjust = 0, vjust = 1.4, size = 3.2) +
  labs(x = "Bin (geometric-mean intensity, mg)", y = "Difference in clr change (later - earlier)",
       title = "B. Paired differences between visits",
       subtitle = "same patients; dark = pointwise 95% CI, light = simultaneous band") +
  theme_myo() + theme(plot.subtitle = element_text(size = 11))
fig_D1 <- patchwork::wrap_plots(fig_D1a, fig_D1b, ncol = 2, widths = c(1, 1.2))

# **********************************************************************
# 7. Save ----
# **********************************************************************
save_figure(fig_C1, "08_C1_F3v2_E3_shape", width = 14, height = 11)
save_figure(fig_C2, "08_C2_E3_variants", width = 13, height = 10)
save_figure(fig_C3, "08_C3_df_simulation", width = 13, height = 8)
save_figure(fig_D1, "08_D1_shape_contrasts", width = 15, height = 10)
save_tables(list(C_invariance = invariance, C_closure_summary = closure |> summarise(max_abs_sum_clr = max(abs(sum_clr)),
                   max_abs_raw_minus_1440 = max(abs(sum_raw - 1440)), max_abs_repl_minus_1440 = max(abs(sum_repl - 1440)),
                   recordings = n(), recordings_with_zero = sum(n_zero > 0)),
                 C_zeros_by_part = zeros_by_part, C_variant_decisions = var_decisions, C_variant_profiles = var_profiles,
                 C_engine_vs_primary = engine_vs_primary, C_df_simulation = df_sim_summary,
                 D_shape_summary = shape_summary, D_shape_bins = shape_bins, D_smooth_curves = smooth_curves,
                 D_timing_season = timing, D_season_by_visit = season_by_visit), "08_REV_C_D_E3_shape")
rev_CD_objects <- list(invariance = invariance, same_basis = same_basis, grouped_native_err = grouped_native_err,
                       closure = closure |> summarise(max_abs_sum_clr = max(abs(sum_clr)), max_abs_repl_minus_1440 = max(abs(sum_repl - 1440)),
                                                      recordings = n(), recordings_with_zero = sum(n_zero > 0)),
                       zeros_by_part = zeros_by_part, variants = VARIANTS, var_decisions = var_decisions,
                       var_profiles = var_profiles, engine_vs_primary = engine_vs_primary,
                       df_sim_summary = df_sim_summary, shape_summary = shape_summary, shape_bins = shape_bins,
                       smooth_curves = smooth_curves, c1_check = c1_check, timing = timing,
                       season_by_visit = season_by_visit, n_all4 = length(all4),
                       settings = list(B = B_BOOT, seed = SEED, N_SIM = N_SIM, B_SIM = B_SIM),
                       figs = list(C1 = fig_C1, C2 = fig_C2, C3 = fig_C3, D1 = fig_D1))
save_public("rev_CD_objects", file = file.path(PATHS$rdata, "rev_CD_objects.RData"))
write_run_meta("08_REV_C_D_E3_shape_BOTH.R", list(B = B_BOOT, seed = SEED, N_SIM = N_SIM, B_SIM = B_SIM))
message("08 REV C/D finished.")
