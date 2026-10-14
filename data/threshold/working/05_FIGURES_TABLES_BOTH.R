# ======================================================================
# Script:  05_FIGURES_TABLES_BOTH.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-13
# ======================================================================
#
# PURPOSE
# ---
# SCOPE: BOTH. Builds publication figures (B/W printable; M0 dot-dash, M3 dotted,
# M6 dashed, M18 solid) and presentation tables from the saved analysis objects.
# No model is fitted here. Also writes the figure legends.
#
# INPUTS
# ---
#   output/RData/qc_objects.RData, primary_objects.RData,
#   sensitivity_objects.RData (optional), descriptive_objects.RData (optional)
#
# ANALYSES
# ---
#   F1 E1 bins (min/day per 10 mg)        F1c E1 cumulative (min/day above T)
#   F2 E2 bins (ratio of marginal means)  F2c E2 cumulative
#   F3 E3 shape: clr change + smooth null; residuals with simultaneous band
#   F4 three dimensions: change, precision, coverage per visit
#   F5 confidence sets (inner / outer) for E1 and E2
#   F6 C3 heterogeneity (% patients with gain, top-10 % share)
#   F7 subtype profiles (descriptive)
#
# OUTPUTS
# ---
#   Figures: output/figures/YYMMDD_05_F*.png/.tiff
#   RData:   output/RData/figure_objects.RData (ggplot objects for the report)
#   Text:    output/YYMMDD_05_figure_legends.md
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
for (o in c("qc_objects", "primary_objects", "sensitivity_objects", "descriptive_objects")) {
  f <- file.path(PATHS$rdata, paste0(o, ".RData"))
  if (file.exists(f)) load(f) else if (o %in% c("qc_objects", "primary_objects")) stop("Missing ", f)
}
S <- primary_objects$summary
B_used <- primary_objects$boot_settings$B

lt_follow <- VISIT_LINETYPES[FOLLOWUP_VISITS]
ttl <- function(t, s = NULL) labs(title = t, subtitle = s)
caption_n <- function() {
  n <- qc_objects$pairs_n
  sprintf("Pairs with M0: %s. Cluster bootstrap B = %d.", paste(n$visit, n$n_pairs, sep = " n=", collapse = ", "), B_used)
}

# **********************************************************************
# 2. Profile plot helper ----
# **********************************************************************
plot_change_profile <- function(dd, ylab, title, subtitle, log_y = FALSE, hline = 0, x_lim = c(20, 400)) {
  p <- ggplot(dd, aes(intensity, estimate, linetype = visit, group = visit)) +
    geom_hline(yintercept = hline, colour = "grey40", linewidth = .4) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "grey40", alpha = .12, colour = NA) +
    geom_line(aes(y = band_low), linewidth = .3, colour = "grey30") +
    geom_line(aes(y = band_high), linewidth = .3, colour = "grey30") +
    geom_line(linewidth = .9) + geom_point(size = 1.4) +
    facet_wrap(~visit, ncol = 1) +
    scale_linetype_manual(values = lt_follow, drop = TRUE) +
    x_axis_mg(x_lim) + ttl(title, subtitle) +
    labs(y = ylab, linetype = "Visit", caption = caption_n()) +
    theme_myo() + theme(legend.position = "none")
  if (log_y) p <- p + scale_y_log10()
  p
}

# **********************************************************************
# 3. F1 / F1c — E1 ----
# **********************************************************************
f1 <- plot_change_profile(S |> filter(estimand == "E1", block == "bins"),
  "Change vs M0 (min/day per 10 mg)", "Absolute change across exclusive bins (E1)",
  "MMRM; shaded = pointwise 95% CI, thin lines = simultaneous 95% band") +
  labs(x = "Exclusive bin midpoint (mg)")
f1c <- plot_change_profile(S |> filter(estimand == "E1", block == "cumulative"),
  "Change vs M0 (min/day above T)", "Absolute change in cumulative time above threshold",
  "MMRM; overlapping exposures; shaded = pointwise 95% CI, thin lines = simultaneous band") +
  labs(x = "Cumulative threshold T (mg)")

# **********************************************************************
# 4. F2 / F2c — E2 (ratio scale; >=400 plotted at 412.5 mg, grid edge marked) ----
# **********************************************************************
e2 <- S |> filter(estimand == "E2", block == "bins") |>
  mutate(intensity = ifelse(descriptor == ">=400", 412.5, intensity),
         across(c(estimate, ci_low, ci_high, band_low, band_high), exp))
f2 <- plot_change_profile(e2, "Ratio of means, follow-up / M0", "Relative change across exclusive bins (E2)",
  "MMRM marginal means; point at 412.5 mg = time >= 400 mg (grid edge, dotted line)",
  log_y = TRUE, hline = 1, x_lim = c(20, 425)) +
  geom_vline(xintercept = 400, linetype = "dotted") + labs(x = "Exclusive bin midpoint (mg)")
e2c <- S |> filter(estimand == "E2", block == "cumulative") |>
  mutate(across(c(estimate, ci_low, ci_high, band_low, band_high), exp))
f2c <- plot_change_profile(e2c, "Ratio of means, follow-up / M0", "Relative change in cumulative time above threshold",
  "MMRM marginal means; overlapping exposures", log_y = TRUE, hline = 1) + labs(x = "Cumulative threshold T (mg)")

# **********************************************************************
# 5. F3 — E3 shape ----
# **********************************************************************
e3c <- S |> filter(estimand == "E3_clr", is.finite(upper_mg), lower_mg > 0) |>
  mutate(intensity = sqrt(lower_mg * upper_mg))
e3r <- S |> filter(estimand == "E3_residual") |> mutate(intensity = sqrt(lower_mg * upper_mg))
f3a <- ggplot(e3c, aes(intensity, estimate, linetype = visit, group = visit)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = .4) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "grey40", alpha = .12, colour = NA) +
  geom_smooth(method = "lm", formula = y ~ splines::ns(log(x), df = 3), se = FALSE, colour = "grey50", linewidth = .5) +
  geom_line(linewidth = .9) + geom_point(size = 1.4) + facet_wrap(~visit, ncol = 1) +
  scale_linetype_manual(values = lt_follow) + x_axis_mg(c(20, 300)) +
  labs(x = "Bin (geometric-mean intensity, mg)", y = "Mean paired clr change",
       title = "A. Shape of change (E3)", subtitle = "grey = smooth null (natural spline, df = 3)") +
  theme_myo() + theme(legend.position = "none")
f3b <- ggplot(e3r, aes(intensity, estimate, linetype = visit, group = visit)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = .4) +
  geom_ribbon(aes(ymin = band_low, ymax = band_high), fill = "grey40", alpha = .15, colour = NA) +
  geom_line(linewidth = .9) + geom_point(size = 1.4) + facet_wrap(~visit, ncol = 1) +
  scale_linetype_manual(values = lt_follow) + x_axis_mg(c(20, 300)) +
  labs(x = "Bin (geometric-mean intensity, mg)", y = "Residual from smooth null",
       title = "B. Localised departure?", subtitle = "shaded = simultaneous 95% band (max-t)") +
  theme_myo() + theme(legend.position = "none")
f3 <- patchwork::wrap_plots(f3a, f3b, ncol = 2) + patchwork::plot_annotation(caption = caption_n())

# **********************************************************************
# 6. F4 — three dimensions (change, precision, coverage per visit) ----
# **********************************************************************
e1b <- S |> filter(estimand == "E1", block == "bins")
cov <- qc_objects$coverage |> filter(kind == "bins", limit == COVERAGE_PRIMARY)
d4 <- bind_rows(
  e1b |> transmute(visit = as.character(visit), intensity, panel = "A. Change vs M0 (min/day per 10 mg)", y = estimate, lo = ci_low, hi = ci_high),
  e1b |> transmute(visit = as.character(visit), intensity, panel = "B. Bootstrap SE of change", y = boot_se, lo = NA, hi = NA),
  cov |> transmute(visit = as.character(visit_code), intensity, panel = sprintf("C. Recordings with > %.2f min/day per 10 mg (%%)", COVERAGE_PRIMARY), y = pct_above, lo = NA, hi = NA)
) |> mutate(visit = factor(visit, levels = MAIN_VISITS))
f4 <- ggplot(d4, aes(intensity, y, linetype = visit, group = visit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey40", alpha = .10, colour = NA, na.rm = TRUE) +
  geom_line(linewidth = .8) + facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_visit_linetype() + x_axis_mg() +
  labs(x = "Exclusive bin midpoint (mg)", y = NULL, linetype = "Visit",
       title = "Three dimensions of longitudinal information",
       subtitle = "Coverage = share of all recordings at that visit (denominator: recordings, not pairs)",
       caption = caption_n()) +
  theme_myo() + theme(legend.position = "bottom")

# **********************************************************************
# 7. F5 — confidence sets ----
# **********************************************************************
set_tiles <- S |> filter(estimand %in% c("E1", "E2"), block == "bins") |>
  mutate(state = case_when(band_low > 0 ~ "inner set (certainly > 0)",
                           band_high > 0 ~ "outer set only (plausibly > 0)",
                           TRUE ~ "not in set"),
         xmin = lower_mg, xmax = ifelse(descriptor == ">=400", 425, upper_mg),
         y = as.numeric(factor(visit, levels = FOLLOWUP_VISITS)),
         estimand = recode(estimand, E1 = "E1 absolute", E2 = "E2 relative"))
f5 <- ggplot(set_tiles, aes(xmin = xmin, xmax = xmax, ymin = y - .4, ymax = y + .4, fill = state)) +
  geom_rect(colour = "white", linewidth = .3) +
  scale_y_continuous(breaks = seq_along(FOLLOWUP_VISITS), labels = FOLLOWUP_VISITS) +
  facet_wrap(~estimand, ncol = 1) +
  scale_fill_manual(values = c("inner set (certainly > 0)" = "grey15", "outer set only (plausibly > 0)" = "grey65", "not in set" = "white")) +
  x_axis_mg(c(15, 425)) + geom_vline(xintercept = 400, linetype = "dotted") +
  labs(x = "Exclusive bin (mg); 400-425 = time >= 400 mg", y = NULL, fill = NULL,
       title = "Where does change lie?", subtitle = "Inner / outer simultaneous confidence sets for {effect > 0}") +
  theme_myo() + theme(legend.position = "bottom")

# **********************************************************************
# 8. F6 — heterogeneity (C3) ----
# **********************************************************************
f6 <- NULL
if (exists("descriptive_objects")) {
  d6 <- descriptive_objects$c3 |>
    tidyr::pivot_longer(c(pct_positive, top10_share_of_gain), names_to = "metric", values_to = "y") |>
    mutate(y = ifelse(metric == "top10_share_of_gain", 100 * y, y),
           metric = recode(metric, pct_positive = "A. Patients with a gain (%)", top10_share_of_gain = "B. Share of total gain from top-10 % patients (%)"))
  f6 <- ggplot(d6, aes(intensity, y, linetype = visit, group = visit)) +
    geom_line(linewidth = .8) + geom_point(size = 1.2) + facet_wrap(~metric, ncol = 1, scales = "free_y") +
    scale_linetype_manual(values = lt_follow) + x_axis_mg() +
    labs(x = "Exclusive bin midpoint (mg)", y = NULL, linetype = "Visit",
         title = "Is change a broad shift or a minority effect?", subtitle = "Complete pairs; descriptive") +
    theme_myo() + theme(legend.position = "bottom")
}

# **********************************************************************
## F7 — subtype profiles (descriptive; full spectrum, no single-bin selection) ----
# **********************************************************************
f7 <- NULL
if (exists("descriptive_objects")) {
  f7 <- ggplot(descriptive_objects$subtype_profile,
               aes(intensity, ratio_of_means, linetype = visit, group = visit)) +
    geom_hline(yintercept = 1, colour = "grey40", linewidth = .4) +
    geom_line(linewidth = .8) + geom_point(size = 1.1) +
    facet_wrap(~subtype, ncol = 1) + scale_y_log10() +
    scale_linetype_manual(values = lt_follow) + x_axis_mg() +
    labs(x = "Exclusive bin midpoint (mg)", y = "Ratio of means, follow-up / M0 (complete pairs)", linetype = "Visit",
         title = "Relative change by subtype (descriptive)",
         subtitle = "DM (codes 0, 4), IMNM (2, 3), ASyS (5); PM and UCTD/SSc excluded; no CIs, 12-18 pairs per cell") +
    theme_myo() + theme(legend.position = "bottom")
}

# **********************************************************************
# 9. Save ----
# **********************************************************************
figs <- list(F1_E1_bins = list(f1, 11, 10), F1c_E1_cumulative = list(f1c, 11, 10),
             F2_E2_bins = list(f2, 11, 10), F2c_E2_cumulative = list(f2c, 11, 10),
             F3_E3_shape = list(f3, 14, 10), F4_three_dimensions = list(f4, 11, 11),
             F5_confidence_sets = list(f5, 12, 6))
if (!is.null(f6)) figs$F6_heterogeneity <- list(f6, 11, 8)
if (!is.null(f7)) figs$F7_subtype_profiles <- list(f7, 11, 10)
purrr::iwalk(figs, \(x, nm) save_figure(x[[1]], paste0("05_", nm), width = x[[2]], height = x[[3]]))
figure_objects <- purrr::map(figs, 1)
save_public("figure_objects", file = file.path(PATHS$rdata, "figure_objects.RData"))

legends <- c(
  "# Figure legends (generated by 05_FIGURES_TABLES_BOTH.R)", "",
  "**Population (all figures):** QC-eligible GENEActiv recordings (non-dominant wrist, GGIR) of patients with idiopathic inflammatory myopathies; main visits M0, M3, M6, M18. Line types: M0 dot-dash, M3 dotted, M6 dashed, M18 solid. Intensity in mg ENMO. 'MVPA' in GGIR variable names is legacy naming; no validated intensity categories are implied.", "",
  sprintf("**F1 — Absolute change across exclusive bins (E1).** Mean change from M0 in time spent within each non-overlapping ENMO bin, expressed as min/day per 10 mg of bin width (bins are 10, 20 or 25 mg wide). Estimates from per-bin mixed models with unstructured visit covariance (MMRM; all available measurements, valid under missing-at-random). Shaded: pointwise 95%% percentile intervals; thin lines: simultaneous 95%% band over all 23 bins of a visit (max-t). Uncertainty from a patient-level cluster bootstrap (B = %d; patients resampled with all their visits).", B_used), "",
  "**F1c — Absolute change in cumulative time above threshold.** As F1 for minutes/day above each cumulative threshold T. Thresholds overlap (T100 contains all time above 100 mg), so adjacent estimates are not independent and the profile is not a localisation.", "",
  "**F2 — Relative change across exclusive bins (E2).** Ratio of marginal mean minutes/day at follow-up vs M0 from the same MMRM (log scale). The point at 412.5 mg is time >= 400 mg; the dotted line marks the upper edge of the threshold grid — any upper limit of a relative-change region is censored by the grid.", "",
  "**F2c — Relative change in cumulative time above threshold.** As F2 for cumulative thresholds.", "",
  "**F3 — Shape of change (E3).** (A) Mean paired change in centred log-ratio (clr) of a 21-part 24-h composition (<20 mg, 19 bins 20–300 mg, >=300 mg), complete pairs M0–Mx; zeros replaced multiplicatively (0.65 x one 5-s epoch per recording). Grey: smooth null (natural spline, df = 3, in log intensity). (B) Residuals from the smooth null with simultaneous 95% band. A localised departure requires >= 2 contiguous bins excluding zero, bounded on both sides and away from the grid edges, invariant to zero replacement, >= 5 valid days and exclusion of the UCTD/SSc patient.", "",
  sprintf("**F4 — Three dimensions of longitudinal information.** (A) E1 with pointwise 95%% CI; (B) bootstrap standard error of E1 (group-level precision; not test-retest reliability); (C) percentage of all recordings at each visit with > %.2f min/day per 10 mg in the bin (descriptive, not clinically validated; values below are valid low activity, not missing data).", COVERAGE_PRIMARY), "",
  "**F5 — Confidence sets.** Bins where the simultaneous lower bound exceeds zero (inner set: change certainly > 0) and bins where only the upper bound does (outer set: change plausibly > 0), for absolute (E1) and relative (E2) change.", "",
  "**F6 — Heterogeneity of change (C3).** Percentage of patients with a gain and the share of the total gain contributed by the 10% of patients with the largest gains, per bin; complete pairs; descriptive.", "",
  "**F7 — Relative change by subtype.** Ratio of mean minutes/day (follow-up / M0) per exclusive bin in complete pairs, separately for DM, IMNM and ASyS (PM and UCTD/SSc excluded). Descriptive only: no confidence intervals and no tests; 12-18 pairs per subtype and visit."
)
writeLines(legends, file.path(WORKING_DIR, "output", paste0(date_prefix(), "_05_figure_legends.md")))
message("05 FIGURES finished: ", length(figs), " figures.")
