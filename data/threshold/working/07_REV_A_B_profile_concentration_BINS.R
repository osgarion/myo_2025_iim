# ======================================================================
# Script:  07_REV_A_B_profile_concentration_BINS.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-16
# ======================================================================
#
# PURPOSE
# ---
# SCOPE: BINS (+ cumulative T100 / T200 aggregates in block B).
# A07 revision, blocks A and B (docs/261016_A07_critical_review.md).
#   A — descriptive spectrum profile normalised per 10 mg of bin width
#       (the original 01 figure plots raw min/day within bins of 10/20/25 mg).
#   B — heterogeneity: % of patients with a gain, concentration of POSITIVE gains in
#       the top 10 % (new), and the original top-10 % share of the NET gain (kept,
#       with its denominator problem flagged). Cluster-bootstrap CIs for the new metrics.
# Original scripts 01-06 and their outputs are not modified.
#
# INPUTS
# ---
#   MYOII_ACC_OUTCOMES_MASTER_V1.rds — via R/00_config.R
#   output/RData/qc_objects.RData, descriptive_objects.RData — originals, for comparison
#
# ANALYSES
# ---
#   1. A: per-recording normalisation, median and P10-P90 per visit x bin; original vs
#         normalised comparison; local-bump check at the 10 -> 20 mg width change
#   2. B: concentration metrics per visit x bin (+ T100, T200), bootstrap CIs
#   3. Figures
#
# OUTPUTS
# ---
#   RData:   output/RData/rev_AB_objects.RData
#   Tables:  output/tables/YYMMDD_07_REV_A_B_profile_concentration.xlsx
#   Figures: output/figures/YYMMDD_07_A1_profile_bins_per10mg.png/.tiff
#            output/figures/YYMMDD_07_A2_profile_raw_vs_per10mg.png/.tiff
#            output/figures/YYMMDD_07_B1_F6v2_heterogeneity.png/.tiff
#            output/figures/YYMMDD_07_B2_original_net_share_flagged.png/.tiff
#   Logs:    output/logs/YYMMDD_07_REV_A_B_profile_concentration_BINS_sessionInfo.txt
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
for (f in c("01_helpers_data.R", "02_helpers_bootstrap.R", "03_helpers_estimands.R", "06_helpers_revision.R"))
  source(file.path(WORKING_DIR, "R", f))
for (o in c("qc_objects", "descriptive_objects")) {
  f <- file.path(PATHS$rdata, paste0(o, ".RData"))
  if (!file.exists(f)) stop("Missing ", f, " — run 01 and 04 first.")
  load(f)
}

acc <- load_master()
dl  <- build_analysis_long(acc)
bin_levels <- descriptor_levels(dl, "bins")
WIDTH_BREAKS <- c(120, 300)   # native GGIR width changes: 10 -> 20 mg at 120, 20 -> 25 mg at 300

# **********************************************************************
# 2. A — normalised descriptive profile ----
# **********************************************************************
spec_bins <- prepare_spectrum(acc, "bins")
stopifnot(max(abs(spec_bins$density - spec_bins$value * 10 / spec_bins$width_mg)) < 1e-8)
profile_norm <- describe_profile_norm(spec_bins)

# Raw medians must equal the original 01 table exactly (same population, same recordings).
orig_prof <- qc_objects$profiles |> filter(kind == "bins") |>
  transmute(visit_code, descriptor = as.character(descriptor), orig_median = median, orig_p10 = p10, orig_p90 = p90, orig_n = n)
a_compare <- profile_norm |> mutate(descriptor = as.character(descriptor)) |>
  inner_join(orig_prof, by = c("visit_code", "descriptor")) |>
  mutate(raw_equals_original = abs(raw_median - orig_median) < 1e-9 & abs(raw_p10 - orig_p10) < 1e-9 &
           abs(raw_p90 - orig_p90) < 1e-9 & n == orig_n,
         # Quantiles are scale-equivariant, so normalising per recording or after summarising is identical.
         norm_equals_scaled_raw = abs(norm_median - raw_median * factor_10mg) < 1e-9 &
           abs(norm_p10 - raw_p10 * factor_10mg) < 1e-9 & abs(norm_p90 - raw_p90 * factor_10mg) < 1e-9)

# Local-bump check: step from the last 10-mg bin (110-120) to the first 20-mg bin (120-140)
# and from the last 20-mg bin (280-300) to the first 25-mg bin (300-325).
step_check <- profile_norm |> arrange(visit_code, lower_mg) |> group_by(visit_code) |>
  mutate(raw_step = raw_median - lag(raw_median), norm_step = norm_median - lag(norm_median)) |>
  ungroup() |> filter(lower_mg %in% WIDTH_BREAKS) |>
  select(visit_code, descriptor, width_mg, raw_median, raw_step, norm_median, norm_step)
mono_check <- profile_norm |> arrange(visit_code, lower_mg) |> group_by(visit_code) |>
  summarise(raw_increases = sum(diff(raw_median) > 0), norm_increases = sum(diff(norm_median) > 0),
            raw_increase_at = paste(descriptor[-1][diff(raw_median) > 0], collapse = ", "),
            norm_increase_at = paste(descriptor[-1][diff(norm_median) > 0], collapse = ", "), .groups = "drop")

ybreaks <- c(0, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)
fig_A1 <- ggplot(profile_norm, aes(intensity, norm_median, linetype = visit_code, group = visit_code)) +
  geom_vline(xintercept = WIDTH_BREAKS, colour = "grey60", linetype = "dotted") +
  geom_ribbon(aes(ymin = norm_p10, ymax = norm_p90), fill = "grey50", alpha = .08, colour = NA) +
  geom_line(linewidth = .8) + geom_point(aes(shape = visit_code), size = 1.6) +
  scale_visit_linetype() + scale_shape_manual(values = VISIT_SHAPES) +
  scale_y_continuous(trans = "log1p", breaks = ybreaks) + x_axis_mg() +
  labs(title = "Descriptive spectrum profile: BINS, per 10 mg",
       subtitle = "Median (line) and P10-P90 (band) of min/day per 10 mg of bin width; y on log(1+x) scale",
       x = "Exclusive bin midpoint (mg)", y = "Minutes/day per 10 mg", linetype = "Visit", shape = "Visit",
       caption = "Each recording normalised before summarising: min/day in bin x 10 / bin width (10, 20 or 25 mg). Dotted lines: bin width changes at 120 and 300 mg.") +
  theme_myo()

a_long <- profile_norm |>
  select(visit_code, intensity, raw_median, raw_p10, raw_p90, norm_median, norm_p10, norm_p90) |>
  tidyr::pivot_longer(-c(visit_code, intensity), names_to = c("scale", ".value"), names_pattern = "(raw|norm)_(.*)") |>
  mutate(scale = factor(scale, levels = c("raw", "norm"),
                        labels = c("A. Original: min/day within bin (raw)", "B. Revised: min/day per 10 mg")))
fig_A2 <- ggplot(a_long, aes(intensity, median, linetype = visit_code, group = visit_code)) +
  geom_vline(xintercept = WIDTH_BREAKS, colour = "grey60", linetype = "dotted") +
  geom_ribbon(aes(ymin = p10, ymax = p90), fill = "grey50", alpha = .08, colour = NA) +
  geom_line(linewidth = .7) + geom_point(aes(shape = visit_code), size = 1.3) +
  facet_wrap(~scale, ncol = 2) +
  scale_visit_linetype() + scale_shape_manual(values = VISIT_SHAPES) +
  scale_y_continuous(trans = "log1p", breaks = ybreaks) + x_axis_mg() +
  labs(title = "Original vs per-10-mg descriptive profile",
       subtitle = "The rise at 120-140 mg in A is created by the doubling of bin width (10 -> 20 mg)",
       x = "Exclusive bin midpoint (mg)", y = "Median and P10-P90", linetype = "Visit", shape = "Visit") +
  theme_myo() + theme(legend.position = "bottom")

# **********************************************************************
# 3. B — concentration of change ----
# **********************************************************************
pairs <- make_pairs(dl[dl$block == "bins", ] |> mutate(descriptor = factor(descriptor, levels = bin_levels))) |>
  mutate(delta10 = (value - baseline) * 10 / width_mg)
agg <- make_pairs(dl[dl$block == "cumulative" & dl$descriptor %in% c("T100", "T200"), ] |>
                    mutate(descriptor = factor(descriptor, levels = c("T100", "T200")))) |>
  mutate(delta10 = value - baseline)

conc_stat <- function(d) {
  d |> group_by(set, visit = as.character(visit), descriptor = as.character(descriptor)) |>
    group_modify(\(x, k) concentration_metrics(x$delta10)) |> ungroup() |>
    tidyr::pivot_longer(-c(set, visit, descriptor), names_to = "metric", values_to = "estimate")
}
b_data <- bind_rows(pairs |> transmute(patient_id, set = "bin", visit, descriptor = as.character(descriptor), delta10),
                    agg   |> transmute(patient_id, set = "aggregate", visit, descriptor = as.character(descriptor), delta10))
b_obs <- conc_stat(b_data)

# Original C3 values must be reproduced exactly (estimand preserved).
c3_orig <- descriptive_objects$c3 |> transmute(visit = as.character(visit), descriptor = as.character(descriptor),
                                              orig_pct_positive = pct_positive, orig_share = 100 * top10_share_of_gain)
b_repro <- b_obs |> filter(set == "bin", metric %in% c("pct_positive", "share_net_orig")) |>
  tidyr::pivot_wider(names_from = metric, values_from = estimate) |>
  inner_join(c3_orig, by = c("visit", "descriptor")) |>
  mutate(same_pct = abs(pct_positive - orig_pct_positive) < 1e-9,
         same_share = (is.na(share_net_orig) & is.na(orig_share)) | abs(share_net_orig - orig_share) < 1e-9)

boot_metrics <- c("pct_positive", "pct_negative", "share_pos", "share_net_orig", "net_to_gross", "mean_pos_part", "mean_neg_part")
b_boot <- run_cluster_bootstrap(b_data, \(d) conc_stat(d) |> filter(metric %in% boot_metrics),
                                B = B_BOOT, seed = SEED, workers = N_WORKERS,
                                checkpoint_dir = file.path(PATHS$rdata, "boot_chunks"), tag = "rev_B", chunk_size = 500L)
b_summ <- summarise_bootstrap(b_obs |> filter(metric %in% boot_metrics), b_boot, keys = c("set", "visit", "descriptor", "metric"))
conc <- b_obs |> tidyr::pivot_wider(names_from = metric, values_from = estimate) |>
  left_join(b_summ |> select(set, visit, descriptor, metric, ci_low, ci_high) |>
              tidyr::pivot_wider(names_from = metric, values_from = c(ci_low, ci_high), names_glue = "{metric}_{.value}"),
            by = c("set", "visit", "descriptor")) |>
  left_join(dl |> filter(block == "bins") |> distinct(descriptor, lower_mg, upper_mg, width_mg, intensity), by = "descriptor") |>
  mutate(visit = factor(visit, levels = FOLLOWUP_VISITS),
         net_denominator_stable = !is.na(net_to_gross) & net_to_gross >= NET_GROSS_STABLE,
         share_net_orig_status = case_when(is.na(share_net_orig) ~ "undefined (net change <= 0)",
                                           !net_denominator_stable ~ "unstable (net/gross < 0.20)",
                                           TRUE ~ "interpretable"),
         region = factor(ifelse(lower_mg < 120, "20-120 mg (10-mg bins)", "120-400 mg (20/25-mg bins)"),
                         levels = c("20-120 mg (10-mg bins)", "120-400 mg (20/25-mg bins)")))
undefined_cells <- conc |> filter(set == "bin", is.na(share_net_orig)) |> select(visit, descriptor, net_to_gross)

# **********************************************************************
# 4. B — figures ----
# **********************************************************************
cb <- conc |> filter(set == "bin")
d_b1 <- bind_rows(
  cb |> transmute(visit, intensity, region, panel = "A. Patients with a gain (%)", y = pct_positive,
                  lo = pct_positive_ci_low, hi = pct_positive_ci_high, ref = 50),
  cb |> transmute(visit, intensity, region, panel = "B. Top-10 % share of positive gains (%)", y = share_pos,
                  lo = share_pos_ci_low, hi = share_pos_ci_high, ref = share_pos_ref),
  cb |> transmute(visit, intensity, region, panel = "C. Net / gross change (cancellation index)", y = net_to_gross,
                  lo = net_to_gross_ci_low, hi = net_to_gross_ci_high, ref = 0)
) |> mutate(panel = factor(panel, levels = unique(panel)))
fig_B1 <- ggplot(d_b1, aes(intensity, y, linetype = visit, group = visit)) +
  geom_line(aes(y = ref), colour = "grey60", linewidth = .35) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey40", alpha = .10, colour = NA) +
  geom_line(linewidth = .8) + geom_point(size = 1.1) +
  facet_grid(panel ~ region, scales = "free", switch = "y") +
  scale_linetype_manual(values = VISIT_LINETYPES[FOLLOWUP_VISITS]) +
  scale_x_continuous(breaks = seq(20, 400, 20), guide = guide_axis(angle = 45)) +
  labs(x = "Exclusive bin midpoint (mg)", y = NULL, linetype = "Visit",
       title = "Is change a broad shift or a minority effect? (revised)",
       subtitle = "Complete pairs; shaded = 95% cluster-bootstrap CI (pointwise); grey line = reference",
       caption = paste("Grey reference: A 50 %; B share expected if all positive gains were equal (k / number of patients with a gain);",
                       "C 0 (gains and losses cancel).\nTop 10 % = ceiling(0.1 x pairs) patients with the largest gain, selected separately per bin and visit.")) +
  theme_myo() + theme(legend.position = "bottom", strip.placement = "outside", strip.text.y = element_text(size = 9))

fig_B2 <- ggplot(cb |> filter(!is.na(share_net_orig)),
                 aes(intensity, share_net_orig, linetype = visit, group = visit)) +
  geom_hline(yintercept = 100, colour = "grey55", linewidth = .4) +
  geom_line(linewidth = .7) + geom_point(aes(shape = share_net_orig_status), size = 2) +
  scale_y_log10() + scale_shape_manual(values = c("interpretable" = 16, "unstable (net/gross < 0.20)" = 4)) +
  scale_linetype_manual(values = VISIT_LINETYPES[FOLLOWUP_VISITS]) + x_axis_mg() +
  labs(x = "Exclusive bin midpoint (mg)", y = "Top-10 % share of NET gain (%, log scale)",
       linetype = "Visit", shape = "Denominator",
       title = "Original estimand: share of the net cohort gain (kept for reproducibility)",
       subtitle = "Values > 100 % and the missing points arise when gains and losses nearly cancel",
       caption = sprintf("Undefined (net change <= 0, not plotted): %s.",
                         paste(paste(undefined_cells$visit, undefined_cells$descriptor), collapse = ", "))) +
  theme_myo() + theme(legend.position = "bottom", legend.box = "vertical")

# **********************************************************************
# 5. Save ----
# **********************************************************************
save_figure(fig_A1, "07_A1_profile_bins_per10mg", width = 11, height = 7)
save_figure(fig_A2, "07_A2_profile_raw_vs_per10mg", width = 14, height = 7)
save_figure(fig_B1, "07_B1_F6v2_heterogeneity", width = 12, height = 11)
save_figure(fig_B2, "07_B2_original_net_share_flagged", width = 11, height = 7)
save_tables(list(A_profile_per10mg = profile_norm, A_compare_original = a_compare, A_width_steps = step_check,
                 A_monotonicity = mono_check, B_concentration = conc, B_undefined_original = undefined_cells,
                 B_reproduce_original = b_repro, bootstrap_settings = attr(b_boot, "settings")),
            "07_REV_A_B_profile_concentration")
rev_AB_objects <- list(profile_norm = profile_norm, a_compare = a_compare, step_check = step_check,
                       mono_check = mono_check, conc = conc, undefined_cells = undefined_cells, b_repro = b_repro,
                       boot_settings = attr(b_boot, "settings"),
                       figs = list(A1 = fig_A1, A2 = fig_A2, B1 = fig_B1, B2 = fig_B2))
save_public("rev_AB_objects", file = file.path(PATHS$rdata, "rev_AB_objects.RData"))
write_run_meta("07_REV_A_B_profile_concentration_BINS.R", list(B = B_BOOT, seed = SEED))
message("07 REV A/B finished: raw profile reproduced = ", all(a_compare$raw_equals_original),
        "; original C3 reproduced = ", all(b_repro$same_pct & b_repro$same_share))
