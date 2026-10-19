# ======================================================================
# Script:  10_REV_VALIDATION_BOTH.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-16
# ======================================================================
#
# PURPOSE
# ---
# SCOPE: BOTH. Validation of the A07 revision (section 4 of MYO_A07_narhy_na_zmeny.txt).
# Re-derives key numbers by independent (base-R) code paths, checks that unchanged
# analyses reproduce the original results, verifies bootstrap bands, sample sizes and
# patient sets, influence of single patients, and agreement of figures, tables and the
# statements used in the report. Every row is PASS / WARN / FAIL / INFO.
#
# INPUTS
# ---
#   MYOII_ACC_OUTCOMES_MASTER_V1.rds, klinicka_data_12102026.xlsx
#   output/RData/{qc,primary,descriptive}_objects.RData, primary_boot_replicates.rds (originals)
#   output/RData/rev_{AB,CD,E}_objects.RData (07-09)
#
# ANALYSES
# ---
#   R1-R3  bins vs cumulative arithmetic, widths, normalisation factor, impossible values
#   R4     independent top-10 % formulas
#   R5     clr closure, zero handling, invariance
#   R6     bootstrap reproducibility and band re-derivation
#   R7     leave-one-patient-out influence on the new estimands
#   R8     sample sizes and patient sets (hashes only)
#   R9     unchanged analyses reproduce the original numbers
#   R10    figures vs tables vs report statements
#
# OUTPUTS
# ---
#   RData:  output/RData/rev_validation_objects.RData
#   Tables: output/tables/YYMMDD_10_REV_validation.xlsx
#   Logs:   output/logs/YYMMDD_10_REV_VALIDATION_BOTH_sessionInfo.txt
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
for (o in c("qc_objects", "primary_objects", "descriptive_objects", "rev_AB_objects", "rev_CD_objects", "rev_E_objects")) {
  f <- file.path(PATHS$rdata, paste0(o, ".RData"))
  if (!file.exists(f)) stop("Missing ", f, " — run 01-09 first.")
  load(f)
}
acc <- load_master()
dl  <- build_analysis_long(acc)
val <- list()
vrow <- function(id, check, status, detail) tibble::tibble(id = id, check = check, status = status, detail = detail)
pf <- function(ok, warn = FALSE) if (ok) "PASS" else if (warn) "WARN" else "FAIL"

# **********************************************************************
# 2. R1-R3 arithmetic, widths, impossible values ----
# **********************************************************************
qc <- qc_master(acc)
r1 <- qc |> filter(id %in% c("D8", "D9", "D10", "D11", "D12"))
val$R1 <- vrow("R1", "Exclusive bins reconstructed from cumulative thresholds (D8-D12)", pf(all(r1$status == "PASS")),
               paste(r1$id, r1$status, r1$detail, sep = " ", collapse = "; "))

bw <- acc$bin_long |> distinct(lower_mg, upper_mg, bin_width_mg) |> arrange(lower_mg)
exp_w <- ifelse(bw$lower_mg < 120, 10, ifelse(bw$lower_mg < 300, 20, 25))
fac_ok <- max(abs(rev_AB_objects$profile_norm$factor_10mg - 10 / rev_AB_objects$profile_norm$width_mg)) < 1e-12
val$R2 <- vrow("R2", "Bin widths 10 (20-120), 20 (120-300), 25 (300-400) mg; factor 10/width; normalised = raw x factor",
               pf(all(bw$bin_width_mg == exp_w) && all(bw$upper_mg - bw$lower_mg == bw$bin_width_mg) && fac_ok &&
                    all(rev_AB_objects$a_compare$norm_equals_scaled_raw)),
               sprintf("%d bins; widths %s; quantiles per recording = scaled raw quantiles: %s",
                       nrow(bw), paste(names(table(bw$bin_width_mg)), table(bw$bin_width_mg), sep = " mg x", collapse = ", "),
                       all(rev_AB_objects$a_compare$norm_equals_scaled_raw)))

num_cols <- function(df) df[vapply(df, is.numeric, logical(1))]
nonfinite <- function(df, allow = character(0)) {
  x <- num_cols(df); x <- x[setdiff(names(x), allow)]
  sum(vapply(x, \(v) sum(!is.finite(v)), numeric(1)))
}
cb <- rev_AB_objects$conc
imp <- c(
  profile_negative = sum(unlist(num_cols(rev_AB_objects$profile_norm |> select(starts_with("raw_"), starts_with("norm_")))) < 0),
  profile_above_1440 = sum(unlist(num_cols(rev_AB_objects$profile_norm |> select(starts_with("raw_")))) > 1440),
  pct_out_of_0_100 = sum(cb$pct_positive < 0 | cb$pct_positive > 100 | cb$share_pos < 0 | cb$share_pos > 100 + 1e-9, na.rm = TRUE),
  net_gross_out_of_range = sum(abs(cb$net_to_gross) > 1 + 1e-12, na.rm = TRUE),
  nonfinite_profile = nonfinite(rev_AB_objects$profile_norm),
  nonfinite_conc_except_documented = nonfinite(cb |> select(-starts_with("share_net_orig")), c("lower_mg", "upper_mg", "width_mg", "intensity")),
  undefined_net_share_cells = sum(is.na(cb$share_net_orig[cb$set == "bin"])),
  nonfinite_variants = nonfinite(rev_CD_objects$var_profiles),
  nonfinite_shape = nonfinite(rev_CD_objects$shape_bins) + nonfinite(rev_CD_objects$shape_summary),
  nonfinite_models = nonfinite(rev_E_objects$models |> filter(!is.na(estimate)), c("ci_low", "ci_high", "boot_se", "p_positive", "n_boot")),
  failed_model_fits = sum(is.na(rev_E_objects$models$estimate))
)
val$R3 <- vrow("R3", "Impossible / non-finite values in revised outputs",
               pf(all(imp[setdiff(names(imp), "undefined_net_share_cells")] == 0)),
               paste(names(imp), imp, sep = "=", collapse = "; "))

# **********************************************************************
# 3. R4 independent top-10 % formulas (base R from the master) ----
# **********************************************************************
el <- acc$record_level[acc$record_level$qc_analysis_eligible %in% TRUE, c("patient_id", "visit_code")]
bl <- merge(as.data.frame(acc$bin_long)[, c("patient_id", "visit_code", "lower_mg", "minutes_bin", "bin_width_mg")], el)
bl$visit_code <- as.character(bl$visit_code)
b0 <- bl[bl$visit_code == "M0", ]; names(b0)[names(b0) == "minutes_bin"] <- "m0"
bf <- merge(bl[bl$visit_code %in% FOLLOWUP_VISITS, ], b0[, c("patient_id", "lower_mg", "m0")])
bf$d <- (bf$minutes_bin - bf$m0) * 10 / bf$bin_width_mg
indep <- do.call(rbind, lapply(split(bf, list(bf$visit_code, bf$lower_mg), drop = TRUE), \(z) {
  x <- z$d; n <- length(x); k <- ceiling(0.1 * n); o <- x[order(-x)]
  p <- x[x > 0]; po <- p[order(-p)]
  data.frame(visit = z$visit_code[1], lower_mg = z$lower_mg[1], n = n,
             i_pct_pos = 100 * mean(x > 0), i_share_net = if (sum(x) > 0) 100 * sum(o[1:k]) / sum(x) else NA,
             i_share_pos = if (length(p)) 100 * sum(po[seq_len(min(k, length(p)))]) / sum(p) else NA)
}))
cb_bin <- as.data.frame(cb[cb$set == "bin", c("visit", "lower_mg", "n_pairs", "pct_positive", "share_net_orig", "share_pos")])
cb_bin$visit <- as.character(cb_bin$visit)
cmp4 <- merge(indep, cb_bin, by = c("visit", "lower_mg"))
d4 <- c(n = max(abs(cmp4$n - cmp4$n_pairs)), pct = max(abs(cmp4$i_pct_pos - cmp4$pct_positive)),
        share_pos = max(abs(cmp4$i_share_pos - cmp4$share_pos), na.rm = TRUE),
        share_net = max(abs(cmp4$i_share_net - cmp4$share_net_orig), na.rm = TRUE),
        na_pattern = sum(is.na(cmp4$i_share_net) != is.na(cmp4$share_net_orig)))
val$R4 <- vrow("R4", "Top-10 % shares re-computed independently (base R, master bin_long)", pf(all(d4 < 1e-9) && nrow(cmp4) == 69),
               paste(sprintf("%d cells", nrow(cmp4)), paste(names(d4), signif(d4, 2), sep = " max diff ", collapse = "; "), sep = "; "))
val$R4b <- vrow("R4b", "Original C3 estimand reproduced exactly by the revision code", pf(all(rev_AB_objects$b_repro$same_pct & rev_AB_objects$b_repro$same_share)),
                sprintf("%d cells; undefined (net <= 0): %s", nrow(rev_AB_objects$b_repro),
                        paste(rev_AB_objects$undefined_cells$visit, rev_AB_objects$undefined_cells$descriptor, collapse = ", ")))

# **********************************************************************
# 4. R5 clr closure, zeros, invariance ----
# **********************************************************************
cl <- rev_CD_objects$closure
inv <- rev_CD_objects$invariance
val$R5 <- vrow("R5", "clr closure (sum clr = 0), replaced composition = 1440 min, zero handling",
               pf(cl$max_abs_sum_clr < 1e-10 && cl$max_abs_repl_minus_1440 < 1e-8),
               sprintf("max |sum clr| %.1e; max |sum - 1440| %.1e; %d of %d recordings had >= 1 zero part (max share of zeros in a part %.1f%%)",
                       cl$max_abs_sum_clr, cl$max_abs_repl_minus_1440, cl$recordings_with_zero, cl$recordings,
                       max(rev_CD_objects$zeros_by_part$pct_zero)))
val$R5b <- vrow("R5b", "E3 residual invariant to clr centre and to per-10-mg normalisation; same spline basis per visit",
                pf(max(inv$max_diff_clr_vs_saved, inv$max_diff_log_minutes_vs_clr, inv$max_diff_per10mg_vs_clr) < 1e-10 &&
                     rev_CD_objects$same_basis && rev_CD_objects$grouped_native_err < 1e-12),
                sprintf("max diff vs saved %.1e, log-minutes vs clr %.1e, per-10 mg vs clr %.1e; identical basis: %s",
                        max(inv$max_diff_clr_vs_saved), max(inv$max_diff_log_minutes_vs_clr), max(inv$max_diff_per10mg_vs_clr),
                        rev_CD_objects$same_basis))

# **********************************************************************
# 5. R6 bootstrap reproducibility and bands ----
# **********************************************************************
boot <- readRDS(file.path(PATHS$rdata, "primary_boot_replicates.rds"))
S <- primary_objects$summary
band_err <- purrr::map_dbl(FOLLOWUP_VISITS, \(v) {
  o  <- S[S$estimand == "E3_residual" & S$visit == v, ]
  bb <- boot[!boot$failed & boot$estimand == "E3_residual" & boot$visit == v, ]
  M  <- tapply(bb$estimate, list(bb$replicate, bb$descriptor), identity)
  M  <- M[, o$descriptor, drop = FALSE]
  se <- apply(M, 2, sd)
  t  <- apply(abs(sweep(M, 2, o$estimate)) / rep(se, each = nrow(M)), 1, max)
  cr <- unname(quantile(t, .95))
  max(abs(o$estimate - cr * se - o$band_low), abs(o$estimate + cr * se - o$band_high))
})
val$R6 <- vrow("R6", "Primary E3 simultaneous bands re-derived from saved replicates (base R)", pf(max(band_err) < 1e-10),
               sprintf("max |difference| %.1e; B = %d", max(band_err), primary_objects$boot_settings$B))
comp <- e3_composition(dl); Dm6 <- e3_pair_matrix(comp, "M6"); R3 <- e3_projector_w(attr(Dm6, "x"), 3)
r_a <- e3_band_decision(Dm6, R3, 300, SEED); r_b <- e3_band_decision(Dm6, R3, 300, SEED)
stat_b <- function(d) d |> group_by(visit = as.character(visit), descriptor = as.character(descriptor)) |>
  summarise(estimate = concentration_metrics(delta10)$share_pos, .groups = "drop")
pairs_b <- make_pairs(dl[dl$block == "bins", ]) |> mutate(delta10 = (value - baseline) * 10 / width_mg) |> select(patient_id, visit, descriptor, delta10)
b1 <- run_cluster_bootstrap(pairs_b, stat_b, B = 24, workers = 1)
bN <- run_cluster_bootstrap(pairs_b, stat_b, B = 24, workers = min(4L, N_WORKERS))
eng <- rev_CD_objects$engine_vs_primary
val$R6b <- vrow("R6b", "Matrix engine and revision bootstrap reproducible (same seed; 1 vs N workers)",
                pf(identical(r_a, r_b) && isTRUE(all.equal(b1$estimate, bN$estimate, tolerance = 1e-12))),
                sprintf("engine identical: %s; 1 vs %d workers identical: %s", identical(r_a, r_b), min(4L, N_WORKERS),
                        isTRUE(all.equal(b1$estimate, bN$estimate, tolerance = 1e-12))))
val$R6c <- vrow("R6c", "Matrix engine (complete pairs) vs primary bootstrap (all patients): E3 residual and band half-width",
                pf(max(eng$max_diff_estimate) < 1e-10 && all(abs(eng$median_halfwidth_engine / eng$median_halfwidth_saved - 1) < 0.15), warn = TRUE),
                paste(sprintf("%s: est diff %.0e, half-width %.3f vs %.3f", eng$visit, eng$max_diff_estimate,
                              eng$median_halfwidth_engine, eng$median_halfwidth_saved), collapse = "; "))

# **********************************************************************
# 6. R7 leave-one-patient-out influence ----
# **********************************************************************
if (N_WORKERS > 1L) future::plan(future::multisession, workers = N_WORKERS) else future::plan(future::sequential)
fo <- furrr::furrr_options(seed = NULL, packages = c("dplyr", "tidyr", "purrr", "tibble", "lme4", "splines"))
# (a) positive-gain concentration for T100 / T200 aggregates, percentage points.
agg <- make_pairs(dl[dl$block == "cumulative" & dl$descriptor %in% c("T100", "T200"), ]) |> mutate(delta10 = value - baseline)
loo_b <- furrr::future_map_dfr(unique(agg$patient_id), \(p) agg |> filter(patient_id != p) |>
  group_by(visit = as.character(visit), descriptor = as.character(descriptor)) |>
  summarise(share_pos = concentration_metrics(delta10)$share_pos, pct_pos = concentration_metrics(delta10)$pct_positive, .groups = "drop"), .options = fo) |>
  inner_join(cb |> filter(set == "aggregate") |> transmute(visit = as.character(visit), descriptor, full = share_pos, full_pct = pct_positive),
             by = c("visit", "descriptor")) |>
  group_by(visit, descriptor) |> summarise(max_shift_share_pp = max(abs(share_pos - full)), max_shift_pct_pp = max(abs(pct_pos - full_pct)), .groups = "drop")
# (a2) bin-level positive-gain share in the high-intensity bins (>= 240 mg), percentage points.
hb <- make_pairs(dl[dl$block == "bins" & dl$lower_mg >= 240, ]) |> mutate(delta10 = (value - baseline) * 10 / width_mg)
loo_b2 <- furrr::future_map_dfr(unique(hb$patient_id), \(p) hb |> filter(patient_id != p) |>
  group_by(visit = as.character(visit), descriptor = as.character(descriptor)) |>
  summarise(share_pos = concentration_metrics(delta10)$share_pos, .groups = "drop"), .options = fo) |>
  inner_join(cb |> filter(set == "bin") |> transmute(visit = as.character(visit), descriptor, full = share_pos), by = c("visit", "descriptor")) |>
  group_by(visit, descriptor) |> summarise(max_shift_share_pp = max(abs(share_pos - full)), .groups = "drop")
# (b) shape contrasts: linear-tilt difference, shift in units of the bootstrap SE (approx. CI width / 3.92).
D_list <- purrr::set_names(purrr::map(FOLLOWUP_VISITS, \(v) e3_pair_matrix(comp, v)), FOLLOWUP_VISITS)
x <- attr(D_list[[1]], "x"); Pl <- solve(crossprod(cbind(1, x)), t(cbind(1, x)))
ss <- rev_CD_objects$shape_summary |> filter(grepl("^all", sample))
loo_d <- purrr::map_dfr(list(c("M3", "M6"), c("M3", "M18"), c("M6", "M18")), \(p) {
  j <- intersect(rownames(D_list[[p[1]]]), rownames(D_list[[p[2]]]))
  Dl <- D_list[[p[2]]][j, , drop = FALSE] - D_list[[p[1]]][j, , drop = FALSE]
  full <- drop(Pl %*% colMeans(Dl))[2]
  sl <- vapply(seq_along(j), \(i) drop(Pl %*% colMeans(Dl[-i, , drop = FALSE]))[2], numeric(1))
  s <- ss[ss$contrast == paste(p[2], "minus", p[1]), ]
  tibble::tibble(contrast = paste(p[2], "minus", p[1]), n = length(j), max_shift_slope = max(abs(sl - full)),
                 max_shift_in_se = max(abs(sl - full)) / ((s$slope_ci_high - s$slope_ci_low) / 3.92),
                 sign_changes = sum(sign(sl) != sign(full)))
})
# (c) adjusted subtype contrasts (E_b, ph_VAS; mean over visits), in bootstrap SE units.
clin <- load_clinical()
chg  <- patient_change_summaries(acc, dl)
info <- record_info(acc) |> distinct(patient_id, subtype) |> mutate(subtype = factor(subtype, levels = c("DM", "IMNM", "ASyS"))) |>
  left_join(clin |> filter(visit_code == "M0") |> select(patient_id, ph_vas, mmt8), by = "patient_id")
mdat <- function(o, drop_p = NULL) {
  d <- chg |> filter(!patient_id %in% drop_p) |> inner_join(info, by = "patient_id") |> filter(!is.na(subtype), !is.na(ph_vas))
  y <- if (o == "tilt") d$tilt else d[[paste0("chg_", o)]]
  d <- d |> mutate(chg = y, base = 0) |> filter(!is.na(chg))
  pm <- d |> distinct(patient_id, ph_vas)
  d |> mutate(sev = (ph_vas - mean(pm$ph_vas)) / 10, subtype = droplevels(subtype))
}
ref_e <- rev_E_objects$models |> filter(model == "E_b_sev", severity == "ph_vas", visit == "mean M3-M18", contrast %in% c("IMNM - DM", "ASyS - DM"))
loo_e <- purrr::imap_dfr(c(dT100 = "t100", tilt = "tilt", dIG = "ig"), \(o, onm) {
  ps <- unique(mdat(o)$patient_id)
  furrr::future_map_dfr(ps, \(p) fit_subtype_models(mdat(o, p), SUBTYPE_MODELS["E_b_sev"]) |>
                          filter(visit == "mean M3-M18", contrast %in% c("IMNM - DM", "ASyS - DM")) |> select(contrast, estimate), .options = fo) |>
    inner_join(ref_e |> filter(outcome == onm) |> select(contrast, full = estimate, boot_se, ci_low, ci_high), by = "contrast") |>
    group_by(contrast) |>
    summarise(max_shift_in_se = max(abs(estimate - full)) / first(boot_se),
              ci_excludes_0 = first(ci_low > 0 | ci_high < 0),
              loo_sign_changes = sum(sign(estimate) != sign(first(full))), .groups = "drop") |>
    mutate(outcome = onm)
})
future::plan(future::sequential)
val$R7 <- vrow("R7", "Leave-one-patient-out influence on the new estimands",
               pf(max(loo_d$max_shift_in_se) <= 1 && max(loo_e$max_shift_in_se) <= 1 && max(loo_b$max_shift_share_pp) < 10 &&
                    max(loo_b2$max_shift_share_pp) < 10, warn = TRUE),
               sprintf("T100/T200 positive-gain share: max shift %.1f pp; bins >= 240 mg: max %.1f pp; shape tilt difference: max %.2f SE; adjusted subtype contrasts: max %.2f SE (IMNM-DM dT100 sign changes: %d)",
                       max(loo_b$max_shift_share_pp), max(loo_b2$max_shift_share_pp), max(loo_d$max_shift_in_se), max(loo_e$max_shift_in_se),
                       loo_e$loo_sign_changes[loo_e$outcome == "dT100" & loo_e$contrast == "IMNM - DM"]))
# R7b: fit diagnostics of the observed subtype models (singular fits, convergence messages).
fitdiag <- purrr::imap_dfr(c(dT100 = "t100", tilt = "tilt", dIG = "ig"), \(o, onm) purrr::imap_dfr(SUBTYPE_MODELS[c("E_a_unadjusted", "E_b_sev", "E_c_sev_by_visit")], \(f, nm) {
  d <- mdat(o)
  m <- suppressMessages(lme4::lmer(stats::as.formula(paste(f, "+ (1 | patient_id)")), data = d, REML = FALSE))
  tibble::tibble(outcome = onm, model = nm, singular = lme4::isSingular(m),
                 conv_messages = length(m@optinfo$conv$lme4$messages),
                 icc = { v <- as.data.frame(lme4::VarCorr(m))$vcov; v[1] / sum(v) })
}))
val$R7b <- vrow("R7b", "Subtype models (ph_VAS sample): singular fits and convergence", pf(!any(fitdiag$singular) && all(fitdiag$conv_messages == 0), warn = TRUE),
                paste(sprintf("%s/%s: singular %s, ICC %.2f", fitdiag$outcome, fitdiag$model, fitdiag$singular, fitdiag$icc), collapse = "; "))

# **********************************************************************
# 7. R8 sample sizes and patient sets ----
# **********************************************************************
hash_ids <- function(ids) substr(digest_key(sort(unique(ids))), 1, 12)
orig_pairs <- qc_objects$pairs_n
rev_pairs <- cb |> filter(set == "bin") |> group_by(visit) |> summarise(n = first(n_pairs), n_unique = n_distinct(n_pairs))
e3_pairs <- purrr::map_int(D_list, nrow)
r8 <- tibble::tibble(visit = FOLLOWUP_VISITS,
                     original_pairs = orig_pairs$n_pairs[match(FOLLOWUP_VISITS, orig_pairs$visit)],
                     revision_B_pairs = rev_pairs$n[match(FOLLOWUP_VISITS, rev_pairs$visit)],
                     revision_E3_pairs = unname(e3_pairs[FOLLOWUP_VISITS]),
                     revision_variant_pairs = purrr::map_int(FOLLOWUP_VISITS, \(v) unique(rev_CD_objects$var_decisions$n_pairs[rev_CD_objects$var_decisions$visit == v])[1]))
ids_primary <- hash_ids(dl$patient_id); ids_e3 <- hash_ids(unlist(lapply(D_list, rownames)))
ids_e3_orig <- hash_ids(unlist(lapply(FOLLOWUP_VISITS, \(v) rownames(e3_pair_matrix(e3_composition(dl), v)))))
val$R8 <- vrow("R8", "Sample sizes unchanged vs the ORIGINAL 01 output (pairs per visit); variants use the same pairs; E linkage (hash is descriptive only)",
               pf(all(r8$original_pairs == r8$revision_B_pairs) && all(r8$original_pairs == r8$revision_E3_pairs) &&
                    all(r8$revision_variant_pairs == r8$revision_E3_pairs) && all(rev_pairs$n_unique == 1) && ids_e3 == ids_e3_orig &&
                    rev_E_objects$linkage$matched_in_clinical == rev_E_objects$linkage$master_patients),
               sprintf("pairs %s; patients in analysis = %d (hash %s); E: %d/%d master patients linked, subtype codes agree %d/%d; ph_VAS model n = %s",
                       paste(r8$visit, r8$original_pairs, sep = "=", collapse = ", "), n_distinct(dl$patient_id), ids_primary,
                       rev_E_objects$linkage$matched_in_clinical, rev_E_objects$linkage$master_patients,
                       rev_E_objects$linkage$subtype_code_agreement, rev_E_objects$linkage$subtype_code_compared,
                       paste(unique(rev_E_objects$models$n_patients[rev_E_objects$models$severity == "ph_vas"]), collapse = "/")))

# **********************************************************************
# 8. R9 unchanged analyses reproduce the original numbers ----
# **********************************************************************
message("[10] re-computing the primary observed statistic (MMRM per descriptor) ...")
obs_new <- estimate_all(dl)
keys <- c("estimand", "block", "visit", "descriptor")
cmp9 <- obs_new |> select(all_of(keys), new = estimate) |>
  inner_join(primary_objects$observed |> select(all_of(keys), old = estimate), by = keys)
err9 <- cmp9 |> group_by(estimand) |> summarise(max_abs_diff = max(abs(new - old), na.rm = TRUE), n = n(), .groups = "drop")
val$R9 <- vrow("R9", "Primary observed estimates (E1-E4, all blocks) reproduce the saved originals",
               pf(nrow(cmp9) == nrow(primary_objects$observed) && max(err9$max_abs_diff) < 1e-6, warn = max(err9$max_abs_diff) < 1e-3),
               paste(sprintf("%s %.1e", err9$estimand, err9$max_abs_diff), collapse = "; "))
val$R9b <- vrow("R9b", "Original descriptive outputs reproduced (raw profile 01, C3 04, C1 04)",
                pf(all(rev_AB_objects$a_compare$raw_equals_original) && all(rev_AB_objects$b_repro$same_pct & rev_AB_objects$b_repro$same_share) &&
                     rev_CD_objects$c1_check$max_diff < 1e-12),
                sprintf("raw profile %d/%d cells; C3 %d cells; C1 max diff %.1e (%d bounded bins)",
                        sum(rev_AB_objects$a_compare$raw_equals_original), nrow(rev_AB_objects$a_compare), nrow(rev_AB_objects$b_repro),
                        rev_CD_objects$c1_check$max_diff, rev_CD_objects$c1_check$n))

# **********************************************************************
# 9. R10 figures vs tables vs statements ----
# **********************************************************************
fA <- rev_AB_objects$figs$A1$data
fB <- rev_AB_objects$figs$B1$data
fE <- rev_E_objects$figs$E1$data
ms <- rev_E_objects$models
fig_ok <- c(
  A1 = max(abs(fA$norm_median - rev_AB_objects$profile_norm$norm_median)) < 1e-12,
  B1 = isTRUE(all.equal(sort(fB$y[fB$panel == levels(fB$panel)[2]]), sort(cb$share_pos[cb$set == "bin"]))),
  C1 = isTRUE(all.equal(sort(rev_CD_objects$figs$C1[[2]]$data$estimate), sort(S$estimate[S$estimand == "E3_residual"]))),
  D1 = isTRUE(all.equal(rev_CD_objects$figs$D1[[2]]$data$estimate, rev_CD_objects$shape_bins$estimate[grepl("^all", rev_CD_objects$shape_bins$sample)])),
  E1 = all(vapply(seq_len(nrow(fE)), \(i) any(abs(ms$estimate - fE$estimate[i]) < 1e-12), logical(1)))
)
val$R10 <- vrow("R10", "Figure data equal the tabulated results", pf(all(fig_ok)), paste(names(fig_ok), fig_ok, sep = "=", collapse = "; "))

# Statements used in the revised report (each must be TRUE; the report text depends on them).
vd <- rev_CD_objects$var_decisions; dfs <- rev_CD_objects$df_sim_summary; shp <- rev_CD_objects$shape_summary
lrt <- rev_E_objects$lrt
claims <- tibble::tribble(
  ~claim, ~holds,
  "A: per-10-mg medians decrease monotonically with intensity at every visit", all(rev_AB_objects$mono_check$norm_increases == 0),
  "A: the raw profile rises at 120-140 mg at every visit (width artefact)", all(grepl("120-140", rev_AB_objects$mono_check$raw_increase_at)),
  "B: original net share undefined in some low bins and > 100 % in others", any(is.na(cb$share_net_orig)) && any(cb$share_net_orig > 100, na.rm = TRUE),
  "B: positive-gain share stays within 0-100 %", all(cb$share_pos >= 0 & cb$share_pos <= 100 + 1e-9, na.rm = TRUE),
  "C: no bump under any df >= 3, weighted or alternative-grid variant",
  all(vd$verdict[!vd$variant %in% c("linear null (df 1)", "df 2")] == "no localised departure"),
  "C: df 3 type-I rate <= 5 % (MC CI lower bound) at every visit under the df-3 true null",
  all(as.numeric(sub("-.*", "", dfs$any_bump_mc95[dfs$scenario == "null" & dfs$df_true == 3 & dfs$df_test == 3])) <= 5),
  "D: no Holm-adjusted global shape difference between visits (p >= 0.05)", all(shp$p_wald_boot_holm[grepl("^all", shp$sample)] >= 0.05),
  "D: no tilt difference between visits after Holm (p >= 0.05)", all(shp$slope_p_holm[grepl("^all", shp$sample)] >= 0.05),
  "E: no subtype x visit interaction at p < 0.05 (ph_VAS-adjusted, any outcome)",
  all(lrt$p[lrt$test == "subtype x visit interaction" & lrt$severity == "ph_vas"] >= 0.05)
)
val$R10b <- vrow("R10b", "Report statements checked against the results", if (all(claims$holds)) "PASS" else "WARN",
                 paste0(sum(claims$holds), "/", nrow(claims), " hold; not holding: ",
                        paste(claims$claim[!claims$holds], collapse = " | ")))

# **********************************************************************
# 10. Save ----
# **********************************************************************
validation_table <- bind_rows(val)
print(validation_table, n = Inf, width = 220)
save_tables(list(validation = validation_table, claims = claims, R4_independent = cmp4, R7_loo_concentration = loo_b,
                 R7_loo_high_bins = loo_b2, R7b_fit_diagnostics = fitdiag, R7_loo_shape = loo_d, R7_loo_subtype = loo_e, R8_samples = r8, R9_reproduction = err9,
                 impossible_values = tibble::enframe(imp)), "10_REV_validation")
rev_validation_objects <- list(validation_table = validation_table, claims = claims, loo_b = loo_b, loo_b2 = loo_b2,
                               fitdiag = fitdiag, loo_d = loo_d,
                               loo_e = loo_e, r8 = r8, err9 = err9, imp = imp, fig_ok = fig_ok)
save_public("rev_validation_objects", file = file.path(PATHS$rdata, "rev_validation_objects.RData"))
write_run_meta("10_REV_VALIDATION_BOTH.R")
message("10 REV VALIDATION finished: ", sum(validation_table$status == "FAIL"), " FAIL, ",
        sum(validation_table$status == "WARN"), " WARN")
