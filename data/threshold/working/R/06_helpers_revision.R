# ======================================================================
# Script:  06_helpers_revision.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-16
# ======================================================================
#
# PURPOSE
# ---
# Functions for the A07 revision (MYO_A07_narhy_na_zmeny.txt; review in
# docs/261016_A07_critical_review.md). The original pipeline (01-06 and helpers
# 00-05) is not modified; everything new lives here. Sourced only.
#
# INPUTS
# ---
#   Analytic long table from build_analysis_long() (03_helpers_estimands.R)
#   klinicka_data_12102026.xlsx (only for block E; see locate_clinical())
#
# ANALYSES
# ---
#   A. describe_profile_norm()   — per-10-mg normalised descriptive profile
#   B. concentration_metrics()   — % with gain, positive-gain concentration,
#                                  original top-10 % share of net gain
#   C. e3_comp_grouped(), e3_projector_w(), e3_variant_decision(),
#      simulate_e3_df()          — width-aware / df sensitivity of E3
#   D. shape_contrast()          — paired M3 / M6 / M18 shape comparison
#   E. load_clinical(), fit_subtype_models() — exploratory adjusted subtype models
#
# OUTPUTS
# ---
#   Functions only.
#
# ======================================================================

# **********************************************************************
# A. Descriptive profile normalised per 10 mg ----
# **********************************************************************
# density = minutes_bin x 10 / bin_width_mg, computed per recording (master column
# minutes_per_10mg, checked against the formula in QC D11). Because the width is
# constant within a bin, every quantile of density equals the raw quantile x 10/width.
describe_profile_norm <- function(d) {
  d |> group_by(kind, visit_code, descriptor, lower_mg, upper_mg, width_mg, intensity) |>
    summarise(n = n(),
              raw_median = median(value), raw_p10 = quantile(value, .10), raw_p90 = quantile(value, .90),
              norm_median = median(density), norm_p10 = quantile(density, .10), norm_p90 = quantile(density, .90),
              norm_p25 = quantile(density, .25), norm_p75 = quantile(density, .75), norm_mean = mean(density),
              .groups = "drop") |>
    mutate(factor_10mg = 10 / width_mg)
}

# **********************************************************************
# B. Concentration of change ----
# **********************************************************************
# x = paired changes of one bin x visit (any linear unit; all shares are scale-free).
# k = ceiling(q * n) patients with the LARGEST changes, n = all complete pairs of the cell
# (selected separately within each bin and visit). Ties at the k-th value do not change
# the sum of the k largest values, so the order among tied patients is irrelevant.
#   share_net_orig : original estimand, sum(top k) / sum(all)           (NA if sum <= 0)
#   share_pos      : sum of the k largest POSITIVE changes / sum of all positive changes
#   share_pos_ref  : the same quantity if all positive changes were equal = min(k, n_pos) / n_pos
#   net_to_gross   : sum(x) / sum(|x|) in [-1, 1]; near 0 = gains and losses cancel
top_k <- function(n, q = 0.10) max(1L, ceiling(q * n))
concentration_metrics <- function(x, q = 0.10) {
  n <- length(x); k <- top_k(n, q)
  pos <- x[x > 0]; neg <- x[x < 0]
  s <- sort(x, decreasing = TRUE)
  net <- sum(x); gross <- sum(abs(x))
  tibble::tibble(
    n_pairs = n, k_top = k, top_fraction_actual = k / n,
    n_pos = length(pos), n_neg = length(neg), n_zero = sum(x == 0),
    pct_positive = 100 * length(pos) / n, pct_negative = 100 * length(neg) / n, pct_zero = 100 * mean(x == 0),
    mean_change = mean(x), median_change = median(x),
    mean_pos_part = sum(pos) / n, mean_neg_part = -sum(neg) / n,
    net_to_gross = if (gross > 0) net / gross else NA_real_,
    share_net_orig = if (net > 0) 100 * sum(s[seq_len(k)]) / net else NA_real_,
    share_pos = if (length(pos)) 100 * sum(sort(pos, decreasing = TRUE)[seq_len(min(k, length(pos)))]) / sum(pos) else NA_real_,
    share_pos_ref = if (length(pos)) 100 * min(k, length(pos)) / length(pos) else NA_real_
  )
}
# Original denominators are flagged unstable when gains and losses largely cancel.
NET_GROSS_STABLE <- 0.20

# **********************************************************************
# C. E3 width-aware / df sensitivity ----
# **********************************************************************
# Composition on an arbitrary grid of NATIVE GGIR edges (merging only; no interpolation):
# "<20" = 1440 - T20; bounded parts between consecutive `edges` (sums of native bins);
# ">=tail" = T_tail. Zero replacement identical to e3_composition().
e3_comp_grouped <- function(dl, edges = c(seq(20, 120, 10), seq(140, 300, 20)), zero = "epoch", frac = ZERO_FRAC) {
  tail_mg <- max(edges)
  stopifnot(all(edges %in% c(THRESHOLD_GRID_MG)), tail_mg %in% THRESHOLD_GRID_MG)
  b <- dl[dl$block == "bins" & dl$upper_mg <= tail_mg, c("patient_id", "visit_code", "lower_mg", "upper_mg", "value", "n_valid_days")]
  b$grp <- findInterval(b$lower_mg, edges)
  bg <- b |> group_by(patient_id, visit_code, n_valid_days, grp) |>
    summarise(lower_mg = min(lower_mg), upper_mg = max(upper_mg), value = sum(value), .groups = "drop") |>
    mutate(descriptor = paste0(lower_mg, "-", upper_mg))
  cu <- dl[dl$block == "cumulative" & dl$lower_mg %in% c(20, tail_mg), c("patient_id", "visit_code", "lower_mg", "value", "n_valid_days")]
  lo <- cu[cu$lower_mg == 20, ];  lo$value <- MINUTES_PER_DAY - lo$value
  hi <- cu[cu$lower_mg == tail_mg, ]
  comp <- bind_rows(
    lo |> transmute(patient_id, visit_code, descriptor = "<20", lower_mg = 0, upper_mg = 20, value, n_valid_days),
    bg |> select(patient_id, visit_code, descriptor, lower_mg, upper_mg, value, n_valid_days),
    hi |> transmute(patient_id, visit_code, descriptor = paste0(">=", tail_mg), lower_mg = tail_mg, upper_mg = Inf, value, n_valid_days)
  ) |>
    mutate(x_log_mg = ifelse(is.finite(upper_mg) & lower_mg > 0, log(sqrt(lower_mg * upper_mg)), NA_real_))
  dlim <- if (identical(zero, "epoch")) EPOCH_MIN / pmax(comp$n_valid_days, 1) else as.numeric(zero)
  comp |>
    mutate(dlim = dlim) |>
    group_by(patient_id, visit_code) |>
    mutate(is0 = value <= 0, repl = ifelse(is0, frac * dlim, 0),
           value_r = ifelse(is0, repl, value * (1 - sum(repl) / MINUTES_PER_DAY)),
           clr = log(value_r) - mean(log(value_r))) |>
    ungroup()
}

# Residual operator of a (weighted) natural-spline smooth null; df = 1 is a straight line.
# Weighted least squares: R = I - X (X'WX)^-1 X'W.
e3_projector_w <- function(x, df = 3, w = NULL) {
  X <- if (df == 1) cbind(1, x) else cbind(1, splines::ns(x, df = df))
  W <- diag(if (is.null(w)) rep(1, length(x)) else w / mean(w))
  diag(length(x)) - X %*% solve(t(X) %*% W %*% X, t(X) %*% W)
}

# Observed residuals, bootstrap max-t band and bump decision for one variant x visit.
e3_variant_decision <- function(D, R, B, seed, level = 0.95) {
  dec <- e3_band_decision(D, R, B, seed, level)
  lab <- attr(D, "labels")
  fmt_b <- function(ix) {
    if (!length(ix)) return("")
    r <- split(ix, cumsum(c(1, diff(ix) != 1)))
    paste(vapply(r, \(z) paste0(sub("-.*", "", lab[min(z)]), "-", sub(".*-", "", lab[max(z)]), " mg"), ""), collapse = "; ")
  }
  list(table = tibble::tibble(descriptor = lab, x_log_mg = attr(D, "x"), estimate = dec$res,
                              band_low = dec$lo, band_high = dec$hi),
       positive_bumps = fmt_b(dec$pos), negative_bumps = fmt_b(dec$neg), n_pairs = nrow(D))
}

# Plasmode for the df choice: the true-null population profile is the df_true spline fit
# of the observed mean (residual from df_true removed from every patient); an optional bump
# is added to bins `hit`; the decision is taken with the df_test smooth null.
simulate_e3_df <- function(D, x, df_true, df_test, hit = integer(0), size = 0, n_sim, b_sim, seed) {
  R_true <- e3_projector_w(x, df_true); R_test <- e3_projector_w(x, df_test)
  D0 <- sweep(D, 2, drop(R_true %*% colMeans(D)))
  if (length(hit)) D0[, hit] <- D0[, hit] + size
  n <- nrow(D0)
  purrr::map_dfr(seq_len(n_sim), \(i) {
    set.seed(seed + i)
    Di <- D0[sample.int(n, n, replace = TRUE), , drop = FALSE]
    dec <- e3_band_decision(Di, R_test, b_sim, seed + 1e6 + i)
    tibble::tibble(sim = i, any_bump = length(dec$pos) + length(dec$neg) > 0,
                   detected = length(hit) > 0 && any(dec$pos %in% hit))
  })
}

# **********************************************************************
# D. Shape contrasts between follow-up visits ----
# **********************************************************************
# Patients with M0, visit a and visit b. Delta_i = (clr_b - clr_0) - (clr_a - clr_0) per bin.
# (i) bin-wise mean difference with max-t band; (ii) difference of the linear tilt;
# (iii) global test that the smooth (ns df = 3) difference curve is flat: bootstrap Wald
# on the 3 non-intercept spline coefficients. The intercept is excluded by design (it is
# identified: it reflects exchange with the <20 / >=300 parts); the test asks whether the
# SHAPE within 20-300 mg differs. Slope p = two-sided centred-bootstrap tail probability.
# Patients need M0 (consistency with E3 / C1), although M0 cancels in Delta_i.
shape_contrast <- function(Da, Db, B, seed, keep = NULL, level = 0.95) {
  j <- intersect(rownames(Da), rownames(Db))
  if (!is.null(keep)) j <- intersect(j, keep)
  Dl <- Db[j, , drop = FALSE] - Da[j, , drop = FALSE]
  x <- attr(Da, "x"); n <- nrow(Dl)
  Xs <- cbind(1, splines::ns(x, df = 3)); Ps <- solve(crossprod(Xs), t(Xs))
  Xl <- cbind(1, x);                       Pl <- solve(crossprod(Xl), t(Xl))
  m  <- colMeans(Dl)
  beta <- drop(Ps %*% m)[-1]; slope <- drop(Pl %*% m)[2]
  set.seed(seed)
  idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = B)
  W <- t(apply(idx, 1, tabulate, nbins = n)) / n
  Mb <- W %*% Dl                                         # B x bins
  beta_b  <- t(Ps %*% t(Mb))[, -1, drop = FALSE]
  slope_b <- drop(t(Pl %*% t(Mb))[, 2])
  V <- stats::cov(beta_b); Vi <- solve(V)
  wald <- drop(t(beta) %*% Vi %*% beta)
  wald_b <- rowSums((sweep(beta_b, 2, beta) %*% Vi) * sweep(beta_b, 2, beta))
  se <- apply(Mb, 2, sd)
  tmax <- apply(abs(sweep(Mb, 2, m)) / rep(se, each = B), 1, max)
  crit <- unname(stats::quantile(tmax, level))
  a <- (1 - level) / 2
  list(
    summary = tibble::tibble(n_patients = n,
                             slope_diff = slope, slope_ci_low = unname(stats::quantile(slope_b, a)),
                             slope_ci_high = unname(stats::quantile(slope_b, 1 - a)),
                             slope_p_boot = mean(abs(slope_b - slope) >= abs(slope)),
                             wald_3df = wald, p_wald_chisq = stats::pchisq(wald, 3, lower.tail = FALSE),
                             p_wald_boot = mean(wald_b >= wald),
                             n_bins_band_excl0 = sum(m - crit * se > 0 | m + crit * se < 0)),
    bins = tibble::tibble(descriptor = attr(Da, "labels"), x_log_mg = x, estimate = m,
                          ci_low = apply(Mb, 2, stats::quantile, a), ci_high = apply(Mb, 2, stats::quantile, 1 - a),
                          band_low = m - crit * se, band_high = m + crit * se)
  )
}

# **********************************************************************
# E. Clinical covariates and exploratory subtype models ----
# **********************************************************************
CLIN_NAME <- "klinicka_data_12102026.xlsx"
locate_clinical <- function() {
  env  <- Sys.getenv("MYOII_CLINICAL_XLSX")
  cand <- c(if (nzchar(env)) env, file.path(WORKING_DIR, "data", CLIN_NAME), file.path(dirname(WORKING_DIR), CLIN_NAME))
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("Clinical file not found. Tried:\n  ", paste(cand, collapse = "\n  "),
                         "\nSet Sys.setenv(MYOII_CLINICAL_XLSX = '<path>').", call. = FALSE)
  normalizePath(hit[1], winslash = "/")
}
# Baseline (M0) clinical covariates; "x" / blank cells -> NA; IDs "IMMET_01" -> "IMMET01".
load_clinical <- function(path = locate_clinical()) {
  num <- function(v) suppressWarnings(as.numeric(ifelse(trimws(v) %in% c("", "x", "X", "NA"), NA, v)))
  readxl::read_excel(path, sheet = 1, col_types = "text") |>
    transmute(patient_id = gsub("_", "", trimws(`IMMET ID`)), visit_code = trimws(`Poradie vyšetrenia`),
              subtype_clin = trimws(`Podtyp nemoci`),
              ph_vas = num(`Physician activity VAS`), muscle_vas = num(`Muscle Disease Activity`),
              extramuscular_vas = num(`Extramuscular Global Assessment`),
              mmt8 = num(`MMT8 total`), fi2 = num(`FI-2`), haq = num(HAQ),
              disease_duration = num(`Trvání nemoci`), gc_mg_day = num(`Dávka GK (mg/den) equivalent Prednisonu`))
}

# Patient-level summaries of change (pre-specified in 04): dT100 (min/day), clr tilt, dIG.
patient_change_summaries <- function(acc, dl) {
  comp <- e3_composition(dl)
  tilt <- comp |> filter(is.finite(x_log_mg)) |>
    select(patient_id, visit_code, descriptor, x_log_mg, clr) |>
    tidyr::pivot_wider(names_from = visit_code, values_from = clr) |>
    tidyr::pivot_longer(any_of(FOLLOWUP_VISITS), names_to = "visit", values_to = "clr_fu") |>
    filter(!is.na(M0), !is.na(clr_fu)) |>
    group_by(patient_id, visit) |>
    summarise(tilt = unname(coef(lm(I(clr_fu - M0) ~ x_log_mg))[2]), .groups = "drop")
  wide_change <- function(x, nm) x |>
    tidyr::pivot_wider(names_from = visit_code, values_from = v) |>
    tidyr::pivot_longer(any_of(FOLLOWUP_VISITS), names_to = "visit", values_to = "fu") |>
    filter(!is.na(M0), !is.na(fu)) |> transmute(patient_id, visit, base = M0, chg = fu - M0) |>
    rename_with(\(z) paste0(z, "_", nm), c(base, chg))
  t100 <- dl |> filter(block == "cumulative", descriptor == "T100") |>
    transmute(patient_id, visit_code = as.character(visit_code), v = value) |> wide_change("t100")
  ig <- acc$record_level |> tibble::as_tibble() |>
    filter(as.character(visit_code) %in% MAIN_VISITS, qc_analysis_eligible %in% TRUE) |>
    transmute(patient_id = as.character(patient_id), visit_code = as.character(visit_code),
              v = .data[["part2_summary_AD_ig_gradient_ENMO_0-24hr"]]) |> wide_change("ig")
  t100 |> full_join(ig, by = c("patient_id", "visit")) |> full_join(tilt, by = c("patient_id", "visit")) |>
    mutate(visit = factor(visit, levels = FOLLOWUP_VISITS))
}

# Model set for one outcome (random-intercept LMM on change from M0, ML for LRTs).
# sev = centred baseline severity (ph_VAS per 10 mm or MMT8 per 5 points). In E_c the
# 'severity (per unit)' row is the effect at M3 (reference visit of the sev x visit term).
SUBTYPE_MODELS <- list(
  E_a_unadjusted   = "chg ~ visit * subtype",
  E_b_sev          = "chg ~ visit * subtype + sev",
  E_c_sev_by_visit = "chg ~ visit * subtype + sev * visit",
  E_e_sev_baseline = "chg ~ visit * subtype + sev + base"
)
fit_subtype_models <- function(d, models = SUBTYPE_MODELS) {
  grid <- tidyr::crossing(visit = factor(FOLLOWUP_VISITS, levels = FOLLOWUP_VISITS),
                          subtype = factor(c("DM", "IMNM", "ASyS"), levels = c("DM", "IMNM", "ASyS"))) |>
    mutate(sev = 0, base = 0)
  purrr::imap_dfr(models, \(fml, nm) {
    f <- tryCatch(suppressMessages(lme4::lmer(stats::as.formula(paste(fml, "+ (1 | patient_id)")), data = d, REML = FALSE)),
                  error = function(e) NULL)
    if (is.null(f)) return(NULL)
    X <- stats::model.matrix(stats::as.formula(paste("~", sub("^chg ~", "", fml))), grid)
    b <- lme4::fixef(f); X <- X[, names(b), drop = FALSE]
    mm <- drop(X %*% b)
    g <- grid |> mutate(mean_change = mm)
    ctr <- function(s1, s0) g |> filter(subtype == s1) |> select(visit, a = mean_change) |>
      inner_join(g |> filter(subtype == s0) |> select(visit, b = mean_change), by = "visit") |>
      transmute(visit = as.character(visit), contrast = paste(s1, "-", s0), estimate = a - b)
    cc <- bind_rows(ctr("IMNM", "DM"), ctr("ASyS", "DM"), ctr("ASyS", "IMNM"))
    cc <- bind_rows(cc, cc |> group_by(contrast) |> summarise(estimate = mean(estimate), .groups = "drop") |> mutate(visit = "mean M3-M18"))
    sev_eff <- if ("sev" %in% names(b)) tibble::tibble(visit = "all", contrast = "severity (per unit)", estimate = unname(b["sev"])) else NULL
    bind_rows(g |> transmute(visit = as.character(visit), contrast = paste("mean change", subtype), estimate = mean_change), cc, sev_eff) |>
      mutate(model = nm, n_obs = nrow(d), n_patients = dplyr::n_distinct(d$patient_id))
  })
}

# Likelihood-ratio tests (ML) for subtype x visit, global subtype and severity x visit.
subtype_lrt <- function(d) {
  fit <- function(f) suppressMessages(lme4::lmer(stats::as.formula(paste(f, "+ (1 | patient_id)")), data = d, REML = FALSE))
  lrt <- function(f0, f1, label) {
    a <- stats::anova(fit(f0), fit(f1))
    tibble::tibble(test = label, chisq = a$Chisq[2], df = a$Df[2], p = a$`Pr(>Chisq)`[2])
  }
  bind_rows(
    lrt("chg ~ visit + subtype + sev", "chg ~ visit * subtype + sev", "subtype x visit interaction"),
    lrt("chg ~ visit + sev", "chg ~ visit * subtype + sev", "any subtype effect (main + interaction)"),
    lrt("chg ~ visit * subtype + sev", "chg ~ visit * subtype + sev * visit", "severity x visit interaction"),
    lrt("chg ~ visit * subtype", "chg ~ visit * subtype + sev", "severity main effect")
  ) |> mutate(n_obs = nrow(d), n_patients = dplyr::n_distinct(d$patient_id))
}
