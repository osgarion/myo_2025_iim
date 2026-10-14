# ======================================================================
# Script:  03_helpers_estimands.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-13
# ======================================================================
#
# PURPOSE
# ---
# Estimators for the locked plan v2 (docs/261013_phase3_disposition_and_locked_plan.md):
#   E1  absolute change, MMRM (MAR), min/day per 10 mg (bins) / min/day above T (cumulative)
#   E2  log ratio of marginal means from the same MMRM (MAR): log((b0 + bx) / b0);
#       all-available ratio (MCAR, E2_avail) and Tweedie GLMM (E2_tweedie) are sensitivities
#   E3  paired clr change of a 21-part 24-h composition (<20, 19 bins 20–300, >=300);
#       residuals from a smooth monotone-like null (natural spline df = 3 in log intensity)
#   E4  dz = mean paired change / SD of paired change
# E3 and E4 are complete-case (M0 + Mx pairs). Sourced only.
#
# INPUTS
# ---
#   Stacked analytic table from build_analysis_long()
#
# ANALYSES
# ---
#   1. build_analysis_long()  — blocks "bins", "cumulative", "composition" (bins + >=400)
#   2. fit_mmrm_visit(), fit_tweedie_visit() — per-descriptor models with logged method
#   3. est_E1_E2 / est_E2_avail / est_E2_tweedie / est_E3 / est_E4
#   4. estimate_all()         — the statistic for observed data and every replicate
#   5. confidence_sets(), find_bumps(), decide_localisation()
#
# OUTPUTS
# ---
#   Functions only.
#
# ======================================================================

# **********************************************************************
# 1. Stacked analytic table ----
# **********************************************************************
MINUTES_PER_DAY <- 1440
EPOCH_MIN       <- 5 / 60          # one 5-s epoch in minutes (GGIR windowsizes[1] = 5 s)
E3_TAIL_MG      <- 300             # sparse tail merged into one ">=300" part for E3
ZERO_FRAC       <- 0.65            # multiplicative replacement factor (Martín-Fernández et al.)

build_analysis_long <- function(acc, visits = MAIN_VISITS) {
  bins <- prepare_spectrum(acc, "bins", visits)
  cum  <- prepare_spectrum(acc, "cumulative", visits)
  t400 <- cum |> filter(lower_mg == 400) |> select(patient_id, visit_code, t400 = value)
  keep <- c("patient_id", "visit_code", "season", "month_num", "months_from_m0", "n_valid_days", "subtype", "subtype_code")
  top  <- t400 |> inner_join(bins |> distinct(across(all_of(keep))), by = c("patient_id", "visit_code"))
  comp <- bind_rows(
    bins |> select(all_of(keep), descriptor, lower_mg, upper_mg, width_mg, intensity, value) |>
      mutate(descriptor = as.character(descriptor)),
    top |> transmute(across(all_of(keep)), descriptor = ">=400", lower_mg = 400, upper_mg = Inf,
                     width_mg = NA_real_, intensity = NA_real_, value = t400)
  ) |> mutate(density = value, block = "composition")
  out <- bind_rows(
    bins |> mutate(descriptor = as.character(descriptor), block = "bins"),
    cum  |> mutate(descriptor = as.character(descriptor), block = "cumulative"),
    comp
  ) |>
    select(block, patient_id, visit_code, descriptor, lower_mg, upper_mg, width_mg, intensity,
           value, density, all_of(setdiff(keep, c("patient_id", "visit_code"))))
  attr(out, "levels") <- list(bins = levels(bins$descriptor), cumulative = levels(cum$descriptor),
                              composition = c(levels(bins$descriptor), ">=400"))
  out
}
descriptor_levels <- function(dl, block) attr(dl, "levels")[[block]]

# **********************************************************************
# 2. Models ----
# **********************************************************************
# MMRM on y / sd(y) (numerical stability, BFGS); the same constant rescales the estimate.
# Fallback (logged in `method`): random-intercept LMM.
fit_mmrm_visit <- function(x, response = "density", structure = c("us", "ri")) {
  structure <- match.arg(structure)
  s <- sd(x[[response]])
  if (!is.finite(s) || s <= 0) return(NULL)
  x$y <- x[[response]] / s
  x$visit_code <- droplevels(x$visit_code)
  ctrl <- glmmTMB::glmmTMBControl(optimizer = stats::optim, optArgs = list(method = "BFGS"))
  ok <- function(f) !is.null(f) && isTRUE(f$sdr$pdHess) && isTRUE(f$fit$convergence == 0)
  fit <- NULL; method <- "MMRM_us"
  if (structure == "us")
    fit <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(y ~ visit_code + us(visit_code + 0 | patient_id),
                    data = x, dispformula = ~0, REML = TRUE, control = ctrl)), error = function(e) NULL)
  if (!ok(fit)) {
    fit <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(y ~ visit_code + (1 | patient_id), data = x, REML = TRUE)),
                    error = function(e) NULL)
    method <- "LMM_ri"
  }
  if (!ok(fit)) return(NULL)
  b <- glmmTMB::fixef(fit)$cond; v <- as.matrix(vcov(fit)$cond)
  nm <- intersect(paste0("visit_code", FOLLOWUP_VISITS), names(b))
  # Marginal means under MAR: mu_M0 = b0, mu_Mx = b0 + bx -> log ratio (E2, MAR).
  b0 <- unname(b["(Intercept)"]); bx <- unname(b[nm])
  tibble::tibble(visit = sub("visit_code", "", nm), estimate = bx * s,
                 se_model = unname(sqrt(diag(v)[nm])) * s, method = method,
                 log_ratio = ifelse(b0 > 0 & b0 + bx > 0, log((b0 + bx) / b0), NA_real_))
}

# Tweedie GLMM, log link, random intercept: exp(visit coefficient) is the ratio of
# marginal means (log link + random intercept only => conditional = marginal ratio).
fit_tweedie_visit <- function(x, response = "value") {
  s <- mean(x[[response]])
  if (!is.finite(s) || s <= 0) return(NULL)
  x$y <- x[[response]] / s
  x$visit_code <- droplevels(x$visit_code)
  fit <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(y ~ visit_code + (1 | patient_id), data = x,
                  family = glmmTMB::tweedie(link = "log"))), error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$sdr$pdHess) || !isTRUE(fit$fit$convergence == 0)) return(NULL)
  b <- glmmTMB::fixef(fit)$cond; v <- as.matrix(vcov(fit)$cond)
  nm <- intersect(paste0("visit_code", FOLLOWUP_VISITS), names(b))
  tibble::tibble(visit = sub("visit_code", "", nm), estimate = unname(b[nm]),
                 se_model = unname(sqrt(diag(v)[nm])), method = "Tweedie_ri")
}

per_descriptor <- function(d, fitter, ...) {
  purrr::map_dfr(split(d, d$descriptor, drop = TRUE), \(x) {
    r <- fitter(x, ...)
    if (is.null(r)) return(tibble::tibble(descriptor = x$descriptor[1], visit = FOLLOWUP_VISITS,
                                          estimate = NA_real_, se_model = NA_real_, method = "failed"))
    mutate(r, descriptor = x$descriptor[1])
  })
}

# **********************************************************************
# 3. Estimands ----
# **********************************************************************
# E1 and E2 come from the same MMRM fit: E1 = b_x (absolute), E2 = log((b0 + b_x) / b0).
# The ratio is invariant to the per-10-mg normalisation. For the ">=400" part, E2 is taken
# from the cumulative T400 fit (T400 = time >= 400 mg).
est_E1_E2 <- function(dl, block = c("bins", "cumulative"), structure = "us") {
  block <- match.arg(block)
  r <- per_descriptor(dl[dl$block == block, , drop = FALSE], fit_mmrm_visit, response = "density", structure = structure)
  if (!"log_ratio" %in% names(r)) r$log_ratio <- NA_real_
  bind_rows(
    r |> transmute(estimand = "E1", block, visit, descriptor, estimate, se_model, method),
    r |> transmute(estimand = "E2", block, visit, descriptor, estimate = log_ratio, method)
  )
}

# Sensitivity S3: Tweedie GLMM (log link, random intercept) — observed data only (slow).
est_E2_tweedie <- function(dl, block = c("composition", "cumulative")) {
  block <- match.arg(block)
  per_descriptor(dl[dl$block == block, , drop = FALSE], fit_tweedie_visit, response = "value") |>
    mutate(estimand = "E2_tweedie", block = block)
}

# Sensitivity: ratio of means of all available measurements (independence GEE; MCAR).
est_E2_avail <- function(dl, block = c("composition", "cumulative")) {
  block <- match.arg(block)
  d <- dl[dl$block == block, , drop = FALSE]
  m <- d |> group_by(descriptor, visit_code) |> summarise(mu = mean(value), .groups = "drop")
  m0 <- m |> filter(visit_code == "M0") |> select(descriptor, mu0 = mu)
  m |> filter(visit_code != "M0") |> inner_join(m0, by = "descriptor") |>
    transmute(descriptor, visit = as.character(visit_code),
              estimate = ifelse(mu0 > 0 & mu > 0, log(mu / mu0), NA_real_),
              estimand = "E2_avail", block = if (block == "composition") "bins" else block)
}

# 21-part composition: "<20" = 1440 - T20, 19 bins 20–300, ">=300" = T300.
# Zero replacement: multiplicative, value = ZERO_FRAC x detection limit; the detection limit
# is one 5-s epoch per recording ("epoch") or a fixed min/day value (sensitivity grid).
e3_composition <- function(dl, zero = "epoch", frac = ZERO_FRAC) {
  b <- dl[dl$block == "bins" & dl$upper_mg <= E3_TAIL_MG, c("patient_id", "visit_code", "descriptor", "lower_mg", "upper_mg", "value", "n_valid_days")]
  cu <- dl[dl$block == "cumulative" & dl$lower_mg %in% c(20, E3_TAIL_MG), c("patient_id", "visit_code", "lower_mg", "value", "n_valid_days")]
  lo <- cu[cu$lower_mg == 20, ];  lo$value <- MINUTES_PER_DAY - lo$value
  hi <- cu[cu$lower_mg == E3_TAIL_MG, ]
  comp <- bind_rows(
    lo |> transmute(patient_id, visit_code, descriptor = "<20", lower_mg = 0, upper_mg = 20, value, n_valid_days),
    b |> mutate(descriptor = as.character(descriptor)),
    hi |> transmute(patient_id, visit_code, descriptor = paste0(">=", E3_TAIL_MG), lower_mg = E3_TAIL_MG, upper_mg = Inf, value, n_valid_days)
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

# E3: per part mean paired clr change; for the 19 bounded bins, residual from a natural
# spline (df = 3) in log geometric-mean intensity — the smooth null. The linear slope
# (IG-like tilt) is reported as a descriptive "is IG enough?" check.
est_E3 <- function(dl, zero = "epoch", frac = ZERO_FRAC) {
  comp <- e3_composition(dl, zero, frac)
  base <- comp |> filter(visit_code == "M0") |> select(patient_id, descriptor, clr0 = clr)
  m <- comp |> filter(visit_code != "M0") |> inner_join(base, by = c("patient_id", "descriptor")) |>
    group_by(visit = as.character(visit_code), descriptor, lower_mg, x_log_mg) |>
    summarise(estimate = mean(clr - clr0), .groups = "drop")
  purrr::map_dfr(split(m, m$visit), \(x) {
    x <- x[order(x$lower_mg), ]
    bnd <- is.finite(x$x_log_mg)
    xx <- x$x_log_mg[bnd]; yy <- x$estimate[bnd]
    fit_s <- stats::lm(yy ~ splines::ns(xx, df = 3))
    fit_l <- stats::lm(yy ~ xx)
    bind_rows(
      tibble::tibble(visit = x$visit, descriptor = x$descriptor, estimate = x$estimate, estimand = "E3_clr"),
      tibble::tibble(visit = x$visit[bnd], descriptor = x$descriptor[bnd], estimate = unname(stats::residuals(fit_s)), estimand = "E3_residual"),
      tibble::tibble(visit = x$visit[bnd], descriptor = x$descriptor[bnd], estimate = unname(stats::residuals(fit_l)), estimand = "E3_residual_linear"),
      tibble::tibble(visit = x$visit[1], descriptor = "slope_per_log_mg", estimate = unname(stats::coef(fit_l)[2]), estimand = "E3_slope")
    )
  }) |> mutate(block = "composition")
}

est_E4 <- function(dl) {
  d <- dl[dl$block == "bins", , drop = FALSE]
  base <- d[d$visit_code == "M0", c("patient_id", "descriptor", "value")]
  names(base)[3] <- "v0"
  d[d$visit_code != "M0", ] |> inner_join(base, by = c("patient_id", "descriptor")) |>
    group_by(descriptor, visit = as.character(visit_code)) |>
    summarise(estimate = mean(value - v0) / sd(value - v0), .groups = "drop") |>
    mutate(estimand = "E4", block = "bins")
}

# **********************************************************************
# 4. Full statistic ----
# **********************************************************************
estimate_all <- function(dl, structure = "us") {
  bins <- est_E1_E2(dl, "bins", structure)
  cum  <- est_E1_E2(dl, "cumulative", structure)
  top  <- cum |> filter(estimand == "E2", descriptor == "T400") |> mutate(block = "bins", descriptor = ">=400")
  bind_rows(
    bins, top, cum,
    est_E2_avail(dl, "composition"), est_E2_avail(dl, "cumulative"),
    est_E3(dl), est_E4(dl)
  ) |> select(estimand, block, visit, descriptor, estimate, any_of(c("se_model", "method")))
}

# **********************************************************************
# 5. Region inference ----
# **********************************************************************
# Inner set: simultaneous lower band > delta (region certainly includes);
# outer set: simultaneous upper band > delta (region plausibly extends to).
confidence_sets <- function(summ, deltas) {
  tidyr::crossing(summ, delta = deltas) |>
    mutate(inner = band_low > delta, outer = band_high > delta)
}

# A "bump" = run of >= min_len flagged bins, flanked on both sides by unflagged bins and
# not touching the first or last bin of the grid (edges cannot be localised).
find_bumps <- function(flag, labels, min_len = 2L) {
  flag[is.na(flag)] <- FALSE
  r <- rle(flag); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  keep <- r$values & r$lengths >= min_len & starts > 1 & ends < length(flag)
  if (!any(keep)) return(character(0))
  paste0(sub("-.*", "", labels[starts[keep]]), "-", sub(".*-", "", labels[ends[keep]]), " mg")
}

# Locked decision rule (plan v2 §3.4). Invariance across zero replacement, S7 and S8 is
# checked in 03_SENSITIVITY_BOTH.R.
decide_localisation <- function(e3_resid, e3_slope, bin_levels) {
  purrr::map_dfr(FOLLOWUP_VISITS, \(v) {
    x <- e3_resid |> filter(visit == v) |> mutate(descriptor = factor(descriptor, levels = bin_levels)) |>
      filter(!is.na(descriptor)) |> arrange(descriptor)
    pos <- find_bumps(x$band_low > 0, as.character(x$descriptor))
    neg <- find_bumps(x$band_high < 0, as.character(x$descriptor))
    sl  <- e3_slope |> filter(visit == v)
    slope_nonzero <- isTRUE(sl$ci_low > 0 | sl$ci_high < 0)
    verdict <- if (length(pos) || length(neg)) "candidate localised departure (confirm invariance in 03_SENSITIVITY)"
               else if (slope_nonzero) "no evidence of a localised departure; smooth redistribution towards higher intensities"
               else "no evidence of a localised departure; shape change not identifiable"
    tibble::tibble(visit = v, positive_bumps = paste(pos, collapse = "; "), negative_bumps = paste(neg, collapse = "; "),
                   slope = sl$estimate, slope_ci_low = sl$ci_low, slope_ci_high = sl$ci_high, verdict = verdict)
  })
}
