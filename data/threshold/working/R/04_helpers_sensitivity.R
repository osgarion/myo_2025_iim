# ======================================================================
# Script:  04_helpers_sensitivity.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-13
# ======================================================================
#
# PURPOSE
# ---
# Model functions for the pre-specified sensitivity analyses S1–S10 of plan v2.
# Model-based sensitivities report Wald 95% CI (observed data); E3-rule invariance
# checks are bootstrapped in 03_SENSITIVITY_BOTH.R. Sourced only.
#
# INPUTS
# ---
#   Stacked analytic table from build_analysis_long()
#
# ANALYSES
# ---
#   S1  complete-pair E1 / E2                     — sens_complete_pairs()
#   S3  Tweedie GLMM E2                           — est_E2_tweedie() (03_helpers)
#   S4  LMM on log(y + c), c grid                 — sens_log_c()
#   S5  MMRM + season factor / harmonics; anchor  — sens_season(), season_anchor_test()
#   S6  LMM with natural spline in months from M0 — sens_time_spline()
#   S9  delta-adjusted tipping point for M18      — sens_tipping_m18()
#   S10 M3 change by later M18 availability       — sens_dropout_m3()
#
# OUTPUTS
# ---
#   Functions only.
#
# ======================================================================

wald <- function(est, se) tibble::tibble(estimate = est, se = se, ci_low = est - 1.96 * se, ci_high = est + 1.96 * se)

# **********************************************************************
# S1. Complete pairs ----
# **********************************************************************
sens_complete_pairs <- function(dl, block = "bins") {
  d <- dl[dl$block == block, , drop = FALSE]
  base <- d[d$visit_code == "M0", c("patient_id", "descriptor", "density", "value")]
  names(base)[3:4] <- c("d0", "v0")
  d[d$visit_code != "M0", ] |> inner_join(base, by = c("patient_id", "descriptor")) |>
    group_by(descriptor, visit = as.character(visit_code)) |>
    summarise(n_pairs = n(), E1_pairs = mean(density - d0),
              E2_pairs = ifelse(mean(v0) > 0 & mean(value) > 0, log(mean(value) / mean(v0)), NA_real_),
              .groups = "drop") |>
    mutate(block = block)
}

# **********************************************************************
# S4. log(y + c) ----
# **********************************************************************
sens_log_c <- function(dl, block = "bins", cs = c(0.01, 0.1, 1)) {
  d <- dl[dl$block == block, , drop = FALSE]
  tidyr::crossing(descriptor = unique(d$descriptor), c = cs) |>
    purrr::pmap_dfr(\(descriptor, c) {
      x <- d[d$descriptor == descriptor, ]; x$y <- log(x$value + c)
      f <- tryCatch(glmmTMB::glmmTMB(y ~ visit_code + (1 | patient_id), data = x, REML = TRUE), error = function(e) NULL)
      if (is.null(f)) return(NULL)
      b <- glmmTMB::fixef(f)$cond; v <- diag(as.matrix(vcov(f)$cond))
      nm <- intersect(paste0("visit_code", FOLLOWUP_VISITS), names(b))
      bind_cols(tibble::tibble(descriptor = descriptor, c = c, visit = sub("visit_code", "", nm)),
                wald(unname(b[nm]), sqrt(unname(v[nm]))))
    }) |> mutate(block = block, scale = "log(minutes + c), geometric-mean ratio")
}

# **********************************************************************
# S5. Season ----
# **********************************************************************
sens_season <- function(dl, block = "bins", adjust = c("factor", "harmonic")) {
  adjust <- match.arg(adjust)
  d <- dl[dl$block == block & !is.na(dl$season), , drop = FALSE]
  d$m_sin <- sin(2 * pi * d$month_num / 12); d$m_cos <- cos(2 * pi * d$month_num / 12)
  ctrl <- glmmTMB::glmmTMBControl(optimizer = stats::optim, optArgs = list(method = "BFGS"))
  purrr::map_dfr(split(d, d$descriptor, drop = TRUE), \(x) {
    s <- sd(x$density); x$y <- x$density / s
    fml <- if (adjust == "factor") y ~ visit_code + season + us(visit_code + 0 | patient_id)
           else y ~ visit_code + m_sin + m_cos + us(visit_code + 0 | patient_id)
    f <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(fml, data = x, dispformula = ~0, REML = TRUE, control = ctrl)),
                  error = function(e) NULL)
    if (is.null(f) || !isTRUE(f$sdr$pdHess)) return(NULL)
    b <- glmmTMB::fixef(f)$cond; v <- diag(as.matrix(vcov(f)$cond))
    nm <- intersect(paste0("visit_code", FOLLOWUP_VISITS), names(b))
    bind_cols(tibble::tibble(descriptor = x$descriptor[1], visit = sub("visit_code", "", nm)),
              wald(unname(b[nm]) * s, sqrt(unname(v[nm])) * s))
  }) |> mutate(block = block, adjustment = adjust)
}

# Corrected anchor test (A04/A07 bug K3: `Chi Df` -> `Df`): joint 2-df LRT of annual
# harmonics, ML fits on identical rows, log1p(minutes), random intercept.
season_anchor_test <- function(dl, block, descriptor) {
  x <- dl[dl$block == block & dl$descriptor == descriptor & !is.na(dl$month_num), , drop = FALSE]
  x$y <- log1p(x$value); x$m_sin <- sin(2 * pi * x$month_num / 12); x$m_cos <- cos(2 * pi * x$month_num / 12)
  red  <- lme4::lmer(y ~ visit_code + (1 | patient_id), data = x, REML = FALSE)
  full <- lme4::lmer(y ~ visit_code + m_sin + m_cos + (1 | patient_id), data = x, REML = FALSE)
  a <- stats::anova(red, full)
  df_col <- intersect(c("Df", "Chi Df"), names(a))[1]
  tibble::tibble(block = block, descriptor = descriptor, n_obs = nrow(x), n_patients = n_distinct(x$patient_id),
                 months_observed = n_distinct(x$month_num), LR_chisq = a$Chisq[2], LR_df = a[[df_col]][2],
                 LR_p = a$`Pr(>Chisq)`[2], sin_beta = unname(lme4::fixef(full)["m_sin"]),
                 cos_beta = unname(lme4::fixef(full)["m_cos"]))
}

# **********************************************************************
# S6. Actual time since M0 ----
# **********************************************************************
# Random-intercept LMM with ns(months, df = 3) on patients with M0; contrasts are the
# predicted differences at 3, 6 and 18 months vs 0 months (Wald via the model vcov).
sens_time_spline <- function(dl, block = "bins", at = c(M3 = 3, M6 = 6, M18 = 18)) {
  d <- dl[dl$block == block & is.finite(dl$months_from_m0), , drop = FALSE]
  kn <- stats::quantile(d$months_from_m0, c(1/3, 2/3)); bk <- range(d$months_from_m0)
  basis <- function(t) splines::ns(t, knots = kn, Boundary.knots = bk)
  purrr::map_dfr(split(d, d$descriptor, drop = TRUE), \(x) {
    s <- sd(x$density); x$y <- x$density / s
    X <- basis(x$months_from_m0); colnames(X) <- paste0("ns", seq_len(ncol(X)))
    x <- cbind(x, X)
    f <- tryCatch(glmmTMB::glmmTMB(y ~ ns1 + ns2 + ns3 + (1 | patient_id), data = x, REML = TRUE), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    b <- glmmTMB::fixef(f)$cond; V <- as.matrix(vcov(f)$cond)
    L <- cbind(0, basis(at) - basis(rep(0, length(at)))[rep(1, length(at)), , drop = FALSE])
    est <- drop(L %*% b) * s; se <- sqrt(diag(L %*% V %*% t(L))) * s
    bind_cols(tibble::tibble(descriptor = x$descriptor[1], visit = names(at), months = unname(at)), wald(est, se))
  }) |> mutate(block = block)
}

# **********************************************************************
# S9. Delta-adjusted tipping point for M18 (E1) ----
# **********************************************************************
# Missing M18 among patients with M0 is imputed by its conditional mean given the
# patient's observed visits under the fitted unstructured MMRM (MAR), then shifted by
# -delta (min/day per 10 mg). The M0->M18 mean change is linear in delta, so the
# tipping point is delta* = E(0) / share_missing.
sens_tipping_m18 <- function(dl, block = "bins", deltas = c(0, 0.25, 0.5, 1, 2)) {
  d <- dl[dl$block == block, , drop = FALSE]
  ctrl <- glmmTMB::glmmTMBControl(optimizer = stats::optim, optArgs = list(method = "BFGS"))
  purrr::map_dfr(split(d, d$descriptor, drop = TRUE), \(x) {
    s <- sd(x$density); x$y <- x$density / s
    f <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(y ~ visit_code + us(visit_code + 0 | patient_id), data = x,
                  dispformula = ~0, REML = TRUE, control = ctrl)), error = function(e) NULL)
    if (is.null(f) || !isTRUE(f$sdr$pdHess)) return(NULL)
    b <- glmmTMB::fixef(f)$cond
    mu <- setNames(b[1] + c(0, b[paste0("visit_code", FOLLOWUP_VISITS)]), MAIN_VISITS)
    idx <- paste0("visit_code", MAIN_VISITS)
    S <- glmmTMB::VarCorr(f)$cond$patient_id[idx, idx]
    dimnames(S) <- list(MAIN_VISITS, MAIN_VISITS)
    w <- tidyr::pivot_wider(x[, c("patient_id", "visit_code", "y")], names_from = visit_code, values_from = y)
    w <- w[!is.na(w$M0), ]
    for (v in setdiff(MAIN_VISITS, names(w))) w[[v]] <- NA_real_
    miss18 <- is.na(w$M18)
    y18 <- w$M18
    for (i in which(miss18)) {
      o <- MAIN_VISITS[!is.na(unlist(w[i, MAIN_VISITS]))]
      yo <- unlist(w[i, o])
      y18[i] <- mu["M18"] + S["M18", o, drop = FALSE] %*% solve(S[o, o, drop = FALSE], yo - mu[o])
    }
    e0 <- mean((y18 - w$M0) * s); p <- mean(miss18)
    tibble::tibble(descriptor = x$descriptor[1], delta = deltas, n_with_m0 = nrow(w), n_m18_imputed = sum(miss18),
                   E1_M18_delta_adjusted = e0 - deltas * p, tipping_delta = if (p > 0) e0 / p else NA_real_)
  }) |> mutate(block = block)
}

# **********************************************************************
# S10. Does M3 change differ by later M18 availability? ----
# **********************************************************************
sens_dropout_m3 <- function(dl, block = "bins") {
  d <- dl[dl$block == block, , drop = FALSE]
  has18 <- unique(d$patient_id[d$visit_code == "M18"])
  base <- d[d$visit_code == "M0", c("patient_id", "descriptor", "density")]; names(base)[3] <- "d0"
  d[d$visit_code == "M3", ] |> inner_join(base, by = c("patient_id", "descriptor")) |>
    mutate(later_M18 = ifelse(patient_id %in% has18, "with M18", "without M18"), delta = density - d0) |>
    group_by(descriptor) |>
    summarise(n_with = sum(later_M18 == "with M18"), n_without = sum(later_M18 == "without M18"),
              mean_with = mean(delta[later_M18 == "with M18"]), mean_without = mean(delta[later_M18 == "without M18"]),
              diff = mean_without - mean_with,
              se = sqrt(var(delta[later_M18 == "with M18"]) / n_with + var(delta[later_M18 == "without M18"]) / n_without),
              .groups = "drop") |>
    mutate(ci_low = diff - 1.96 * se, ci_high = diff + 1.96 * se, block = block)
}
