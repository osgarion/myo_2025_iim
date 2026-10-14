# ======================================================================
# Script:  05_helpers_e3fast.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-14
# ======================================================================
#
# PURPOSE
# ---
# Matrix implementation of the E3 shape estimand for large simulation studies
# (round-3 review, items 5–6). E3 residuals are a linear map of the patients'
# paired clr changes, so a bootstrap replicate is a weighted column mean followed
# by a fixed projection. Mathematically identical to est_E3() for one visit
# (checked in 06_VALIDATION, V7a). Sourced only.
#
# INPUTS
# ---
#   e3_composition() output (03_helpers_estimands.R)
#
# ANALYSES
# ---
#   1. e3_pair_matrix()   — patients x bounded-bin matrix of paired clr changes per visit
#   2. e3_projector()     — residual-from-spline (df = 3) projection matrix
#   3. e3_band_decision() — bootstrap max-t band of residuals + bump flags (bin indices)
#   4. simulate_e3()      — true-null plasmode with optional injected bump
#
# OUTPUTS
# ---
#   Functions only.
#
# ======================================================================

# Paired clr changes (follow-up minus M0) for the 19 bounded bins of the E3 composition.
e3_pair_matrix <- function(comp, visit) {
  b <- comp[is.finite(comp$x_log_mg), c("patient_id", "visit_code", "descriptor", "lower_mg", "x_log_mg", "clr")]
  lev <- unique(b$descriptor[order(b$lower_mg)])
  w <- tidyr::pivot_wider(b[, c("patient_id", "visit_code", "descriptor", "clr")],
                          names_from = descriptor, values_from = clr)
  m0 <- w[w$visit_code == "M0", c("patient_id", lev)]
  mv <- w[w$visit_code == visit, c("patient_id", lev)]
  j  <- intersect(m0$patient_id, mv$patient_id)
  D  <- as.matrix(mv[match(j, mv$patient_id), lev]) - as.matrix(m0[match(j, m0$patient_id), lev])
  rownames(D) <- j
  attr(D, "x") <- b$x_log_mg[match(lev, b$descriptor)]
  attr(D, "labels") <- lev
  D
}

# Residual operator R = I - X (X'X)^-1 X' for the natural-spline (df = 3) null.
e3_projector <- function(x) {
  X <- cbind(1, splines::ns(x, df = 3))
  diag(length(x)) - X %*% solve(crossprod(X), t(X))
}

# Bootstrap max-t band of the mean residual profile; bumps as in find_bumps()
# (>= 2 contiguous flagged bins, flanked by unflagged bins, not at the grid edges).
bump_index <- function(flag, min_len = 2L) {
  flag[is.na(flag)] <- FALSE
  r <- rle(flag); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  keep <- r$values & r$lengths >= min_len & starts > 1 & ends < length(flag)
  unlist(Map(seq, starts[keep], ends[keep]))
}
e3_band_decision <- function(D, R, B, seed, level = 0.95) {
  n <- nrow(D)
  res_obs <- drop(R %*% colMeans(D))
  set.seed(seed)
  idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = B)
  W <- t(apply(idx, 1, tabulate, nbins = n)) / n          # B x n weights
  res_b <- (W %*% D) %*% t(R)                             # B x bins
  se <- apply(res_b, 2, sd)
  tmax <- apply(abs(sweep(res_b, 2, res_obs)) / rep(se, each = B), 1, max)
  crit <- unname(stats::quantile(tmax, level))
  lo <- res_obs - crit * se; hi <- res_obs + crit * se
  list(res = res_obs, lo = lo, hi = hi, pos = bump_index(lo > 0), neg = bump_index(hi < 0))
}

# One plasmode dataset: resample patients from the TRUE-NULL matrix (observed mean
# residual removed from every patient, so the population residual is exactly 0 while
# the smooth tilt and the patient covariance are kept), optionally add a bump of
# `size` (log) to bins `hit`, then apply the decision rule with B_SIM replicates.
simulate_e3 <- function(D, R, hit = integer(0), size = 0, n_sim, b_sim, seed) {
  r_obs <- drop(R %*% colMeans(D))
  D0 <- sweep(D, 2, r_obs)
  if (length(hit)) D0[, hit] <- D0[, hit] + size
  n <- nrow(D0)
  purrr::map_dfr(seq_len(n_sim), \(i) {
    set.seed(seed + i)
    Di <- D0[sample.int(n, n, replace = TRUE), , drop = FALSE]
    dec <- e3_band_decision(Di, R, b_sim, seed + 1e6 + i)
    tibble::tibble(sim = i, any_bump = length(dec$pos) + length(dec$neg) > 0,
                   detected = length(hit) > 0 && any(dec$pos %in% hit))
  })
}
