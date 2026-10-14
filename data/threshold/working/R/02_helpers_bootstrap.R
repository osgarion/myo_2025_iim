# ======================================================================
# Script:  02_helpers_bootstrap.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-13
# ======================================================================
#
# PURPOSE
# ---
# Patient-level cluster bootstrap engine. All patient draws are generated
# up-front in the main process from one seed, so results are identical for
# sequential and parallel runs and for any number of workers. Sourced only.
#
# INPUTS
# ---
#   Analytic long table with `patient_id` (from 01_helpers_data.R)
#
# ANALYSES
# ---
#   1. make_draws()          — B x N_patients matrix of resampled patient indices
#   2. resample_patients()   — rebuild a replicate; duplicated patients get new IDs
#   3. run_cluster_bootstrap() — apply a statistic to every replicate (furrr)
#   4. summarise_bootstrap() — pointwise percentile CI, bootstrap SE,
#                              simultaneous max-t band, failure counts
#
# OUTPUTS
# ---
#   Functions only.
#
# ======================================================================

# **********************************************************************
# 1. Draws ----
# **********************************************************************
make_draws <- function(patients, B = B_BOOT, seed = SEED) {
  withr_seed <- function(code) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    code
  }
  np <- length(patients)
  withr_seed(matrix(sample.int(np, np * B, replace = TRUE), nrow = B, ncol = np, byrow = TRUE))
}

# **********************************************************************
# 2. One replicate ----
# **********************************************************************
# Each selected copy of a patient becomes a distinct cluster (patient_id = "<id>__b<k>"),
# so LMM random intercepts are not shared between copies.
resample_patients <- function(d, patients, idx) {
  split_rows <- split(seq_len(nrow(d)), factor(d$patient_id, levels = patients))
  rows <- split_rows[idx]
  out  <- d[unlist(rows, use.names = FALSE), , drop = FALSE]
  out$patient_id_orig <- out$patient_id
  out$patient_id <- rep(paste0(patients[idx], "__b", seq_along(idx)), lengths(rows))
  out
}

# **********************************************************************
# 3. Run ----
# **********************************************************************
# stat_fun(d) must return a tibble with key columns + numeric column `estimate`
# (additional columns are kept). Errors in a replicate are captured.
# Checkpointing: replicates are processed in chunks; each finished chunk is saved to
# `checkpoint_dir` and reused on a re-run with identical B, seed, patients and `tag`,
# so an interrupted run resumes where it stopped.
run_cluster_bootstrap <- function(d, stat_fun, B = B_BOOT, seed = SEED, workers = N_WORKERS,
                                  packages = c("dplyr", "tidyr", "purrr", "tibble", "glmmTMB", "splines"),
                                  checkpoint_dir = NULL, tag = "boot", chunk_size = 100L) {
  patients <- sort(unique(d$patient_id))
  draws <- make_draws(patients, B, seed)
  one <- function(b) {
    rep_d <- resample_patients(d, patients, draws[b, ])
    res <- tryCatch(stat_fun(rep_d), error = function(e) NULL)
    if (is.null(res)) return(tibble::tibble(replicate = b, failed = TRUE))
    dplyr::mutate(res, replicate = b, failed = FALSE)
  }
  t0 <- Sys.time()
  if (workers > 1L) {
    old_plan <- future::plan(future::multisession, workers = workers)
    on.exit(future::plan(old_plan), add = TRUE)
  }
  run_ids <- function(ids) {
    if (workers > 1L)
      furrr::future_map(ids, one, .options = furrr::furrr_options(seed = NULL, packages = packages,
                                                                  globals = TRUE, scheduling = 2))
    else purrr::map(ids, one)
  }
  chunks <- split(seq_len(B), ceiling(seq_len(B) / chunk_size))
  key <- substr(digest_key(c(tag, B, seed, patients)), 1, 10)
  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  out <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    f <- if (!is.null(checkpoint_dir)) file.path(checkpoint_dir, sprintf("%s_%s_chunk%03d.rds", tag, key, i))
    if (!is.null(f) && file.exists(f)) { out[[i]] <- readRDS(f); next }
    out[[i]] <- dplyr::bind_rows(run_ids(chunks[[i]]))
    if (!is.null(f)) {
      dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)  # survives folder removal mid-run
      saveRDS(out[[i]], f)
    }
    message(sprintf("[bootstrap %s] chunk %d/%d done (%.1f min)", tag, i, length(chunks),
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  out <- dplyr::bind_rows(out)
  attr(out, "settings") <- tibble::tibble(B = B, seed = seed, workers = workers,
                                          n_patients = length(patients),
                                          n_failed_replicates = sum(dplyr::distinct(out, replicate, failed)$failed),
                                          runtime_min = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2))
  out
}

# Stable short key for checkpoint file names (base R only).
digest_key <- function(x) {
  f <- tempfile(); on.exit(unlink(f))
  writeLines(as.character(x), f)
  unname(tools::md5sum(f))
}

# **********************************************************************
# 4. Summaries ----
# **********************************************************************
# observed: tibble with keys + `estimate`; boot: output of run_cluster_bootstrap().
# The max-t band is simultaneous over all rows sharing the same `band_by` group
# (e.g. all 23 bins within one visit and one model).
summarise_bootstrap <- function(observed, boot, keys, band_by = NULL, level = 0.95) {
  a <- (1 - level) / 2
  b_ok <- boot |> filter(!failed, is.finite(estimate))
  pt <- b_ok |>
    group_by(across(all_of(keys))) |>
    summarise(n_boot = n(), boot_se = sd(estimate), ci_low = unname(quantile(estimate, a)),
              ci_high = unname(quantile(estimate, 1 - a)), p_positive = mean(estimate > 0),
              .groups = "drop")
  res <- observed |> left_join(pt, by = keys)
  if (!is.null(band_by)) {
    tstat <- b_ok |>
      inner_join(res |> select(all_of(keys), est_obs = estimate, boot_se), by = keys) |>
      filter(boot_se > 0) |>
      group_by(across(all_of(c(band_by, "replicate")))) |>
      summarise(t_max = max(abs(estimate - est_obs) / boot_se), .groups = "drop") |>
      group_by(across(all_of(band_by))) |>
      summarise(crit_simultaneous = unname(quantile(t_max, level)), .groups = "drop")
    res <- res |> left_join(tstat, by = band_by) |>
      mutate(band_low = estimate - crit_simultaneous * boot_se,
             band_high = estimate + crit_simultaneous * boot_se)
  }
  failed <- boot |> distinct(replicate, failed) |> summarise(n_failed = sum(failed)) |> pull(n_failed)
  attr(res, "n_failed_replicates") <- failed
  res
}
