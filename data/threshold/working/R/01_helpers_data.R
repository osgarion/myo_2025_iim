# ======================================================================
# Script:  01_helpers_data.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-13
# ======================================================================
#
# PURPOSE
# ---
# Data layer shared by all analysis scripts: loads the master RDS, runs the
# locked data-integrity checks (PASS/WARN/FAIL) and builds the analytic
# long tables for CUMULATIVE thresholds and exclusive BINS. Sourced only.
#
# INPUTS
# ---
#   INPUT_RDS (from 00_config.R)
#
# ANALYSES
# ---
#   1. load_master()        — read + structural assertions
#   2. qc_master()          — data checks D1–D20 as a PASS/WARN/FAIL table
#   3. prepare_bins() / prepare_cumulative() — analytic long tables
#   4. make_pairs()         — M0 vs follow-up complete pairs
#   5. output helpers       — theme, figure/table/meta writers
#
# OUTPUTS
# ---
#   Functions only.
#
# ======================================================================

# **********************************************************************
# 1. Load ----
# **********************************************************************
load_master <- function(path = INPUT_RDS) {
  acc <- readRDS(path)
  need <- c("record_level", "threshold_long", "bin_long", "qc", "metadata")
  miss <- setdiff(need, names(acc))
  if (length(miss)) stop("Master RDS lacks components: ", paste(miss, collapse = ", "))
  acc
}

# Canonical record-level covariates (single source for QC, season, valid days).
record_info <- function(acc) {
  acc$record_level |>
    tibble::as_tibble() |>
    transmute(
      patient_id  = as.character(patient_id),
      visit_code  = as.character(visit_code),
      eligible    = as.logical(qc_analysis_eligible),
      start_time,
      month_num   = as.integer(format(start_time, "%m")),
      season      = factor(as.character(season), levels = c("Winter", "Spring", "Summer", "Autumn")),
      n_valid_days = as.numeric(.data[["part5_personsummary_MM_L40M100V400_T5A5_Nvaliddays"]]),
      enmo_mean   = as.numeric(.data[["part2_summary_ENMO_fullRecordingMean"]]),
      subtype_code = as.character(subtype),
      subtype     = unname(SUBTYPE_LABELS[as.character(subtype)]),
      sex         = as.character(sex)
    ) |>
    group_by(patient_id) |>
    mutate(months_from_m0 = if (any(visit_code == "M0"))
      as.numeric(difftime(start_time, start_time[visit_code == "M0"][1], units = "days")) / 30.4375
      else NA_real_) |>
    ungroup()
}

# **********************************************************************
# 2. Data QC (D1–D20) ----
# **********************************************************************
qc_row <- function(id, check, status, detail) tibble::tibble(id = id, check = check, status = status, detail = detail)

qc_master <- function(acc) {
  rec <- record_info(acc)
  thr <- acc$threshold_long |> mutate(patient_id = as.character(patient_id), visit_code = as.character(visit_code))
  bin <- acc$bin_long       |> mutate(patient_id = as.character(patient_id), visit_code = as.character(visit_code))
  main <- rec |> filter(eligible %in% TRUE, visit_code %in% MAIN_VISITS)
  out <- list()

  dup <- anyDuplicated(rec[c("patient_id", "visit_code")]) > 0
  out$D1 <- qc_row("D1", "Unique patient x visit (record_level)", if (dup) "FAIL" else "PASS",
                   sprintf("%d records", nrow(rec)))
  n_inel <- sum(!(rec$eligible %in% TRUE))
  out$D2 <- qc_row("D2", "QC eligibility", if (n_inel) "WARN" else "PASS",
                   sprintf("%d of %d not eligible", n_inel, nrow(rec)))
  out$D3 <- qc_row("D3", "Main sample (M0-M18)", "INFO",
                   sprintf("%d recordings / %d patients; %s", nrow(main), n_distinct(main$patient_id),
                           paste(names(table(main$visit_code)[MAIN_VISITS]), table(main$visit_code)[MAIN_VISITS], sep = "=", collapse = ", ")))
  has <- function(v) main |> group_by(patient_id) |> summarise(ok = all(c("M0", v) %in% visit_code)) |> pull(ok) |> sum()
  out$D4 <- qc_row("D4", "Complete pairs with M0", "INFO",
                   paste(FOLLOWUP_VISITS, vapply(FOLLOWUP_VISITS, has, numeric(1)), sep = "=", collapse = ", "))
  no_m0 <- main |> group_by(patient_id) |> summarise(m0 = any(visit_code == "M0")) |> filter(!m0)
  out$D5 <- qc_row("D5", "Patients without M0", if (nrow(no_m0)) "WARN" else "PASS",
                   sprintf("%d patients (enter LMM, not paired analyses)", nrow(no_m0)))

  grid_ok <- setequal(unique(thr$threshold_mg), THRESHOLD_GRID_MG) &&
    all(table(paste(thr$patient_id, thr$visit_code)) == length(THRESHOLD_GRID_MG))
  out$D6 <- qc_row("D6", "Threshold grid complete", if (grid_ok) "PASS" else "FAIL",
                   paste(sort(unique(thr$threshold_mg)), collapse = ","))
  nb <- n_distinct(bin$lower_mg)
  open_bin <- any(!is.finite(bin$upper_mg))
  out$D7 <- qc_row("D7", "Exclusive bins", if (nb == 23 && !open_bin) "PASS" else "WARN",
                   sprintf("%d bounded bins %g-%g mg; open >400 bin present: %s", nb,
                           min(bin$lower_mg), max(bin$upper_mg, na.rm = TRUE), open_bin))

  mono <- thr |> arrange(patient_id, visit_code, threshold_mg) |>
    group_by(patient_id, visit_code) |> summarise(v = sum(diff(minutes_above) > 1e-8), .groups = "drop")
  out$D8 <- qc_row("D8", "Cumulative monotonicity", if (sum(mono$v)) "FAIL" else "PASS",
                   sprintf("%d violations", sum(mono$v)))

  tw <- thr |> select(patient_id, visit_code, threshold_mg, minutes_above)
  chk <- bin |> select(patient_id, visit_code, lower_mg, upper_mg, minutes_bin) |>
    left_join(tw |> rename(lower_mg = threshold_mg, t_lo = minutes_above), by = c("patient_id", "visit_code", "lower_mg")) |>
    left_join(tw |> rename(upper_mg = threshold_mg, t_hi = minutes_above), by = c("patient_id", "visit_code", "upper_mg"))
  err9 <- max(abs(chk$minutes_bin - (chk$t_lo - chk$t_hi)), na.rm = TRUE)
  out$D9 <- qc_row("D9", "bin = T_lo - T_hi", if (err9 < 1e-8) "PASS" else "FAIL", sprintf("max |error| = %.2e", err9))
  rec_sum <- bin |> group_by(patient_id, visit_code) |> summarise(s = sum(minutes_bin), .groups = "drop") |>
    left_join(tw |> filter(threshold_mg == 20) |> select(patient_id, visit_code, t20 = minutes_above), by = c("patient_id", "visit_code")) |>
    left_join(tw |> filter(threshold_mg == 400) |> select(patient_id, visit_code, t400 = minutes_above), by = c("patient_id", "visit_code"))
  err10 <- max(abs(rec_sum$s - (rec_sum$t20 - rec_sum$t400)))
  out$D10 <- qc_row("D10", "sum(bins) = T20 - T400", if (err10 < 1e-8) "PASS" else "FAIL", sprintf("max |error| = %.2e", err10))
  err11 <- max(abs(bin$minutes_per_10mg - bin$minutes_bin * 10 / bin$bin_width_mg))
  out$D11 <- qc_row("D11", "minutes_per_10mg consistent", if (err11 < 1e-8) "PASS" else "FAIL", sprintf("max |error| = %.2e", err11))
  n_neg <- sum(bin$minutes_bin < -1e-8)
  out$D12 <- qc_row("D12", "Non-negative bins", if (n_neg) "FAIL" else "PASS", sprintf("%d negative", n_neg))

  bm <- bin |> filter(visit_code %in% MAIN_VISITS)
  z <- bm |> group_by(lower_mg) |> summarise(zero = mean(minutes_bin == 0), le1 = mean(minutes_bin <= 1))
  out$D13 <- qc_row("D13", "Exact zeros vs near-zero in bins", "INFO",
                    sprintf("max exact-zero share %.1f%%; bins with >50%% <=1 min/day: %s", 100 * max(z$zero),
                            paste(z$lower_mg[z$le1 > .5], collapse = ",")))
  vd <- table(main$n_valid_days)
  out$D14 <- qc_row("D14", "Valid days (main sample)", if (any(main$n_valid_days < 5, na.rm = TRUE)) "WARN" else "PASS",
                    paste(names(vd), vd, sep = "d=", collapse = ", "))
  out$D15 <- qc_row("D15", "Month/season from start_time", if (all(rec$month_num == acc$record_level$acc_month_num)) "PASS" else "WARN",
                    "season/month taken from local start_time (acc_date is shifted by 1 day in some records; not used)")
  unk <- setdiff(unique(rec$subtype_code), c(names(SUBTYPE_LABELS), "1", "UCTD/SSc?"))
  out$D16 <- qc_row("D16", "Subtype coding", if (length(unk)) "WARN" else "PASS",
                    sprintf("undefined codes (excluded from subtype analyses): %s", paste(unk, collapse = ",")))
  tm <- main |> filter(visit_code != "M0", is.finite(months_from_m0)) |>
    mutate(visit_code = factor(visit_code, levels = FOLLOWUP_VISITS)) |> group_by(visit_code) |>
    summarise(r = sprintf("%s %.1f (%.1f-%.1f)", first(visit_code), median(months_from_m0), min(months_from_m0), max(months_from_m0)))
  out$D18 <- qc_row("D18", "Months from M0: median (range)", "WARN", paste(tm$r, collapse = "; "))
  sv <- main |> count(visit_code, season) |> group_by(visit_code) |> mutate(p = round(100 * n / sum(n))) |> ungroup()
  out$D19 <- qc_row("D19", "Season x visit", if (min(sv$n) < 5) "WARN" else "INFO",
                    paste(sv$visit_code, sv$season, sv$n, sep = ":", collapse = "; "))
  bind_rows(out)
}

# **********************************************************************
# 3. Analytic tables ----
# **********************************************************************
# kind = "bins": value = minutes_bin (min/day within the exclusive bin)
# kind = "cumulative": value = minutes_above (min/day above threshold T)
prepare_spectrum <- function(acc, kind = c("bins", "cumulative"), visits = MAIN_VISITS) {
  kind <- match.arg(kind)
  rec  <- record_info(acc) |> select(patient_id, visit_code, eligible, season, month_num,
                                     months_from_m0, n_valid_days, subtype, subtype_code)
  if (kind == "bins") {
    d <- acc$bin_long |>
      transmute(patient_id = as.character(patient_id), visit_code = as.character(visit_code),
                lower_mg = as.numeric(lower_mg), upper_mg = as.numeric(upper_mg),
                width_mg = as.numeric(bin_width_mg), intensity = (lower_mg + upper_mg) / 2,
                descriptor = paste0(lower_mg, "-", upper_mg),
                value = as.numeric(minutes_bin), density = as.numeric(minutes_per_10mg))
  } else {
    d <- acc$threshold_long |>
      transmute(patient_id = as.character(patient_id), visit_code = as.character(visit_code),
                lower_mg = as.numeric(threshold_mg), upper_mg = NA_real_, width_mg = NA_real_,
                intensity = lower_mg, descriptor = paste0("T", lower_mg),
                value = as.numeric(minutes_above), density = value)
  }
  d <- d |>
    inner_join(rec, by = c("patient_id", "visit_code")) |>
    filter(eligible %in% TRUE, visit_code %in% visits) |>
    mutate(visit_code = factor(visit_code, levels = visits),
           descriptor = factor(descriptor, levels = unique(descriptor[order(lower_mg)])),
           kind = kind)
  stopifnot(!anyDuplicated(d[c("patient_id", "visit_code", "descriptor")]),
            all(is.finite(d$value)), all(d$value >= -1e-8))
  d
}

make_pairs <- function(d) {
  base <- d |> filter(visit_code == "M0") |> select(patient_id, descriptor, baseline = value, baseline_density = density)
  d |> filter(visit_code != "M0") |>
    inner_join(base, by = c("patient_id", "descriptor")) |>
    mutate(visit = factor(as.character(visit_code), levels = FOLLOWUP_VISITS), delta = value - baseline)
}

# **********************************************************************
# 4. Output helpers ----
# **********************************************************************
theme_myo <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
      axis.title    = ggplot2::element_text(size = 12, face = "bold"),
      axis.text     = ggplot2::element_text(size = 11),
      panel.background      = ggplot2::element_rect(fill = "white"),
      plot.background       = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major      = ggplot2::element_line(color = "grey80", linewidth = 0.3),
      panel.grid.minor      = ggplot2::element_line(color = "grey92", linewidth = 0.2),
      legend.background     = ggplot2::element_rect(fill = "white"),
      legend.box.background = ggplot2::element_rect(fill = "white")
    )
}
scale_visit_linetype <- function(...) ggplot2::scale_linetype_manual(values = VISIT_LINETYPES, drop = TRUE, ...)
x_axis_mg <- function(lim = c(20, 400)) ggplot2::scale_x_continuous(breaks = seq(25, 400, 25), limits = lim,
                                                                     guide = ggplot2::guide_axis(angle = 45))

save_figure <- function(p, name, width = 11, height = 7, dir = PATHS$figures) {
  stem <- file.path(dir, paste0(date_prefix(), "_", name))
  ggplot2::ggsave(paste0(stem, ".png"),  p, width = width, height = height, dpi = 300, bg = "white")
  ggplot2::ggsave(paste0(stem, ".tiff"), p, width = width, height = height, dpi = 600, compression = "lzw", bg = "white")
  invisible(stem)
}
save_tables <- function(lst, name, dir = PATHS$tables) {
  f <- file.path(dir, paste0(date_prefix(), "_", name, ".xlsx"))
  rio::export(lst, f)
  invisible(f)
}
write_run_meta <- function(script, extra = list(), dir = PATHS$logs) {
  f <- file.path(dir, paste0(date_prefix(), "_", sub("\\.R$", "", script), "_sessionInfo.txt"))
  writeLines(c(paste("Script:", script), paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
               paste("Input:", INPUT_RDS), paste("Input md5:", INPUT_MD5),
               paste("Master created by:", acc_meta_string()),
               vapply(names(extra), function(n) paste0(n, ": ", paste(extra[[n]], collapse = ",")), ""),
               "", capture.output(sessionInfo())), f)
  invisible(f)
}
acc_meta_string <- function() {
  if (!exists("acc")) return("unknown")
  m <- acc$metadata
  paste0(m$script_version, " (", format(m$created, "%Y-%m-%d %H:%M"), ")")
}
