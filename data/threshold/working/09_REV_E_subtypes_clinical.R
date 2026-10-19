# ======================================================================
# Script:  09_REV_E_subtypes_clinical.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-16
# ======================================================================
#
# PURPOSE
# ---
# SCOPE: patient-level summaries of change (dT100, clr tilt, dIG); no bin selection.
# A07 revision, block E (docs/261016_A07_critical_review.md). EXPLORATORY.
# Checks availability, distributions and collinearity of baseline clinical severity
# indicators (physician global activity VAS = ph_VAS, MMT8, FI2 and the MDAAT muscle /
# extramuscular components) and fits parsimonious random-intercept models of change by
# subtype (DM / IMNM / ASyS), unadjusted and adjusted for baseline ph_VAS (primary, per
# user decision 2026-10-16) or MMT8 (sensitivity), with severity-by-visit interaction and
# baseline-outcome adjustment as further variants. The descriptive F7 is not changed.
#
# INPUTS
# ---
#   MYOII_ACC_OUTCOMES_MASTER_V1.rds — via R/00_config.R
#   klinicka_data_12102026.xlsx      — baseline clinical covariates (locate_clinical())
#
# ANALYSES
# ---
#   1. Linkage and availability by subtype
#   2. Distributions and collinearity (Spearman), construct check of ph_VAS
#   3. Models E_a-E_e per outcome x severity indicator; LRTs (ML)
#   4. Patient-level cluster bootstrap of subtype contrasts (B = B_BOOT)
#   5. Figure F7b (adjusted vs unadjusted contrasts)
#
# OUTPUTS
# ---
#   RData:   output/RData/rev_E_objects.RData   (aggregated only; no patient rows)
#   Tables:  output/tables/YYMMDD_09_REV_E_subtypes_clinical.xlsx
#   Figures: output/figures/YYMMDD_09_E1_F7b_subtype_contrasts.png/.tiff
#   Logs:    output/logs/YYMMDD_09_REV_E_subtypes_clinical_sessionInfo.txt
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
check_packages(c("readxl", "lme4"))
for (f in c("01_helpers_data.R", "02_helpers_bootstrap.R", "03_helpers_estimands.R", "06_helpers_revision.R"))
  source(file.path(WORKING_DIR, "R", f))

acc <- load_master()
dl  <- build_analysis_long(acc)
CLIN_PATH <- locate_clinical()
CLIN_MD5  <- unname(tools::md5sum(CLIN_PATH))
message("[09] clinical file: ", CLIN_PATH, " (md5 ", CLIN_MD5, ")")
clin <- load_clinical(CLIN_PATH)
SEV_VARS <- c("ph_vas", "muscle_vas", "extramuscular_vas", "mmt8", "fi2", "haq", "disease_duration", "gc_mg_day")

# **********************************************************************
# 2. Linkage and availability ----
# **********************************************************************
rec <- record_info(acc) |> filter(eligible %in% TRUE, visit_code %in% MAIN_VISITS)
pts <- rec |> distinct(patient_id, subtype_code, subtype)
stopifnot(!anyDuplicated(clin[c("patient_id", "visit_code")]))
m0c <- clin |> filter(visit_code == "M0") |> select(patient_id, subtype_clin, all_of(SEV_VARS))
pinfo <- pts |> left_join(m0c, by = "patient_id") |>
  mutate(has_m0_acc = patient_id %in% rec$patient_id[rec$visit_code == "M0"],
         subtype = factor(subtype, levels = c("DM", "IMNM", "ASyS")))
linkage <- tibble::tibble(
  master_patients = nrow(pts), matched_in_clinical = sum(pts$patient_id %in% clin$patient_id),
  subtype_code_agreement = sum(pinfo$subtype_code == pinfo$subtype_clin, na.rm = TRUE),
  subtype_code_compared = sum(!is.na(pinfo$subtype_clin)))

chg <- patient_change_summaries(acc, dl)
with_pair <- unique(chg$patient_id)
availability <- pinfo |> mutate(subtype = ifelse(is.na(subtype), "PM / UCTD (excluded)", as.character(subtype))) |>
  group_by(subtype) |>
  summarise(patients = n(), with_M0_recording = sum(has_m0_acc), with_follow_up_pair = sum(patient_id %in% with_pair),
            across(all_of(c("ph_vas", "mmt8", "fi2", "muscle_vas", "extramuscular_vas", "disease_duration", "gc_mg_day")),
                   \(v) sum(!is.na(v) & patient_id %in% with_pair), .names = "pairs_with_{.col}"), .groups = "drop")

# **********************************************************************
# 3. Distributions, collinearity, construct ----
# **********************************************************************
base_set <- pinfo |> filter(patient_id %in% with_pair)
q3 <- function(v) if (sum(!is.na(v)) < 3) NA_character_ else
  sprintf("%.1f [%.1f; %.1f] (n=%d)", median(v, na.rm = TRUE), quantile(v, .25, na.rm = TRUE), quantile(v, .75, na.rm = TRUE), sum(!is.na(v)))
distributions <- base_set |> filter(!is.na(subtype)) |>
  bind_rows(base_set |> mutate(subtype = "All (incl. PM/UCTD)")) |>
  mutate(subtype = as.character(subtype)) |>
  group_by(subtype) |> summarise(across(all_of(SEV_VARS), q3), .groups = "drop")
cor_spearman <- stats::cor(base_set[SEV_VARS], use = "pairwise.complete.obs", method = "spearman")
cor_n <- crossprod(!is.na(as.matrix(base_set[SEV_VARS])))
collinearity <- tibble::as_tibble(cor_spearman, rownames = "variable")
construct <- base_set |> filter(!is.na(subtype)) |> group_by(subtype) |>
  summarise(n = sum(!is.na(ph_vas) & !is.na(muscle_vas)),
            rho_phvas_muscle = suppressWarnings(cor(ph_vas, muscle_vas, use = "complete.obs", method = "spearman")),
            rho_phvas_extramuscular = suppressWarnings(cor(ph_vas, extramuscular_vas, use = "complete.obs", method = "spearman")),
            rho_phvas_mmt8 = suppressWarnings(cor(ph_vas, mmt8, use = "complete.obs", method = "spearman")),
            median_muscle_vas = median(muscle_vas, na.rm = TRUE), median_extramuscular_vas = median(extramuscular_vas, na.rm = TRUE),
            .groups = "drop")
# Variance inflation of the severity term given subtype (R2 of severity ~ subtype).
vif_sev <- purrr::map_dfr(c(ph_vas = "ph_vas", mmt8 = "mmt8"), \(v) {
  d <- base_set |> filter(!is.na(subtype), !is.na(.data[[v]]))
  r2 <- summary(stats::lm(stats::as.formula(paste(v, "~ subtype")), d))$r.squared
  tibble::tibble(severity = v, n_patients = nrow(d), r2_on_subtype = r2, vif = 1 / (1 - r2))
})

# **********************************************************************
# 4. Models ----
# **********************************************************************
OUTCOMES <- c(dT100 = "t100", tilt = "tilt", dIG = "ig")
SEVERITY <- c(ph_vas = 10, mmt8 = 5)          # scaling: per 10 mm VAS, per 5 MMT8 points
model_data <- function(d_chg, info, outcome, sev) {
  y <- if (outcome == "tilt") d_chg$tilt else d_chg[[paste0("chg_", outcome)]]
  b <- if (outcome == "tilt") rep(0, nrow(d_chg)) else d_chg[[paste0("base_", outcome)]]
  out <- d_chg |> transmute(patient_id, visit, chg = y, base = b) |>
    inner_join(info |> select(patient_id, subtype, sevraw = all_of(if (is.null(sev)) "ph_vas" else sev)), by = "patient_id") |>
    filter(!is.na(subtype), !is.na(chg), !is.na(base))
  if (!is.null(sev)) out <- out |> filter(!is.na(sevraw))
  pm <- out |> distinct(patient_id, sevraw, .keep_all = TRUE)
  out |> mutate(sev = if (is.null(sev)) 0 else (sevraw - mean(pm$sevraw)) / SEVERITY[[sev]],
                base = base - mean(pm$base), subtype = droplevels(subtype), visit = factor(visit, levels = FOLLOWUP_VISITS))
}
model_set <- function(d_chg, info) {
  purrr::imap_dfr(OUTCOMES, \(o, onm) bind_rows(
    fit_subtype_models(model_data(d_chg, info, o, NULL), SUBTYPE_MODELS["E_a_unadjusted"]) |>
      mutate(severity = "none", sample = "all subtype patients"),
    purrr::map_dfr(names(SEVERITY), \(s) {
      d <- model_data(d_chg, info, o, s)
      # No baseline value exists for the tilt (a change by construction): E_e is not fitted.
      adj <- SUBTYPE_MODELS[setdiff(names(SUBTYPE_MODELS), c("E_a_unadjusted", if (o == "tilt") "E_e_sev_baseline"))]
      bind_rows(fit_subtype_models(d, SUBTYPE_MODELS["E_a_unadjusted"]) |> mutate(sample = paste("complete", s)),
                fit_subtype_models(d, adj) |> mutate(sample = paste("complete", s))) |>
        mutate(severity = s)
    })
  ) |> mutate(outcome = onm))
}
info_mod <- pinfo |> select(patient_id, subtype, ph_vas, mmt8)
obs_models <- model_set(chg, info_mod)
lrt <- purrr::imap_dfr(OUTCOMES, \(o, onm) purrr::map_dfr(names(SEVERITY), \(s)
  subtype_lrt(model_data(chg, info_mod, o, s)) |> mutate(outcome = onm, severity = s)))
cell_n <- purrr::imap_dfr(OUTCOMES, \(o, onm) bind_rows(
  model_data(chg, info_mod, o, NULL) |> count(subtype, visit) |> mutate(sample = "all subtype patients"),
  purrr::map_dfr(names(SEVERITY), \(s) model_data(chg, info_mod, o, s) |> count(subtype, visit) |> mutate(sample = paste("complete", s)))
) |> mutate(outcome = onm))

# **********************************************************************
# 5. Cluster bootstrap of contrasts ----
# **********************************************************************
# Each replicate resamples patients (with all their visits and their baseline covariates).
boot_in <- chg |> inner_join(info_mod, by = "patient_id") |> filter(!is.na(subtype))
stat_E <- function(d) {
  info <- d |> distinct(patient_id, subtype, ph_vas, mmt8)
  model_set(d |> select(-subtype, -ph_vas, -mmt8), info) |>
    select(outcome, severity, sample, model, visit, contrast, estimate)
}
boot_E <- run_cluster_bootstrap(boot_in, stat_E, B = B_BOOT, seed = SEED, workers = N_WORKERS,
                                packages = c("dplyr", "tidyr", "purrr", "tibble", "lme4", "splines"),
                                checkpoint_dir = file.path(PATHS$rdata, "boot_chunks"), tag = "rev_E", chunk_size = 250L)
keysE <- c("outcome", "severity", "sample", "model", "visit", "contrast")
models_summ <- summarise_bootstrap(obs_models, boot_E, keys = keysE) |>
  mutate(outcome = factor(outcome, levels = names(OUTCOMES)))

# **********************************************************************
# 6. Figure F7b ----
# **********************************************************************
show_models <- tibble::tribble(
  ~severity, ~sample,                  ~model,             ~label,
  "none",    "all subtype patients",   "E_a_unadjusted",   "1 unadjusted, all patients",
  "ph_vas",  "complete ph_vas",        "E_a_unadjusted",   "2 unadjusted, ph_VAS sample",
  "ph_vas",  "complete ph_vas",        "E_b_sev",          "3 + ph_VAS (primary adjusted)",
  "ph_vas",  "complete ph_vas",        "E_c_sev_by_visit", "4 + ph_VAS x visit",
  "ph_vas",  "complete ph_vas",        "E_e_sev_baseline", "5 + ph_VAS + baseline outcome",
  "mmt8",    "complete mmt8",          "E_b_sev",          "6 + MMT8 (sensitivity)")
d_f7b <- models_summ |> filter(grepl(" - ", contrast), contrast != "ASyS - IMNM") |>
  inner_join(show_models, by = c("severity", "sample", "model")) |>
  mutate(visit = factor(visit, levels = c(FOLLOWUP_VISITS, "mean M3-M18")),
         label = factor(label, levels = rev(show_models$label)),
         outcome = factor(outcome, levels = names(OUTCOMES),
                          labels = c("dT100 (min/day)", "clr tilt (per log mg)", "dIG")))
fig_E1 <- ggplot(d_f7b, aes(estimate, label, shape = contrast)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y", position = position_dodge(width = .6)) +
  geom_point(size = 1.8, position = position_dodge(width = .6)) +
  facet_grid(visit ~ outcome, scales = "free_x") +
  scale_shape_manual(values = c("IMNM - DM" = 16, "ASyS - DM" = 1)) +
  labs(x = "Difference in mean change from M0 (95% cluster-bootstrap CI)", y = NULL, shape = "Contrast",
       title = "Subtype differences in change, unadjusted and adjusted (exploratory)",
       subtitle = "Random-intercept models of change; DM = reference; PM and UCTD/SSc excluded",
       caption = sprintf("ph_VAS per 10 mm and MMT8 per 5 points, centred; B = %d. Rows 1 vs 2 separate sample restriction from adjustment.", B_BOOT)) +
  theme_myo() + theme(legend.position = "bottom", strip.text.y = element_text(size = 9), axis.text.y = element_text(size = 9))

# **********************************************************************
# 7. Save ----
# **********************************************************************
save_figure(fig_E1, "09_E1_F7b_subtype_contrasts", width = 14, height = 11)
save_tables(list(linkage = linkage, availability = availability, distributions = distributions,
                 collinearity_spearman = collinearity, collinearity_pairwise_n = tibble::as_tibble(cor_n, rownames = "variable"),
                 construct_ph_vas = construct, vif_severity = vif_sev, cell_n = cell_n,
                 models = models_summ, lrt = lrt, bootstrap_settings = attr(boot_E, "settings")),
            "09_REV_E_subtypes_clinical")
rev_E_objects <- list(linkage = linkage, availability = availability, distributions = distributions,
                      collinearity = collinearity, cor_n = cor_n, construct = construct, vif_sev = vif_sev,
                      cell_n = cell_n, models = models_summ, lrt = lrt, show_models = show_models,
                      boot_settings = attr(boot_E, "settings"),
                      clinical = list(file = basename(CLIN_PATH), md5 = CLIN_MD5), figs = list(E1 = fig_E1))
save_public("rev_E_objects", file = file.path(PATHS$rdata, "rev_E_objects.RData"))
write_run_meta("09_REV_E_subtypes_clinical.R", list(B = B_BOOT, seed = SEED, clinical = CLIN_PATH, clinical_md5 = CLIN_MD5))
message("09 REV E finished.")
