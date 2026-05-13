# **********************************************************************
# Script:  Script_myo_2025_iim_02.R
# Project: myo_2025_iim
# Author:  Jiří Baloun
# Date:    2026-05-07
# **********************************************************************
#
# PURPOSE
# ---
# Reanalysis of myokine-clinical associations using the revised dataset
# (patients with excessive disease/treatment duration excluded, *_revize
# sheets). Runs lmer covariate models M0 vs M6 and post-hoc subgroup
# visualisations stratified by sex, sub-diagnosis, anti-Jo1, anti-HMGCR.
#
# INPUTS
# ---
#   IMMET_terapie a klinická data_01102025.xlsx (*_revize sheet)
#     — revised clinical data
#   Myokiny sérum_pro stat.xlsx (*_revize sheet)
#     — revised myokine concentrations (mstn, fst, fstl3, akta)
#
# ANALYSES
# ---
#   1. Data loading from *_revize sheets and pipeline to d04_sel2_rev
#   2. Covariate models 2 — lmer for all dep x indep combinations
#      - model without confounder (_without)
#      - model with age (_age), GK daily dose (_gk), BCM (_bcm),
#        creatinine (_creatine)
#   3. Post-hoc subgroup visualisation (for significant results p < 0.06)
#      - stratified by: sex, sub-diagnosis, anti-Jo1, anti-HMGCR
#      - separately for each covariate model
#   4. PLS-DA — M0 data, grouping by mda_categories
#      - predictors: all clinical + myokine variables (log-transforms applied)
#      - rename labels from sheet '260507_MDA categories' if available
#      - output: group plot + loadings plot (TIFF)
#
# OUTPUTS
# ---
#   Tables:
#     output/tables/YYMMDD_covariates_revize_01.xlsx
#       — main covariate table (all dep x indep x confounder combinations)
#   Figures (PDF):
#     output/figures/YYMMDD_revize_age_{sex,sub,jo1,hmgcr}.pdf
#     output/figures/YYMMDD_revize_gk_{sex,sub,jo1,hmgcr}.pdf
#     output/figures/YYMMDD_revize_bcm_{sex,sub,jo1,hmgcr}.pdf
#     output/figures/YYMMDD_revize_creatine_{sex,sub,jo1,hmgcr}.pdf
#   Figures (TIFF):
#     output/figures/cluster analyses/YYMMDD_revize_m0_pls_mda_categories_group.tiff
#     output/figures/cluster analyses/YYMMDD_revize_m0_pls_mda_categories_loadings.tiff
#   RData:
#     reports/YYMMDD_covariates_revize_01.RData
#       — pre-computed objects for the report
#
# **********************************************************************


# **********************************************************************
# 1. Libraries & functions ----
# **********************************************************************

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, purrr, conflicted, renv)

options(myo_pacman_update = FALSE, myo_pacman_install = FALSE)

# Source only function files — NOT OBJ_01.R (data loaded fresh below)
list.files("scripts/functions/",
           pattern = "^FUN_|ggplot_the_model|theme_Publication",
           full.names = TRUE) |>
  purrr::walk(source)

# multilevelmod must be loaded explicitly to register the "lmer" engine
# for linear_reg() in tidymodels (FUN_01.R loads it via pacman but
# engine registration can silently fail in some environments)
library(multilevelmod)

back_up("scripts/Script_myo_2025_iim_02.R")


# **********************************************************************
# 2. Variables ----
# **********************************************************************

var_indep_01 <- c("mstn", "fst", "fstl3", "akta")

var_dep_01 <- c("mmt8_total", "mmt10_total", "fi_2", "borg10", "haq",
                "sf36_mcs", "sf36_pf", "sf36_rp", "sf36_bp", "sf36_gh", "sf36_vt",
                "sf36_sf", "sf36_re", "sf36_mh", "sf36_pcs", "crp", "ast", "alt",
                "ck", "ld", "kreatinin_umol_l", "mitax", "myoglobin", "myoact",
                "physician_activity_vas", "muscle_disease_activity",
                "odpoved_na_terapii_m0_vs_m6")

folder_path <- "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Data/"

ncores <- getOption("myo_ncores", default = max(1L, future::availableCores() - 1L))


# **********************************************************************
# 3. Data loading — *_revize sheets ----
# **********************************************************************

## helper: find the *_revize sheet in an Excel file ----
get_revize_sheet <- function(path) {
  sheets <- readxl::excel_sheets(path)
  rev_sheet <- sheets[grepl("_revize$", sheets, ignore.case = TRUE)]
  if (length(rev_sheet) == 0) stop("No *_revize sheet found in: ", path)
  if (length(rev_sheet) > 1) message("Multiple *_revize sheets found — using first: ", rev_sheet[1])
  rev_sheet[1]
}

## myokines ----
path_myo_01 <- file_to_load(phrase = "Myokiny", folder_path = folder_path) |> tail(1)

d01_myokiny_rev <- import(path_myo_01, which = get_revize_sheet(path_myo_01)) |>
  select(!starts_with("d")) |>
  mutate(across(where(is.character), ~na_if(.x, "x")),
         across(-c(1:4), ~as.numeric(.x))) |>
  pivot_longer(
    cols      = matches("^[A-Za-z0-9]+ M\\d+$"),
    names_to  = c("protein", "poradie_vysetrenia"),
    names_sep = " ",
    values_to = "concentration"
  ) |>
  pivot_wider(names_from = protein, values_from = concentration) |>
  janitor::clean_names() |>
  mutate(poradie_vysetrenia = factor(poradie_vysetrenia, levels = c("M0", "M3", "M6", "M18")),
         bbm_id  = factor(bbm_id),
         pohlavi = factor(pohlavi))

pre_filt_rev <- d01_myokiny_rev |>
  select(maju_klinicke_data_pre_odfiltrovanie) |>
  na.omit() |>
  distinct() |>
  pull()

## clinical data ----
path_clin_01 <- file_to_load(phrase = "IMMET_terapie", folder_path = folder_path) |> tail(1)

d02_clinics_rev <- import(path_clin_01, which = get_revize_sheet(path_clin_01)) |>
  mutate(across(where(is.character), ~na_if(.x, "x"))) |>
  janitor::clean_names() |>
  tidyr::separate(col = bbm_id, sep = "_",
                  into = c("bbm_id", "bbm_collection"), extra = "merge") |>
  mutate(
    bbm_id                      = factor(bbm_id),
    datum_vysetreni             = if_else(str_detect(datum_vysetreni, "^x"),
                                          as.Date(NA),
                                          openxlsx::convertToDate(datum_vysetreni)),
    datum_biopsie               = openxlsx::convertToDate(datum_biopsie),
    date_of_bioimpedance        = openxlsx::convertToDate(date_of_bioimpedance),
    datum_nar                   = as.Date(datum_nar),
    bbm_collection              = factor(bbm_collection),
    poradie_vysetrenia          = factor(poradie_vysetrenia, levels = c("M0", "M3", "M6", "M18")),
    pohlavi                     = factor(pohlavi),
    odpoved_na_terapii_m0_vs_m6 = factor(odpoved_na_terapii_m0_vs_m6,
                                         levels = c("No", "Minimal", "Moderate", "Major")),
    odpoved_na_terapii_m0_vs_m18 = factor(odpoved_na_terapii_m0_vs_m18),
    podtyp_nemoci_zjednoduseny  = factor(podtyp_nemoci_zjednoduseny),
    mda_categories              = factor(mda_categories),
    mmt8_categories             = factor(mmt8_categories),
    kurak_nejmene_100_cigaret   = factor(kurak_nejmene_100_cigaret)
  ) |>
  mutate(across(where(is_01_col), as.factor)) |>
  mutate(across((vek:last_col()) & where(is.character), as.numeric)) |>
  group_by(bbm_id) |>
  fill(kurak_nejmene_100_cigaret) |>
  ungroup()

## merge into d03_rev ----
d03_rev <- d02_clinics_rev |>
  full_join(d01_myokiny_rev |>
              select(-maju_klinicke_data_pre_odfiltrovanie),
            by = c("projekt_id", "poradie_vysetrenia", "bbm_id", "pohlavi")) |>
  filter(projekt_id %in% pre_filt_rev) |>
  group_by(projekt_id) |>
  fill(podtyp_nemoci_zjednoduseny) |>
  fill(odpoved_na_terapii_m0_vs_m6, .direction = "downup") |>
  ungroup()

## d04_sel2_rev ----
d04_sel2_rev <- d03_rev |>
  mutate(odpoved_na_terapii_m0_vs_m6 =
           str_replace(odpoved_na_terapii_m0_vs_m6, "No Improvement", "No improvement")) |>
  group_by(projekt_id) |>
  mutate(across(c(jo_1, anti_hmgcr),
                ~ if_else(any(. == 1, na.rm = TRUE), 1L, 0L))) |>
  ungroup() |>
  select(projekt_id, pohlavi, vek, gk, davka_gk_mg_den, denna_davka_gk_mg_kg_bw,
         podtyp_nemoci_zjednoduseny, poradie_vysetrenia, bcm, jo_1, anti_hmgcr,
         all_of(var_dep_01), all_of(var_indep_01)) |>
  filter(poradie_vysetrenia %in% c("M0", "M3", "M6", "M18"))


# **********************************************************************
# 4. Covariate models 2 ----
# **********************************************************************

## functions ----

model_min_2 <- function(data, filter_na = TRUE) {

  if (filter_na) {
    data <- data |>
      select(var_dep_value, var_indep_value, poradie_vysetrenia, projekt_id) |>
      na.omit()
  }

  # Use lmerTest::lmer() directly — avoids multilevelmod engine registration
  # in furrr worker processes (workers are fresh R sessions without side effects)
  mod1 <- lmerTest::lmer(
    var_dep_value ~ poradie_vysetrenia + var_indep_value + (1 | projekt_id),
    data = data, REML = FALSE
  )

  params         <- parameters::model_parameters(mod1)
  omega          <- tryCatch(
    effectsize::omega_squared(mod1, partial = TRUE)$Omega2_partial[[2]],
    error = function(e) NA_real_
  )
  r2_mod         <- performance::r2_nakagawa(mod1)$R2_conditional
  beta_covar     <- params$Coefficient[[3]]
  lower_covar    <- params$CI_low[[3]]
  upper_covar    <- params$CI_high[[3]]
  p_val_covar    <- params$p[[3]]
  std_beta_covar <- tryCatch(
    parameters::standardize_parameters(mod1, method = "refit")$Std_Coefficient[[3]],
    error = function(e)
      parameters::standardize_parameters(mod1, method = "posthoc")$Std_Coefficient[[3]]
  )

  list(table = tibble(omega_sq_p     = omega,
                      r2_mod         = r2_mod,
                      beta_covar     = beta_covar,
                      std_beta_covar = std_beta_covar,
                      lower_covar    = lower_covar,
                      upper_covar    = upper_covar,
                      p_val_covar    = p_val_covar),
       model = mod1)
}

model_covar_2 <- function(data, covar, filter_na = TRUE) {

  if (filter_na) {
    data <- data |>
      select(var_dep_value, var_indep_value, poradie_vysetrenia, all_of(covar), projekt_id) |>
      na.omit()
  }

  rhs_terms <- c("poradie_vysetrenia", "var_indep_value", covar, "(1 | projekt_id)")
  fml <- as.formula(paste("var_dep_value ~", paste(rhs_terms, collapse = " + ")))

  mod1 <- lmerTest::lmer(fml, data = data, REML = FALSE)

  params         <- parameters::model_parameters(mod1)
  omega          <- tryCatch(
    effectsize::omega_squared(mod1, partial = TRUE)$Omega2_partial[[2]],
    error = function(e) NA_real_
  )
  r2_mod         <- performance::r2_nakagawa(mod1)$R2_conditional
  beta_covar     <- params$Coefficient[[3]]
  lower_covar    <- params$CI_low[[3]]
  upper_covar    <- params$CI_high[[3]]
  p_val_covar    <- params$p[[3]]
  std_beta_covar <- tryCatch(
    parameters::standardize_parameters(mod1, method = "refit")$Std_Coefficient[[3]],
    error = function(e)
      parameters::standardize_parameters(mod1, method = "posthoc")$Std_Coefficient[[3]]
  )

  list(table = tibble(omega_sq_p     = omega,
                      r2_mod         = r2_mod,
                      beta_covar     = beta_covar,
                      std_beta_covar = std_beta_covar,
                      lower_covar    = lower_covar,
                      upper_covar    = upper_covar,
                      p_val_covar    = p_val_covar),
       model = mod1)
}

## run models ----
options(renv.config.auto.snapshot = FALSE)
tic(); Sys.time()

future::plan(multisession, workers = ncores)

res_mixMod_covar_02_rev <- d04_sel2_rev |>
  filter(poradie_vysetrenia %in% c("M0", "M6")) |>
  select(any_of(var_dep_01 |> str_subset("odpoved_na_terapii_m0_vs_m6|kreatinin_umol_l", negate = TRUE)),
         any_of(var_indep_01),
         vek, denna_davka_gk_mg_kg_bw, kreatinin_umol_l, bcm, pohlavi,
         podtyp_nemoci_zjednoduseny, jo_1, anti_hmgcr, projekt_id, poradie_vysetrenia) |>
  mutate(across(.cols = c("ast", "alt", "ck", "crp", "haq", "ld", "mitax", "myoact", "myoglobin"),
                function(x) log(x + 1)),
         jo_1       = as.factor(jo_1),
         anti_hmgcr = as.factor(anti_hmgcr)) |>
  pivot_longer(
    cols     = any_of(var_dep_01 |> str_subset("odpoved_na_terapii_m0_vs_m6|kreatinin_umol_l", negate = TRUE)),
    names_to = "var_dep_name", values_to = "var_dep_value") |>
  pivot_longer(
    cols     = any_of(var_indep_01),
    names_to = "var_indep_name", values_to = "var_indep_value") |>
  group_by(var_dep_name, var_indep_name) |>
  nest() |>
  mutate(
    mod_without        = furrr::future_map(data, ~model_min_2(.x, filter_na = FALSE),                              .options = furrr_options(seed = TRUE)),
    mod_covar_gk       = furrr::future_map(data, ~model_covar_2(.x, "denna_davka_gk_mg_kg_bw", filter_na = FALSE), .options = furrr_options(seed = TRUE)),
    mod_covar_age      = furrr::future_map(data, ~model_covar_2(.x, "vek",                     filter_na = FALSE), .options = furrr_options(seed = TRUE)),
    mod_covar_bcm      = furrr::future_map(data, ~model_covar_2(.x, "bcm",                     filter_na = FALSE), .options = furrr_options(seed = TRUE)),
    mod_covar_creatine = furrr::future_map(data, ~model_covar_2(.x, "kreatinin_umol_l",        filter_na = FALSE), .options = furrr_options(seed = TRUE))
  ) |>
  ungroup()

future::plan(sequential); gc()
toc(); beep(4)

## tables ----
extract_tab <- function(nested, col, suffix) {
  nested |>
    select(var_dep_name, var_indep_name, {{ col }}) |>
    mutate(tbl = map({{ col }}, ~.x$table)) |>
    unnest(tbl) |>
    select(-{{ col }}) |>
    arrange(p_val_covar) |>
    rename_with(~ paste0(.x, suffix),
                .cols = c("omega_sq_p", "r2_mod", "beta_covar", "std_beta_covar",
                          "lower_covar", "upper_covar", "p_val_covar"))
}

tab_without  <- extract_tab(res_mixMod_covar_02_rev, mod_without,        "_without")
tab_age      <- extract_tab(res_mixMod_covar_02_rev, mod_covar_age,      "_age")
tab_gk       <- extract_tab(res_mixMod_covar_02_rev, mod_covar_gk,       "_gk")
tab_bcm      <- extract_tab(res_mixMod_covar_02_rev, mod_covar_bcm,      "_bcm")
tab_creatine <- extract_tab(res_mixMod_covar_02_rev, mod_covar_creatine, "_creatine")

res_covar_02_rev_tab <- tab_without |>
  full_join(tab_age,      by = c("var_dep_name", "var_indep_name")) |>
  full_join(tab_gk,       by = c("var_dep_name", "var_indep_name")) |>
  full_join(tab_bcm,      by = c("var_dep_name", "var_indep_name")) |>
  full_join(tab_creatine, by = c("var_dep_name", "var_indep_name"))

## export table ----
export(res_covar_02_rev_tab,
       paste0("output/tables/", format(Sys.Date(), "%y%m%d"), "_covariates_revize_01.xlsx"))


# **********************************************************************
# 5. Post-hoc subgrouping ----
# **********************************************************************

## helper: build subgroup tibble for one covariate ----
build_subgroup <- function(nested, tab_col, covar_var) {
  nested |>
    select(var_dep_name, var_indep_name, data) |>
    inner_join(
      tab_col |> filter(if_any(last_col(), ~ . < 0.06)) |> select(var_dep_name, var_indep_name),
      by = c("var_dep_name", "var_indep_name")
    ) |>
    mutate(
      fig_sex   = pmap(list(data = data, ylab = var_dep_name, xlab = var_indep_name,
                            covar = covar_var, covar_spec = "pohlavi"),
                       plot_lmer_grid_basic_02),
      fig_sub   = pmap(list(data = data, ylab = var_dep_name, xlab = var_indep_name,
                            covar = covar_var, covar_spec = "podtyp_nemoci_zjednoduseny"),
                       plot_lmer_grid_basic_02),
      fig_jo    = pmap(list(data = data, ylab = var_dep_name, xlab = var_indep_name,
                            covar = covar_var, covar_spec = "jo_1"),
                       plot_lmer_grid_basic_02),
      fig_hmgcr = pmap(list(data = data, ylab = var_dep_name, xlab = var_indep_name,
                            covar = covar_var, covar_spec = "anti_hmgcr"),
                       plot_lmer_grid_basic_02)
    )
}

## helper: save one figure column to PDF ----
save_subgroup_pdf <- function(tbl, fig_col, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  pdf(path, width = 14, height = 8)
  walk(tbl[[fig_col]], print)
  dev.off()
  invisible(path)
}

date_prefix <- format(Sys.Date(), "%y%m%d")

## age ----
sub_age <- build_subgroup(res_mixMod_covar_02_rev, tab_age, "vek")

save_subgroup_pdf(sub_age, "fig_sex",   paste0("output/figures/", date_prefix, "_revize_age_sex.pdf"))
save_subgroup_pdf(sub_age, "fig_sub",   paste0("output/figures/", date_prefix, "_revize_age_sub.pdf"))
save_subgroup_pdf(sub_age, "fig_jo",    paste0("output/figures/", date_prefix, "_revize_age_jo1.pdf"))
save_subgroup_pdf(sub_age, "fig_hmgcr", paste0("output/figures/", date_prefix, "_revize_age_hmgcr.pdf"))

## gk ----
sub_gk <- build_subgroup(res_mixMod_covar_02_rev, tab_gk, "denna_davka_gk_mg_kg_bw")

save_subgroup_pdf(sub_gk, "fig_sex",   paste0("output/figures/", date_prefix, "_revize_gk_sex.pdf"))
save_subgroup_pdf(sub_gk, "fig_sub",   paste0("output/figures/", date_prefix, "_revize_gk_sub.pdf"))
save_subgroup_pdf(sub_gk, "fig_jo",    paste0("output/figures/", date_prefix, "_revize_gk_jo1.pdf"))
save_subgroup_pdf(sub_gk, "fig_hmgcr", paste0("output/figures/", date_prefix, "_revize_gk_hmgcr.pdf"))

## bcm ----
sub_bcm <- build_subgroup(res_mixMod_covar_02_rev, tab_bcm, "bcm")

save_subgroup_pdf(sub_bcm, "fig_sex",   paste0("output/figures/", date_prefix, "_revize_bcm_sex.pdf"))
save_subgroup_pdf(sub_bcm, "fig_sub",   paste0("output/figures/", date_prefix, "_revize_bcm_sub.pdf"))
save_subgroup_pdf(sub_bcm, "fig_jo",    paste0("output/figures/", date_prefix, "_revize_bcm_jo1.pdf"))
save_subgroup_pdf(sub_bcm, "fig_hmgcr", paste0("output/figures/", date_prefix, "_revize_bcm_hmgcr.pdf"))

## creatine ----
sub_creatine <- build_subgroup(res_mixMod_covar_02_rev, tab_creatine, "kreatinin_umol_l")

save_subgroup_pdf(sub_creatine, "fig_sex",   paste0("output/figures/", date_prefix, "_revize_creatine_sex.pdf"))
save_subgroup_pdf(sub_creatine, "fig_sub",   paste0("output/figures/", date_prefix, "_revize_creatine_sub.pdf"))
save_subgroup_pdf(sub_creatine, "fig_jo",    paste0("output/figures/", date_prefix, "_revize_creatine_jo1.pdf"))
save_subgroup_pdf(sub_creatine, "fig_hmgcr", paste0("output/figures/", date_prefix, "_revize_creatine_hmgcr.pdf"))


# **********************************************************************
# 6. PLS-DA — M0, MDA kategorie (revidovaná data) ----
# **********************************************************************

library(mixOmics)

## variables to include in PLS-DA — taken from 'abbreviation' column of MDA sheet ----
pls_vars <- tryCatch({
  mda_raw <- import(path_clin_01, which = "260507_MDAcategories ")
  names(mda_raw) <- janitor::make_clean_names(names(mda_raw))
  janitor::make_clean_names(mda_raw$abbreviation)
}, error = function(e) {
  message("Sheet '260507_MDAcategories': ", conditionMessage(e), " — no column filtering applied.")
  character(0)
})

## data ----
#
# VÝBĚR PROMĚNNÝCH A OŠETŘENÍ CHYBĚJÍCÍCH HODNOT
# -----------------------------------------------
# Vstup: d03_rev filtrovaný na M0 → 57 pacientů, všichni mají mda_categories.
#
# Vyřazené proměnné (záměrně, před PLS-DA):
#   - muscle_disease_activity  : zdrojová kontinuální proměnná pro mda_categories
#                                (data leakage — triviálně by separovala skupiny)
#   - odpoved_na_terapii_m0_vs_m6 : outcome léčby, v M0 neznámý
#
# Log-transformace (pravostranně zešikmené klinické markery):
#   ast, alt, ck, crp, haq, ld, mitax, myoact, myoglobin → log(x + 1)
#
# Problém s chybějícími hodnotami:
#   Bez dalšího ošetření by complete.cases() ponechal jen 22 ze 57 pacientů.
#   Hlavní zdroje NA (ze 57 pacientů):
#     sf36_mcs / sf36_pcs : 16 NA (28 %)
#     sf36_gh              : 15 NA (26 %)
#     sf36_vt/re/mh        : 13 NA (23 %)
#     sf36_pf/rp/bp/sf     : 12 NA (21 %)
#     anti_hmgcr           : 10 NA (18 %)
#     fi_2, borg10         :  8 NA (14 %)
#     ld                   :  6 NA (11 %)
#     fstl3                :  5 NA  (9 %)
#   Různí pacienti mají chybějící hodnoty v různých proměnných → kombinace
#   vyřadí 35 pacientů.
#
# Řešení:
#   1. Práh 40 % — odstraní proměnné se strukturální absencí (v tomto datasetu
#      žádná proměnná nepřekračuje 40 %, jde o pojistku pro budoucí datové verze).
#   2. Jednorázová mice imputace (PMM, m = 1, maxit = 10, seed = 42) doplní
#      zbývající NA → finální dataset má 57 kompletních řádků.
#
d_pls_m0_rev <- d03_rev |>
  filter(poradie_vysetrenia == "M0") |>
  mutate(
    across(.cols = c("ast", "alt", "ck", "crp", "haq", "ld", "mitax", "myoact", "myoglobin"),
           function(x) log(x + 1)),
    jo_1       = as.factor(jo_1),
    anti_hmgcr = as.factor(anti_hmgcr)
  ) |>
  select(projekt_id, mda_categories,
         pohlavi, jo_1, anti_hmgcr, podtyp_nemoci_zjednoduseny,
         any_of(var_dep_01 |> str_subset("odpoved_na_terapii_m0_vs_m6|muscle_disease_activity", negate = TRUE)),
         any_of(var_indep_01)) |>
  filter(!is.na(mda_categories)) |>
  tibble::column_to_rownames("projekt_id")

# práh 40 % — pojistka proti strukturální absenci, neodstraní žádnou proměnnou
# v aktuálních datech (max NA rate = 28 % u sf36_mcs / sf36_pcs)
na_rate   <- colMeans(is.na(d_pls_m0_rev |> select(-mda_categories)))
keep_cols <- names(na_rate[na_rate <= 0.40])
d_pls_m0_rev <- d_pls_m0_rev |> select(mda_categories, all_of(keep_cols))

# mice PMM imputace — zachová všech 57 pacientů
set.seed(42)
d_pls_m0_rev_imp <- mice::mice(d_pls_m0_rev, method = "pmm",
                                m = 1, maxit = 10, printFlag = FALSE) |>
  mice::complete(1)

## function ----
run_plsda_rev <- function(df, yvar, prefix,
                          outdir = "output/figures/cluster analyses/",
                          ncomp = 2,
                          sel_vars = character(0)) {
  message("Running PLS-DA for: ", yvar)

  df2 <- as_tibble(df) |> filter(!is.na(.data[[yvar]]))

  predictors <- df2 |> select(-all_of(yvar))

  if (length(sel_vars) > 0) {
    predictors <- predictors |> select(any_of(as.character(sel_vars)))
  }

  complete_idx <- complete.cases(predictors)
  predictors   <- predictors |> slice(which(complete_idx))
  Y            <- droplevels(df2[[yvar]][complete_idx])

  X <- predictors |>
    mutate(across(where(is.factor), as.integer)) |>
    as.matrix()
  if (nrow(X) != length(Y)) stop("X/Y mismatch: nrow(X)=", nrow(X), " length(Y)=", length(Y))

  pls_mod <- plsda(X, Y, ncomp = ncomp, scale = TRUE)

  ellipse_flag <- !any(table(Y) <= 2)
  if (!ellipse_flag) message("Ellipse disabled: some groups have ≤ 2 samples")

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  tiff(paste0(outdir, prefix, "_", yvar, "_group.tiff"),
       compression = "lzw", res = 300, width = 1600, height = 1200)
  plotIndiv(pls_mod, comp = c(1, 2), group = Y, ellipse = ellipse_flag,
            legend = TRUE, pch = 19, ind.names = FALSE,
            title = paste0("PLS-DA — ", yvar))
  dev.off()

  tiff(paste0(outdir, prefix, "_", yvar, "_loadings.tiff"),
       compression = "lzw", res = 300, width = 1600, height = 1200)
  plotLoadings(pls_mod, contrib = "max", method = "median",
               title = paste0("Loadings — ", yvar))
  dev.off()

  list(model = pls_mod, Y = Y)
}

## run ----
pls_m0_rev_mda <- run_plsda_rev(
  df       = d_pls_m0_rev_imp,
  yvar     = "mda_categories",
  prefix   = paste0(date_prefix, "_revize_m0_pls"),
  sel_vars = pls_vars
)


# **********************************************************************
# 7. Save RData for report ----
# **********************************************************************

rdata_path <- paste0("reports/", date_prefix, "_covariates_revize_01.RData")

save(res_covar_02_rev_tab,
     sub_age, sub_gk, sub_bcm, sub_creatine,
     pls_m0_rev_mda,
     d04_sel2_rev,
     file = rdata_path)

message("Done. RData saved to: ", rdata_path)


# **********************************************************************
# 8. Render report ----
# **********************************************************************

rmarkdown::render(
  input       = "reports/260511_iim_revize_02_report.Rmd",
  output_file = paste0(date_prefix, "_iim_revize_02_report.html"),
  output_dir  = "reports",
  quiet       = TRUE
)
message("Report rendered: reports/", date_prefix, "_iim_revize_02_report.html")
beep(4)
