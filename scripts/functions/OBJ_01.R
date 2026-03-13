# *******************************************************
# 1. variables ----
# *******************************************************
var_indep_01 <- c("mstn", "fst", "fstl3", "akta")
var_dep_01 <- c("mmt8_total", "mmt10_total", "fi_2", "borg10", "haq", 
                "sf36_mcs", "sf36_pf", "sf36_rp", "sf36_bp", "sf36_gh", "sf36_vt",
                "sf36_sf", "sf36_re", "sf36_mh", "sf36_pcs", "crp", "ast", "alt", "ck", "ld", 
                "kreatinin_umol_l", "mitax", "myoglobin","myoact", 
                "physician_activity_vas", "muscle_disease_activity", 
                "odpoved_na_terapii_m0_vs_m6")
var_dep_02 <- c("mmt8_total", "mmt10_total", "fi_2")
var_dep_03 <- c("mmt8_total", "mmt10_total", "fi_2", "borg10", "haq", 
                "sf36_mcs", "sf36_pf", "sf36_rp", "sf36_bp", "sf36_gh", "sf36_vt",
                "sf36_sf", "sf36_re", "sf36_mh", "sf36_pcs", "crp", "ast", "alt", "ck", "ld", 
                "mitax", "myoglobin","myoact", 
                "physician_activity_vas", "muscle_disease_activity", 
                "odpoved_na_terapii_m0_vs_m6")
sel_val_01 <- import("data/processed/selected_colnames_01.xlsx")

sel_val_02 <- c("projekt ID", "Poradie vyšetrenia", "Věk", "Pohlaví", "Trvání nemoci", 
                "Podtyp nemoci_zjednodušený", "MITAX", "MYOACT", 
                "Patient activity VAS", "Physician activity VAS", 
                "Muscle Disease Activity", "Extramuscular Global Assessment", 
                "MMT8 total", "MMT10 total", "FI-2", "Borg10", "HAQ", "SF36-PCS", 
                "SF36-MCS", "odpověď na terapii M0 vs M6", "CRP", "AST", "ALT", 
                "CK", "LD", "myoglobin", "kreatinin (umol/l)", "Jo-1", 
                "anti HMGCR", "denná dávka GK (mg/kg BW)", "BMI", "Body Fat in %", 
                "Lean Body Mass", "BCM", "ECM/BCM Index", "% cell share")

sel_val_03 <- c("projekt ID", "Poradie vyšetrenia", "Věk", "Pohlaví", "Trvání nemoci", 
                "Podtyp nemoci_zjednodušený", "MITAX", "MYOACT", 
                "Patient activity VAS", "Physician activity VAS", 
                "Muscle Disease Activity", "Extramuscular Global Assessment",
                "MDA categories", "MMT8 categories",
                "MMT8 total", "MMT10 total", "FI-2", "Borg10", "HAQ", "SF36-PCS", 
                "SF36-MCS", "odpověď na terapii M0 vs M6", "CRP", "AST", "ALT", 
                "CK", "LD", "myoglobin", "kreatinin (umol/l)", "Jo-1", 
                "anti HMGCR", "denná dávka GK (mg/kg BW)", "BMI", "Body Fat in %", 
                "Lean Body Mass", "BCM", "ECM/BCM Index", "% cell share")

ranges <- c("t20_t40", "t40_t100", "t100_t400", "t400_t8000")
bouts  <- c("b1", "b5", "b10")

sel_time_01 <- c("p9931_enmo", "m5_enmo", "ig_intercept_enmo", "ig_gradient_enmo",
                 "enmo_full_recording_mean", "p9861_enmo", "p9792_enmo", "p9688_enmo",
                 "p9583_enmo", "p9375_enmo", "p9167_enmo")

sel_activ_01 <- c("t100_time_min", "moderate", "light", "inactivity",
                  "t40_time_min", "t20_time_min", "t30_time_min")
sel_activ_01b <- c("mvpa_t100_time_min", "moderate", "light", "inactivity",
                   "t40_time_min", "t20_time_min", "t30_time_min")

# Safe worker cap for heavy parallel sections (important for background jobs).
safe_future_workers_01 <- as.integer(getOption(
  "myo_ncores",
  max(1L, min(4L, future::availableCores() - 1L))
))
if (!is.finite(safe_future_workers_01) || is.na(safe_future_workers_01)) {
  safe_future_workers_01 <- max(1L, min(4L, future::availableCores() - 1L))
}
safe_future_workers_01 <- max(1L, safe_future_workers_01)
sel_clinics_02 <- c("mmt8_total", "fi2", "muscle_activity", "pulmonary_activity", 
                    "borg10", "haq", "physician_vas",  "phase_angle", "ast", "b_m_rate",
                    "creatinine", "gc_dose_mg_day")

sel_clinics_01 <- c("mmt8_total", "fi2", "muscle_activity", "pulmonary_activity", 
                    "borg10", "haq", "bcm",  "phase_angle", "ast", "b_m_rate",
                    "creatinine", "gc_dose_mg_day")

sel_time_02 <- c("t100_time_min", "moderate", "light", "inactivity",
                 "t40_time_min", "t20_time_min", "t30_time_min")
sel_time_02b <- c("mvpa_t100_time_min", "moderate", "light", "inactivity",
                 "t40_time_min", "t20_time_min", "t30_time_min")

sel_activ_02 <- c("p9931_enmo", "m5_enmo", "ig_intercept_enmo", "ig_gradient_enmo",
                  "acc_day_spt_wei", "p9861_enmo", "p9792_enmo", "p9688_enmo",
                  "p9583_enmo", "p9375_enmo", "p9167_enmo")

folder_path <- "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Data/"
d17_myokines <- import(paste0(folder_path, "Treatment Response_sérum,expr_pre STAT .xlsx")) |> 
  select(!contains("pg")) |> 
  rename(immet_id = `projekt ID`) |> 
  mutate(immet_id = str_remove_all(immet_id, "_"))

sel_myokines_01 <- d17_myokines[,-1] |> names()

sel_tis_d18_data_01 <- c(
  "extramuscular_global_assessment", "pt_vas", "ph_vas", "mmt8", "haq", "ast", "alt", "ck", "ld",
  "borg", "da_muscle", "crp", "myoglobin",
  "phase_angle", "bcm", "sf36_pcs",  "sf36_mcs", "sf36pf", "sf36rp", "sf36bp",
  "sf36gh", "sf36vt", "sf36sf", "sf36re", "sf36mf"
)

sel_tis_d18_data_02 <- c("extramuscular_global_assessment", "pt_vas", "ph_vas", "mmt8", "haq",
                         "ast", "alt", "ck", "ld")

# *******************************************************
# 2. data ---- 
# *******************************************************
## importing ----
# folder_path <- "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Data/"

path_myo_01 <- file_to_load(phrase = "Myokiny", 
                            folder_path = folder_path)

d01_myokiny <- import(path_myo_01) |>
  # relocate(`věk skupiny`, .before = `věk k M0`) |> 
  select(!starts_with("d")) |> 
  mutate(across(where(is.character), ~na_if(.x, "x")),
         across(-c(1:4), ~as.numeric(.x))) |>
  pivot_longer(
    cols = matches("^[A-Za-z0-9]+ M\\d+$"),       
    names_to = c("protein", "poradie_vysetrenia"),
    names_sep = " ",
    values_to = "concentration"
  ) |>
  pivot_wider(
    names_from = protein,
    values_from = concentration
  ) |> 
  janitor::clean_names() |> 
  mutate(poradie_vysetrenia = factor(poradie_vysetrenia,
                                     levels = c("M0", "M3", "M6", "M18")),
         bbm_id = factor(bbm_id),
         pohlavi = factor(pohlavi))

d01_myokiny_2 <- import("F:/Analysis/Vernerová Lucia/myo_2025_iim/data/processed/Myokiny sérum_pro stat_processed.xlsx") |> 
  mutate(immet_id = str_remove_all(immet_id, "_")) |> 
  mice::futuremice(m = 40, maxit = 5, method = "pmm", 
                   n.core = safe_future_workers_01,
                   parallelseed = 123,            # seed pro paralelní běh (doporučeno)
                   future.plan = "sequential",
                   ridge = 1e-02,
                   remove.collinear = TRUE,
                   remove.constant = TRUE
  ) |> 
  mice::complete() |> 
  dplyr::mutate(across(where(is.character),readr::parse_number))

path_clin_01 <- file_to_load(phrase = "IMMET_terapie", 
                            folder_path = folder_path)


d02_clinics <- import(path_clin_01 |> tail(1)) |>
  mutate(across(where(is.character), ~na_if(.x, "x"))) |> 
  clean_names() |> 
  tidyr::separate(
    col = bbm_id, 
    sep = "_", 
    into = c("bbm_id", "bbm_collection"),
    extra = "merge"
  ) |> 
  mutate(
    bbm_id = factor(bbm_id),
    datum_vysetreni = if_else(str_detect(datum_vysetreni, "^x"), 
                              as.Date(NA), 
                              openxlsx::convertToDate(datum_vysetreni)
                              ),
    datum_biopsie = openxlsx::convertToDate(datum_biopsie),
    date_of_bioimpedance = openxlsx::convertToDate(date_of_bioimpedance),
    datum_nar = as.Date(datum_nar),
    bbm_collection = factor(bbm_collection),
    poradie_vysetrenia = factor(poradie_vysetrenia,
                                levels = c("M0", "M3", "M6", "M18")),
    pohlavi = factor(pohlavi),
    odpoved_na_terapii_m0_vs_m6 = factor(odpoved_na_terapii_m0_vs_m6,
                                         levels = c("No", "Minimal", "Moderate", "Major")),
    odpoved_na_terapii_m0_vs_m18 = factor((odpoved_na_terapii_m0_vs_m18)),
    podtyp_nemoci_zjednoduseny = factor(podtyp_nemoci_zjednoduseny),
    mda_categories = factor(mda_categories),
    mmt8_categories = factor(mmt8_categories),
    kurak_nejmene_100_cigaret = factor(kurak_nejmene_100_cigaret)
  ) |> 
  mutate(across(where(is_01_col), as.factor)) |> 
  mutate(across((vek:last_col()) & where(is.character), as.numeric)) |> 
  group_by(bbm_id) |> 
  fill(kurak_nejmene_100_cigaret) |> 
  ungroup()

var_name_01_mda <- import(path_clin_01 |> tail(1), which = "251118_MDAcategories") |> 
  mutate(abbreviation = janitor::make_clean_names(abbreviation))

var_name_01_mmt8 <- import(path_clin_01 |> tail(1), which = "251118_MMT8categories") |> 
  mutate(abbreviation = janitor::make_clean_names(abbreviation))

var_name <- bind_rows(var_name_01_mda, var_name_01_mmt8) |> distinct()

rename_vec <- var_name$abbreviation
names(rename_vec) <- var_name$name_en


d05_norm <- import("R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Data/MSTN normalizácia na kreat, BCM.xlsx") |> 
  mutate(across(where(is.character), ~na_if(.x, "x"))) |> 
  clean_names() |> 
  pivot_longer(
    # only columns ending with _m0/_m3/_m6/_m18 (case-insensitive), allow stray spaces
    cols = matches("_(?:m0|m3|m6|m18)\\s*$", ignore.case = TRUE),
    names_to = c("protein", "poradie_vysetrenia"),
    names_pattern = "^(.*)_(m0|m3|m6|m18)\\s*$",
    values_to = "concentration",
    values_drop_na = TRUE
  ) |>
  mutate(poradie_vysetrenia = toupper(poradie_vysetrenia),
         bbm_id = factor(bbm_id)) |> 
  pivot_wider(
    names_from = protein,
    values_from = concentration
  ) |> 
  janitor::clean_names() |> 
  mutate(poradie_vysetrenia = factor(poradie_vysetrenia,
                                     levels = c("M0", "M3", "M6", "M18")))
  
pre_filt_01 <- d01_myokiny |> 
  select(maju_klinicke_data_pre_odfiltrovanie) |> 
  na.omit() |> 
  distinct() |> 
  pull()

### ggir ----
# *******************************************************
# file_ggir_01 <- "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Geneactiv/merged_table_220126.xlsx"
file_ggir_01 <- "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/GENEActiv/GGIR/merged_table_final_pro_zaverecnou_zpravu_AZV_100226.xlsx"

d10_legend <- rio::import(file_ggir_01,
                          which="legend")
d10_ostatni <- rio::import(file_ggir_01,
                           which="ostatni") |> 
  rename(any_of(setNames(d10_legend$names_orig, d10_legend$names_new))) |> 
  filter(!is.na(file_name)) |> 
  mutate(
    across(-1, ~ as.numeric(.x))
  ) |> 
  separate(
    file_name,
    into = c("sample", "exam"),
    sep = "_",
    extra = "drop",
    fill = "right"
  ) |> 
  distinct() |> 
  mutate(exam = factor(exam, levels = c("M0", "M3", "M6", "M18", "M30", "M48")),
         sample = str_replace_all(sample, "57IMMET", "IMMET57"),
         sample = as.factor(sample))
d10_cas <- rio::import(file_ggir_01,
                           which="cas") |> 
  janitor::remove_empty("cols") |> 
  rename(any_of(setNames(d10_legend$names_orig, d10_legend$names_new))) |> 
  filter(!is.na(file_name)) |> 
  mutate(
    across(-1, ~ as.numeric(.x))
  ) |> 
  separate(
    file_name,
    into = c("sample", "exam"),
    sep = "_",
    extra = "drop",
    fill = "right"
  ) |> 
  mutate(exam = factor(exam, levels = c("M0", "M3", "M6", "M18", "M30", "M48")),
         sample = str_replace_all(sample, "57IMMET", "IMMET57"),
         sample = as.factor(sample),
         moderate = mvpa_t100_time_min - t400_time_min,
         light = t40_time_min - mvpa_t100_time_min - t400_time_min,
         inactivity = t20_time_min - t40_time_min - mvpa_t100_time_min - t400_time_min)

d10_cas_b <- rio::import(file_ggir_01,
                       which="cas") |> 
  janitor::remove_empty("cols") |> 
  rename(any_of(setNames(d10_legend$names_orig, d10_legend$names_new))) |> 
  filter(!is.na(file_name)) |> 
  mutate(
    across(-1, ~ as.numeric(.x))
  ) |> 
  separate(
    file_name,
    into = c("sample", "exam"),
    sep = "_",
    extra = "drop",
    fill = "right"
  ) |> 
  mutate(exam = factor(exam, levels = c("M0", "M3", "M6", "M18", "M30", "M48")),
         sample = str_replace_all(sample, "57IMMET", "IMMET57"),
         sample = as.factor(sample),
         vigorous = t400plus_time_min,
         moderate = t100_t400_time_min,
         light = t40_t100_time_min,
         inactivity_c = t30_t40_time_min,
         inactivity_b = t20_t30_time_min,
         inactivity_a = t0_t20_time_min)


d10_aktivita <- rio::import(file_ggir_01,
                           which="aktivita")|> 
  rename(any_of(setNames(d10_legend$names_orig, d10_legend$names_new))) |> 
  separate(
    file_name,
    into = c("sample", "exam"),
    sep = "_",
    extra = "drop",
    fill = "right"
  ) |> 
  mutate(exam = factor(exam, levels = c("M0", "M3", "M6", "M18", "M30", "M48")),
         sample = str_replace_all(sample, "57IMMET", "IMMET57"),
         sample = as.factor(sample)) |> 
  left_join(d10_ostatni |> select(-r2_enmo))

d11_legend <- rio::import("R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Geneactiv/pro závěrečnou zprávu AZV/Klinicka_tabulka_13012026.xlsx",
                          which="legend")
d11_clinics <- rio::import("R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Geneactiv/pro závěrečnou zprávu AZV/Klinicka_tabulka_13012026.xlsx",
                            which="data") |> 
  rename(any_of(setNames(d11_legend$names_orig, d11_legend$names_new))) |> 
  select(-contains("birth"), -contains("name"), -contains("date"), -notes,
         -sample_id_bbm) |> 
  mutate(
    disease_subtype = str_replace_all(disease_subtype, "^(1|7)$", NA_character_),
    across(
      everything(),
      \(col) if (is.numeric(col)) col else na_if(as.character(col), "x")
    ),
    immet_id = str_replace_all(immet_id, "_","")
  ) |> 
  relocate(immet_id, exam_order, disease_subtype, contains("response"),
           smoker) |> 
  group_by(immet_id) |>
  tidyr::fill(response_m0_m6, response_m0_m18, smoker, .direction = "downup") |>
  ungroup() |> 
  mutate(
    across(
      where(\(x) {
        x_num <- suppressWarnings(as.numeric(as.character(x)))
        !all(is.na(x_num)) && is_01_col(x_num)
      }),
      \(x) factor(as.numeric(as.character(x)), levels = c(0, 1))
    ),
    across(.cols = immet_id:smoker, as.factor),
    across(where(is.character),readr::parse_number),
    response_m0_m6 = str_to_title(response_m0_m6),
    response_m0_m6 = factor(response_m0_m6, levels = c("No Improvement", "Minimal", "Moderate", 
                                                       "Major")) |> ordered()
  )

# mice
future::plan(future::multisession, workers = safe_future_workers_01)
future::nbrOfWorkers()

d11_clinics_imp <- d11_clinics |> 
  filter(exam_order %in% c("M0", "M6")) |>
  select(
    where(~ mean(is.na(.x)) <= 0.20)
  ) |> 
  mice::futuremice(m = 40, maxit = 5, method = "pmm", 
                   n.core = safe_future_workers_01,
                   parallelseed = 123,            # seed pro paralelní běh (doporučeno)
                   future.plan = "multisession",
                   ridge = 1e-02,
                   remove.collinear = TRUE,
                   remove.constant = TRUE
                   ) |> 
  mice::complete() |> 
  dplyr::mutate(across(where(is.character),readr::parse_number))
future::plan("sequential")

# **********************************************************************
### d16 - fi, mmt, haq ----
# **********************************************************************
import_list("R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/VB_FI_MMT/FI3_MMT_260126.xlsx") |> 
  (\(x) setNames(x, paste0("d16_", make.names(names(x), unique = TRUE))))() |>
  list2env(envir = .GlobalEnv)

## merging  and processing ----
# *******************************************************
d03 <- d02_clinics |> 
  full_join(d01_myokiny |> 
              select(-maju_klinicke_data_pre_odfiltrovanie),
            by = c("projekt_id",
                   "poradie_vysetrenia",
                   "bbm_id",
                   "pohlavi")) |> 
  filter(projekt_id %in% pre_filt_01) |> 
  group_by(projekt_id) |> 
  fill(podtyp_nemoci_zjednoduseny) |> 
  fill(odpoved_na_terapii_m0_vs_m6, .direction = "downup") |> 
  ungroup()

d06_norm_merge <- d05_norm |> 
  full_join(d01_myokiny |> 
              select(-maju_klinicke_data_pre_odfiltrovanie),
            by = c("projekt_id",
                   "poradie_vysetrenia",
                   "bbm_id")) |> 
  left_join(d02_clinics |> select(bbm_id, podtyp_nemoci_zjednoduseny, 
                                  odpoved_na_terapii_m0_vs_m6, vek, poradie_vysetrenia,
                                  gk, davka_gk_mg_den,denna_davka_gk_mg_kg_bw, 
                                  all_of(var_dep_01)),
            by = c("bbm_id",
                   "poradie_vysetrenia")) |> 
  group_by(projekt_id) |> 
  fill(podtyp_nemoci_zjednoduseny) |> 
  fill(odpoved_na_terapii_m0_vs_m6, .direction = "downup") |> 
  ungroup()

d16_statistika_mmt <- d16_statistika_mmt |>   
  mutate(immet_id = str_replace_all(immet_id, "_",""),
         across(
           everything(),
           ~ na_if(as.character(.x), "x"))) |> 
  pivot_longer(cols = -immet_id,
               names_to = "var_name",
               values_to = "var_value") |> 
  separate(
    var_name,
    into = c("exam", "var_name"),
    sep = "_",
    extra = "merge"
  ) |> 
  mutate(var_value = as.numeric(var_value)) |> 
  distinct() |> 
  pivot_wider(names_from = var_name,
              values_from = var_value) |> 
  rename(any_of(setNames(d16_legenda[[1]],d16_legenda[[2]])))


d16_statistika_fi <- d16_statistika_fi |> 
  filter(!is.na(immet_id)) |> select(-2) |> 
  mutate(immet_id = str_replace_all(immet_id, "_",""),
         across(
           everything(),
           ~ na_if(as.character(.x), "x"))) |> 
  pivot_longer(cols = -immet_id,
               names_to = "var_name",
               values_to = "var_value") |> 
  separate(
    var_name,
    into = c("exam", "var_name"),
    sep = "_",
    extra = "merge"
  ) |> 
  mutate(var_value = as.numeric(var_value)) |> 
  distinct() |> 
  pivot_wider(names_from = var_name,
              values_from = var_value) |> 
  rename(any_of(setNames(d16_legenda[[1]],d16_legenda[[2]])))

d16_haq <- d16_haq |> 
  mutate(immet_id = str_replace_all(immet_id, "_",""))

age_m0 <- d11_clinics |>
  filter(exam_order == "M0") |>
  select(immet_id, age) |>
  rename(age_m0 = age)


## cluster analyses ----
path_response_01 <- file_to_load(phrase = "Treatment Response", 
                                 folder_path = folder_path)

d07 <- import(path_response_01 |> tail(1)) |> 
  janitor::clean_names() |> 
  mutate(across(where(is.character), ~na_if(.x, "x")),
         projekt_id = as.factor(projekt_id),
         across(where(is.character), as.numeric)) |> 
  droplevels()

### imputation ----
d08_imp <- d02_clinics |> 
  filter(projekt_id != "IMMET_40") |>     # má jenom jedno měření, v M0
  # select(all_of(janitor::make_clean_names(sel_val_02))) |> 
  select(all_of(janitor::make_clean_names(sel_val_03))) |> # oproti sel_val_doplněno o MDA a MMT8 categories
  select(-odpoved_na_terapii_m0_vs_m6) |> # mohl bych odstranit, ale kvůli komplitnímu výpisu všech proměnných to udělám takhle
  group_by(projekt_id) |> 
  fill(vek, .direction = "downup") |> 
  fill(jo_1, .direction = "downup") |> 
  fill(anti_hmgcr, .direction = "downup") |> 
  ungroup() |> 
  filter(poradie_vysetrenia %in% c("M0", "M6")) |> 
  relocate(projekt_id, poradie_vysetrenia, podtyp_nemoci_zjednoduseny, 
           pohlavi, jo_1, anti_hmgcr) |> 
  filter(!is.na(podtyp_nemoci_zjednoduseny)) |> 
  left_join(d07, by = "projekt_id") |> 
  mutate(across(.cols = projekt_id:anti_hmgcr, ~as.factor(.x))) |> 
  mice::futuremice(method = "pmm", m = 20, n.core = ncores, ridge = 1e-02, remove.collinear = TRUE, remove.constant = TRUE, future.plan = "sequential") |>
  complete() |>
  left_join(d02_clinics |> 
              select(projekt_id,odpoved_na_terapii_m0_vs_m6) |> 
              group_by(projekt_id) |> 
              fill(odpoved_na_terapii_m0_vs_m6, .direction = "downup") |> 
              ungroup() |> 
              unique(), 
            by = "projekt_id") |> 
  relocate(odpoved_na_terapii_m0_vs_m6, mda_categories, mmt8_categories, .before =  pohlavi) 

### M0 ----
d08_m0 <- d08_imp |> 
  filter(poradie_vysetrenia == "M0") |> 
  tibble::column_to_rownames("projekt_id") |>
  select(-poradie_vysetrenia) |> 
  mutate(
    across(where(is.numeric), ~scale(.x)),
    across(where(is.factor),
           ~ forcats::fct_relabel(.x, function(lvl) paste0(cur_column(), "_", lvl)))) |> 
  droplevels()

# PLS-DA 
d08_m0_pls <- d08_imp |> 
  filter(poradie_vysetrenia == "M0") |> 
  tibble::column_to_rownames("projekt_id") |>
  select(-poradie_vysetrenia)

### delta ----
# measured once
sel_remove_01 <- c("mstn", "fst", "fstl3", "acvr1b", "acvr2b", "acvr2a", "smad2",
                   "smad3", "smad4", "foxo1", "trim63", "fbxo32", "akta")

d08_delta <- d08_imp |> 
  select(-mda_categories, -mmt8_categories, -contains(sel_remove_01)) |> 
  left_join(
    d08_imp |> 
      filter(poradie_vysetrenia == "M6") |> 
      select(projekt_id, mda_categories, mmt8_categories) |> 
      droplevels(),
    by = "projekt_id"
  ) |> 
  relocate(mda_categories, mmt8_categories, .before = pohlavi) |> 
  pivot_wider(
    names_from = poradie_vysetrenia,
    values_from = where(is.numeric)
  ) |> 
  mutate(across(ends_with("M6"), 
                .fns = ~ . - get(str_replace(cur_column(), "M6$", "M0")),
                .names = "{str_remove(.col, 'M6$')}_delta")) |>
  select(projekt_id, where(is.factor), ends_with("_delta"), 
         -odpoved_na_terapii_m0_vs_m6_delta, -vek__delta, -trvani_nemoci__delta) |> 
  left_join(d08_imp |> filter(poradie_vysetrenia == "M0") |>
              select(projekt_id, vek, trvani_nemoci) |> 
              distinct(),
            by = "projekt_id") |> 
  relocate(vek, trvani_nemoci, .before = mitax__delta) |> 
  tibble::column_to_rownames("projekt_id") |>
  mutate(
    across(where(is.numeric), ~scale(.x)),
    across(where(is.factor),
           ~ forcats::fct_relabel(.x, function(lvl) paste0(cur_column(), "_", lvl)))) |> 
  droplevels()


### PLS-DA ----
d08_delta_pls <- d08_imp |> 
  select(-mda_categories, -mmt8_categories, -contains(sel_remove_01)) |> 
  left_join(
    d08_imp |> 
      filter(poradie_vysetrenia == "M6") |> 
      select(projekt_id, mda_categories, mmt8_categories) |> 
      droplevels(),
    by = "projekt_id"
  ) |> 
  relocate(mda_categories, mmt8_categories, .before = pohlavi) |> 
  pivot_wider(
    names_from = poradie_vysetrenia,
    values_from = where(is.numeric)
  ) |> 
  mutate(across(ends_with("M6"), 
                .fns = ~ . - get(str_replace(cur_column(), "M6$", "M0")),
                .names = "{str_remove(.col, 'M6$')}_delta")) |>
  select(projekt_id, where(is.factor), ends_with("_delta"), 
         -odpoved_na_terapii_m0_vs_m6_delta, -vek__delta, -trvani_nemoci__delta) |> 
  left_join(d08_imp |> filter(poradie_vysetrenia == "M0") |>
              select(projekt_id, vek, trvani_nemoci) |> 
              distinct(),
            by = "projekt_id") |> 
  relocate(vek, trvani_nemoci, .before = mitax__delta) |> 
  tibble::column_to_rownames("projekt_id")


## cluster analyses 2 ----
### imputation ----
d09_imp <- d02_clinics |> 
  filter(projekt_id != "IMMET_40") |>     # má jenom jedno měření, v M0
  group_by(projekt_id) |> 
  fill(vek, .direction = "downup") |> 
  fill(jo_1, .direction = "downup") |> 
  fill(anti_hmgcr, .direction = "downup") |> 
  ungroup() |> 
  filter(poradie_vysetrenia %in% c("M0", "M6")) |> 
  relocate(projekt_id, poradie_vysetrenia,  
           pohlavi, jo_1, anti_hmgcr) |> 
  select(any_of(var_name$abbreviation)) |> 
  mutate(across(where(is_01_col), ~ as.factor(.x))) |> 
  futuremice_safe(method = "pmm", m = ncores, workers = ncores) |> 
  complete() |>
  relocate( mda_categories, mmt8_categories, .before =  pohlavi) 

### M0 ----
d09_m0 <- d09_imp |> 
  filter(poradie_vysetrenia == "M0") |> 
  left_join(d07, 
            by = "projekt_id") |> 
  futuremice_safe(method = "pmm", m = ncores, workers = ncores) |> 
  complete() |> 
  tibble::column_to_rownames("projekt_id") |>
  select(-poradie_vysetrenia) |> 
  mutate(
    across(where(is.numeric), ~scale(.x)),
    across(where(is.factor),
           ~ forcats::fct_relabel(.x, function(lvl) paste0(cur_column(), "_", lvl)))) |> 
  droplevels() |> 
  mutate(across(everything(), ~ if(is.matrix(.x) && ncol(.x) == 1) .x[,1] else .x))

# measured once
sel_remove_01 <- c("mstn", "fst", "fstl3", "acvr1b", "acvr2b", "acvr2a", "smad2",
                   "smad3", "smad4", "foxo1", "trim63", "fbxo32", "akta")

# PLS-DA 
d09_m0_pls <- d09_imp |> 
  filter(poradie_vysetrenia == "M0") |> 
  left_join(d09_m0 |> 
              rownames_as_column("projekt_id") |> 
              select(projekt_id, contains(sel_remove_01)) |> 
              tibble(), 
            by = "projekt_id") |> 
  tibble::column_to_rownames("projekt_id") |>
  select(-poradie_vysetrenia)

### delta ----
d09_delta <- d09_imp |>
  group_by(projekt_id) |>
  mutate(
    across(
      c(anti_hmgcr, jo_1),
      ~ ifelse(any(as.character(.x) == "1", na.rm = TRUE), "1", "0")
    )
  ) |>
  ungroup() |>
  mutate(
    across(c(anti_hmgcr, jo_1), ~ factor(.x, levels = c("0", "1")))
  )|> 
  select(-vek, -trvani_nemoci, -mda_categories, -mmt8_categories, -contains(sel_remove_01)) |> 
  left_join(
    d09_imp |> 
      filter(poradie_vysetrenia == "M6") |> 
      select(projekt_id, mda_categories, mmt8_categories) |> 
      droplevels(),
    by = "projekt_id"
  ) |> 
  relocate(mda_categories, mmt8_categories, .before = pohlavi) |> 
  pivot_wider(
    names_from = poradie_vysetrenia,
    values_from = where(is.numeric)
  ) |> 
  mutate(across(ends_with("M6"), 
                .fns = ~ . - get(str_replace(cur_column(), "M6$", "M0")),
                .names = "{str_remove(.col, 'M6$')}delta")) |>
  select(projekt_id, where(is.factor), ends_with("_delta")) |> 
  left_join(d09_imp |> filter(poradie_vysetrenia == "M0") |>
              select(projekt_id, vek, trvani_nemoci) |> 
              distinct(),
            by = "projekt_id") |> 
  relocate(vek, trvani_nemoci, .before = mitax_delta) |> 
  tibble::column_to_rownames("projekt_id") |>
  mutate(
    across(where(is.numeric), ~scale(.x)),
    across(where(is.factor),
           ~ forcats::fct_relabel(.x, function(lvl) paste0(cur_column(), "_", lvl)))) |> 
  droplevels()  |> 
  mutate(across(everything(), ~ if(is.matrix(.x) && ncol(.x) == 1) .x[,1] else .x)) |> 
  rename_with(~ str_remove(., "_delta"), everything())


### PLS-DA ----
d09_delta_pls <-  d09_imp |>
  group_by(projekt_id) |>
  mutate(
    across(
      c(anti_hmgcr, jo_1),
      ~ ifelse(any(as.character(.x) == "1", na.rm = TRUE), "1", "0")
    )
  ) |>
  ungroup() |>
  mutate(
    across(c(anti_hmgcr, jo_1), ~ factor(.x, levels = c("0", "1")))
  )|> 
  select(-vek, -trvani_nemoci, -mda_categories, -mmt8_categories, -contains(sel_remove_01)) |> 
  left_join(
    d09_imp |> 
      filter(poradie_vysetrenia == "M6") |> 
      select(projekt_id, mda_categories, mmt8_categories) |> 
      droplevels(),
    by = "projekt_id"
  ) |> 
  relocate(mda_categories, mmt8_categories, .before = pohlavi) |> 
  pivot_wider(
    names_from = poradie_vysetrenia,
    values_from = where(is.numeric)
  ) |> 
  mutate(across(ends_with("M6"), 
                .fns = ~ . - get(str_replace(cur_column(), "M6$", "M0")),
                .names = "{str_remove(.col, 'M6$')}delta")) |>
  select(projekt_id, where(is.factor), ends_with("_delta")) |> 
  left_join(d09_imp |> filter(poradie_vysetrenia == "M0") |>
              select(projekt_id, vek, trvani_nemoci) |> 
              distinct(),
            by = "projekt_id") |> 
  relocate(vek, trvani_nemoci, .before = mitax_delta) |> 
  tibble::column_to_rownames("projekt_id") |>
  droplevels()  |> 
  mutate(across(everything(), ~ if(is.matrix(.x) && ncol(.x) == 1) .x[,1] else .x)) |> 
  rename_with(~ str_remove(., "_delta"), everything()) 


## selection ----
d04_sel1 <- d03 |>
  mutate(
    odpoved_na_terapii_m0_vs_m6 =
      str_replace(odpoved_na_terapii_m0_vs_m6, "No Improvement", "No improvement")
  ) |>
  select(projekt_id, pohlavi, vek, gk, davka_gk_mg_den,denna_davka_gk_mg_kg_bw, 
         podtyp_nemoci_zjednoduseny, poradie_vysetrenia, all_of(var_dep_01), all_of(var_indep_01)) |>
  filter(poradie_vysetrenia %in% c("M0", "M3", "M6", "M18"))

d04_sel2 <- d03 |>
  mutate(
    odpoved_na_terapii_m0_vs_m6 =
      str_replace(odpoved_na_terapii_m0_vs_m6, "No Improvement", "No improvement")
  ) |>
  group_by(projekt_id) |>
  mutate(
    across(c(jo_1, anti_hmgcr),
           ~ if_else(any(. == 1, na.rm = TRUE), 1L, 0L))
  ) |>
  ungroup() |>
  select(projekt_id, pohlavi, vek, gk, davka_gk_mg_den,denna_davka_gk_mg_kg_bw, 
         podtyp_nemoci_zjednoduseny, poradie_vysetrenia, bcm, jo_1, anti_hmgcr,
         all_of(var_dep_01), all_of(var_indep_01)) |>
  filter(poradie_vysetrenia %in% c("M0", "M3", "M6", "M18"))

d04_sel1_nested <- d04_sel1 |> 
  # mutate(poradie_vysetrenia = str_remove(poradie_vysetrenia, "M") |> as.numeric()) |>
  # mutate(
  #   odpoved_na_terapii_m0_vs_m6 = recode(
  #     odpoved_na_terapii_m0_vs_m6,
  #     "No"             = 0,
  #     "No improvement" = 1,
  #     "Minimal"        = 2,
  #     "Moderate"       = 3,
  #     "Major"          = 4,
  #     .default = NA_real_
  #   )
  # ) |> 
  pivot_longer(cols = any_of(var_dep_01 |>  str_subset("^odpoved_na_terapii_m0_vs_m6$", negate = TRUE)),
               names_to = "var_dep_name",
               values_to = "var_dep_value") |> 
  pivot_longer(cols = any_of(var_indep_01),
               names_to = "var_indep_name",
               values_to = "var_indep_value") |> 
  group_by(var_dep_name, var_indep_name) |> 
  nest()

d06_norm_merge_nest <- d06_norm_merge |> 
  mutate(
    odpoved_na_terapii_m0_vs_m6 =
      str_replace(odpoved_na_terapii_m0_vs_m6, "No Improvement", "No improvement")
  ) |>
  select(projekt_id, pohlavi, vek, gk, davka_gk_mg_den,denna_davka_gk_mg_kg_bw,
         mstn_to_kreatinin, mstn_to_bcm,
         podtyp_nemoci_zjednoduseny, poradie_vysetrenia, all_of(var_dep_01)) |>
  filter(poradie_vysetrenia %in% c("M0", "M3", "M6", "M18")) |> 
  pivot_longer(cols = any_of(var_dep_01 |>  str_subset("^odpoved_na_terapii_m0_vs_m6$", negate = TRUE)),
               names_to = "var_dep_name",
               values_to = "var_dep_value") |> 
  pivot_longer(cols = mstn_to_kreatinin:mstn_to_bcm,
               names_to = "var_indep_name",
               values_to = "var_indep_value") |> 
  mutate(var_dep_value = as.numeric(var_dep_value),
         var_indep_value = as.numeric(var_indep_value)) |> 
  group_by(var_dep_name, var_indep_name) |> 
  nest()




# data_set_filtered <- import("data/processed/data_set_filtered_for_R.csv") # nebylo možné analyzovat data přímo s mix_model

# **********************************************************************
## rCCA, DIABLO ----
# **********************************************************************
d12_join <- d10_aktivita |>
  transmute(
    immet_id   = as.character(sample),
    exam_order = as.character(exam),
    across(where(is.numeric), identity)
  ) |>
  inner_join(
    d11_clinics_imp |>
      transmute(
        immet_id   = as.character(immet_id),
        exam_order = as.character(exam_order),
        across(-c(immet_id, exam_order), identity)
      ),
    by = c("immet_id", "exam_order")
  ) |>
  distinct() |>
  na.omit() |> 
  group_by(immet_id) |>
  filter(n() == 2) |>
  ungroup()

d13_join <- d10_cas_b |>
  transmute(
    immet_id   = as.character(sample),
    exam_order = as.character(exam),
    across(where(is.numeric), identity)
  ) |>
  inner_join(
    d11_clinics_imp |>
      transmute(
        immet_id   = as.character(immet_id),
        exam_order = as.character(exam_order),
        across(-c(immet_id, exam_order), identity)
      ),
    by = c("immet_id", "exam_order")
  ) |>
  distinct() |>
  na.omit() |> 
  group_by(immet_id) |>
  filter(n() == 2) |>
  ungroup()

# **********************************************************************
## MOFA ----

d12_mofa <- d12_join |>
  mutate(
    sample  = paste0(immet_id, "_", exam_order),
    subtype = as.character(disease_subtype) # convenience col for color_by
  ) |>
  distinct(sample, .keep_all = TRUE)

d13_mofa <- d13_join |>
  mutate(
    sample  = paste0(immet_id, "_", exam_order),
    subtype = as.character(disease_subtype) # convenience col for color_by
  ) |>
  distinct(sample, .keep_all = TRUE)

# **********************************************************************
## Ceiling-Effect Replacement Analysis ----
# **********************************************************************
folder_path_2 <- "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Geneactiv/20260302_GGIR_klinika/"

d18_ggir_legend <-import(paste0(folder_path_2, "260302_data_geneactive_clinics.xlsx"),
                         which = "legend") 
d18_ggir_data <- import(paste0(folder_path_2, "260302_data_geneactive_clinics.xlsx"),
                        which = "sheet1") |>
  rename(any_of(setNames(d18_ggir_legend$names_orig,
                         d18_ggir_legend$names_new))) |>
  mutate(
    exam = factor(
      exam,
      levels = unique(exam)[
        order(as.numeric(sub("^M", "", unique(exam))))
      ]
    )
  ) |>
  droplevels()

d18_mvpa_legend <-import(paste0(folder_path_2, "260302_data_geneactive_clinics.xlsx"),
                         which = "legenda_prahy")  
d18_mvpa_data <-import(paste0(folder_path_2, "260302_data_geneactive_clinics.xlsx"),
                       which = "mvpa_prahy") |> 
  rename(any_of(setNames(d18_mvpa_legend$names_orig, 
                         d18_mvpa_legend$names_new))) |> 
  mutate(
    exam = factor(
      exam,
      levels = unique(exam)[
        order(as.numeric(sub("^M", "", unique(exam))))
      ]
    )
  ) |>
  droplevels()

# **********************************************************************
## TIS score ----
# **********************************************************************
d18_ggir_tis_legend <-import("R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/GENEActiv/GGIR/merged_table_final_klinická.xlsx",
                             which = "legenda") 

d18_ggir_tis_data <- import(
  "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/GENEActiv/GGIR/merged_table_final_klinická.xlsx",
  which = 2
) |>
  rename(any_of(setNames(d18_ggir_tis_legend$names_orig,
                         d18_ggir_tis_legend$names_new))) |>
  mutate(
    exam = factor(
      exam,
      levels = unique(exam)[
        order(as.numeric(sub("^M", "", unique(exam))))
      ]
    )
  ) |>
  droplevels()

# **********************************************************************
## PCA, heatmap ----
# *******************************************************
# m0
# data_m0 <- d03 |> 
#   filter(poradie_vysetrenia == "M0") |> 
#   select(where(is.numeric)) |> 
#   select(where(~mean(is.na(.x)) < 0.8)) |> 
#   mice::futuremice(m = 40, maxit = 5, method = "pmm", 
#                    n.core = max(1L, future::availableCores() - 1L),
#                    parallelseed = 123,            # seed pro paralelní běh (doporučeno)
#                    future.plan = "sequential",
#                    ridge = 1e-02,
#                    remove.collinear = TRUE,
#                    remove.constant = TRUE
#   ) |> 
#   mice::complete() |> 
#   select(where(~ !anyNA(.x)))
# group_m0 <- d03 |> 
#   filter(poradie_vysetrenia == "M0") |> 
#   pull(podtyp_nemoci_zjednoduseny)
# # m3
# data_m3 <- d03 |> 
#   filter(poradie_vysetrenia == "M3") |> 
#   select(where(is.numeric)) |> 
#   select(where(~mean(is.na(.x)) < 0.8)) |> 
#   mice::futuremice(m = 40, maxit = 5, method = "pmm", 
#                    n.core = max(1L, future::availableCores() - 1L),
#                    parallelseed = 123,            # seed pro paralelní běh (doporučeno)
#                    future.plan = "sequential",
#                    ridge = 1e-02,
#                    remove.collinear = TRUE,
#                    remove.constant = TRUE
#   ) |> 
#   mice::complete() |> 
#   select(where(~ !anyNA(.x)))
# group_m3 <- d03 |> 
#   filter(poradie_vysetrenia == "M3") |> 
#   pull(podtyp_nemoci_zjednoduseny)
# # m6
# data_m6 <- d03 |> 
#   filter(poradie_vysetrenia == "M6") |> 
#   select(where(is.numeric)) |> 
#   select(where(~mean(is.na(.x)) < 0.8)) |> 
#   mice::futuremice(m = 40, maxit = 5, method = "pmm", 
#                    n.core = max(1L, future::availableCores() - 1L),
#                    parallelseed = 123,            # seed pro paralelní běh (doporučeno)
#                    future.plan = "sequential",
#                    ridge = 1e-02,
#                    remove.collinear = TRUE,
#                    remove.constant = TRUE
#   ) |> 
#   mice::complete() |> 
#   select(where(~ !anyNA(.x)))
# group_m6 <- d03 |> 
#   filter(poradie_vysetrenia == "M6") |> 
#   pull(podtyp_nemoci_zjednoduseny)
# # m18
# data_m18 <- d03 |> 
#   filter(poradie_vysetrenia == "M18") |> 
#   select(where(is.numeric)) |> 
#   select(where(~mean(is.na(.x)) < 0.8)) |> 
#   mice::futuremice(m = 40, maxit = 5, method = "pmm", 
#                    n.core = max(1L, future::availableCores() - 1L),
#                    parallelseed = 123,            # seed pro paralelní běh (doporučeno)
#                    future.plan = "sequential",
#                    ridge = 1e-02,
#                    remove.collinear = TRUE,
#                    remove.constant = TRUE
#   ) |> 
#   mice::complete() |> 
#   select(where(~ !anyNA(.x)))
# group_m18 <- d03 |> 
#   filter(poradie_vysetrenia == "M18") |> 
#   pull(podtyp_nemoci_zjednoduseny)








  



