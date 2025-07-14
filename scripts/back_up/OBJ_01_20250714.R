# variables ----
var_indep_01 <- c("mstn", "fst", "fstl3", "akta")
var_dep_01 <- c("mmt8_total", "mmt10_total", "fi_2", "borg10", "haq", 
                "sf36_mcs", "crp", "ast", "alt", "ck", "ld", 
                "kreatinin_umol_l", "mitax", "myoact", "physician_activity_vas", 
                "muscle_disease_activity", "odpoved_na_terapii_m0_vs_m6")
var_dep_02 <- c("mmt8_total", "mmt10_total", "fi_2")


# data ----
## importing ----
folder_path <- "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Data/"

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



path_clin_01 <- file_to_load(phrase = "IMMET_terapie", 
                            folder_path = folder_path)


d02_clinics <- import(path_clin_01) |>
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
    odpoved_na_terapii_m0_vs_m6 = factor(odpoved_na_terapii_m0_vs_m6),
    odpoved_na_terapii_m0_vs_m18 = factor((odpoved_na_terapii_m0_vs_m18)),
    podtyp_nemoci = factor(podtyp_nemoci),
    kurak_nejmene_100_cigaret = factor(kurak_nejmene_100_cigaret)
  ) |> 
  mutate(across(where(is_01_col), as.factor)) |> 
  mutate(across((vek:last_col()) & where(is.character), as.numeric)) |> 
  group_by(bbm_id) |> 
  fill(kurak_nejmene_100_cigaret)


## merging ----
d03 <- d02_clinics |> 
  full_join(d01_myokiny |> 
              select(-maju_klinicke_data_pre_odfiltrovanie),
            by = c("projekt_id",
                   "poradie_vysetrenia",
                   "bbm_id",
                   "pohlavi")) |> 
  group_by(projekt_id) |> 
  fill(podtyp_nemoci) |> 
  fill(odpoved_na_terapii_m0_vs_m6, .direction = "downup") |> 
  ungroup()

# selection ----
d04_sel1 <- d03 |>
  mutate(
    odpoved_na_terapii_m0_vs_m6 =
      str_replace(odpoved_na_terapii_m0_vs_m6, "No Improvement", "No improvement")
  ) |>
  select(projekt_id, vek, gk, davka_gk_mg_den, podtyp_nemoci, poradie_vysetrenia, all_of(var_dep_01), all_of(var_indep_01)) |>
  filter(poradie_vysetrenia %in% c("M0", "M3", "M6", "M18"))

d04_sel1_nested <- d04_sel1 |> 
  # mutate(poradie_vysetrenia = str_remove(poradie_vysetrenia, "M") |> as.numeric()) |>
  dplyr::select(-odpoved_na_terapii_m0_vs_m6) |> 
  pivot_longer(cols = any_of(var_dep_01),
               names_to = "var_dep_name",
               values_to = "var_dep_value") |> 
  pivot_longer(cols = any_of(var_indep_01),
               names_to = "var_indep_name",
               values_to = "var_indep_value") |> 
  group_by(var_dep_name, var_indep_name) |> 
  nest()

# data_set_filtered <- import("data/processed/data_set_filtered_for_R.csv") # nebylo možné analyzovat data přímo s mix_model

## PCA, heatmap ----
# m0
data_m0 <- d03 |> 
  filter(poradie_vysetrenia == "M0") |> 
  select(where(is.numeric)) |> 
  select(where(~mean(is.na(.x)) < 0.8)) |> 
  mice::mice(method = "pmm", m = 5, printFlag = FALSE) |> 
  mice::complete() |> 
  select(where(~ !anyNA(.x)))
group_m0 <- d03 |> 
  filter(poradie_vysetrenia == "M0") |> 
  pull(podtyp_nemoci)
# m3
data_m3 <- d03 |> 
  filter(poradie_vysetrenia == "M3") |> 
  select(where(is.numeric)) |> 
  select(where(~mean(is.na(.x)) < 0.8)) |> 
  mice::mice(method = "pmm", m = 5, printFlag = FALSE) |> 
  mice::complete() |> 
  select(where(~ !anyNA(.x)))
group_m3 <- d03 |> 
  filter(poradie_vysetrenia == "M3") |> 
  pull(podtyp_nemoci)
# m6
data_m6 <- d03 |> 
  filter(poradie_vysetrenia == "M6") |> 
  select(where(is.numeric)) |> 
  select(where(~mean(is.na(.x)) < 0.8)) |> 
  mice::mice(method = "pmm", m = 5, printFlag = FALSE) |> 
  mice::complete() |> 
  select(where(~ !anyNA(.x)))
group_m6 <- d03 |> 
  filter(poradie_vysetrenia == "M6") |> 
  pull(podtyp_nemoci)
# m18
data_m18 <- d03 |> 
  filter(poradie_vysetrenia == "M18") |> 
  select(where(is.numeric)) |> 
  select(where(~mean(is.na(.x)) < 0.8)) |> 
  mice::mice(method = "pmm", m = 5, printFlag = FALSE) |> 
  mice::complete() |> 
  select(where(~ !anyNA(.x)))
group_m18 <- d03 |> 
  filter(poradie_vysetrenia == "M18") |> 
  pull(podtyp_nemoci)








  