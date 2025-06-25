# data ----
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


# merging
d03 <- d02_clinics |> 
  full_join(d01_myokiny |> 
              select(-maju_klinicke_data_pre_odfiltrovanie),
            by = c("projekt_id",
                   "poradie_vysetrenia",
                   "bbm_id",
                   "pohlavi")) |> 
  ungroup()
  