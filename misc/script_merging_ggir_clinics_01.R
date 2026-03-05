test_tab <- import("R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/GENEActiv/GGIR/merged_table_final_pro_zaverecnou_zpravu_AZV_100226.xlsx",
                   which = "Sheet1") |> 
  dplyr::mutate(immet_id = str_extract(filename...1, "^[^_]+") |> as.character(),
                immet_id = str_replace(immet_id, "57IMMET","IMMET57")) |> 
  rename(`Poradie vyšetrenia` = `...2`) |> 
  filter(!is.na(filename...1))
test_tab2 <- import("R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Data/IMMET_terapie a klinická data_01102025.xlsx",
                    which = "IMMET_kopie_2026") |> 
  mutate(immet_id = str_remove_all(`projekt ID`, "_") |> as.character())

test_tab_fin <- test_tab |> 
  full_join(test_tab2,
            by = c("immet_id", "Poradie vyšetrenia"))

export(test_tab_fin,
       "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/GENEActiv/GGIR/merged_table_final_klinická.xlsx")
