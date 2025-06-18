# data ----
d01_myokiny <- import("R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Data 052025/Myokiny sérum_pro stat.xlsx") |>
  relocate(`věk skupiny`, .before = `věk k M0`) |> 
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
  janitor::clean_names() 

d02_clinics <- import("R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Data 052025/IMMET_terapie a klinická data_15052025_pre stat.xlsx") |>
  select(-1) |> 
  mutate(across(where(is.character), ~na_if(.x, "x"))) |> 
  clean_names()

d03 <- d02_clinics |> 
  full_join(d01_myokiny,
            by = c("immet_id" = "projekt_id",
                   "poradie_vysetrenia"))
  