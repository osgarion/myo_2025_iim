# **********************************************************************
# Ceiling/floor: minimal objects (GGIR + clinical merged table)
# **********************************************************************

# Load only the object needed by Script_myo_2025_iim_ceilingEffect_01.R
# Wrappers can override via option("myo_ceil_flo_data_path_01" = "...xlsx")
if (!exists("d18_ggir_data", inherits = FALSE)) {
  data_path_ceil_01 <- getOption(
    "myo_ceil_flo_data_path_01",
    "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Geneactiv/20260302_GGIR_klinika/260302_data_geneactive_clinics.xlsx"
  )

  if (!file.exists(data_path_ceil_01)) {
    stop(
      "Input file for ceiling/floor analysis was not found: ",
      data_path_ceil_01,
      "\nSet option 'myo_ceil_flo_data_path_01' to a valid .xlsx path."
    )
  }

  d18_ggir_legend <- rio::import(data_path_ceil_01, which = "legend")
  d18_ggir_data <- rio::import(data_path_ceil_01, which = "sheet1") |>
    dplyr::rename(any_of(setNames(d18_ggir_legend$names_orig, d18_ggir_legend$names_new))) |>
    dplyr::mutate(
      exam = factor(
        exam,
        levels = unique(exam)[order(as.numeric(sub("^M", "", unique(exam))))]
      )
    ) |>
    droplevels()
}

