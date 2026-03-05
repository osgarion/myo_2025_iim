# **********************************************************************
# 1. společná příprava dat (nové bout proměnné + jen kompletní M0/M6) ----
# **********************************************************************
# nově upravná tabulka čas
d14_join <- d10_aktivita |>
  transmute(
    immet_id   = as.character(sample),
    exam_order = tolower(as.character(exam)),
    across(where(is.numeric), identity)
  ) |>
  inner_join(
    d10_cas_b |>
      transmute(
        immet_id   = as.character(sample),
        exam_order = tolower(as.character(exam)),
        across(where(is.numeric), identity)
      ),
    by = c("immet_id", "exam_order"),
    suffix = c("_activ", "_time")
  ) |>
  inner_join(
    d11_clinics_imp |>
      mutate(exam_order = tolower(as.character(exam_order))) |>
      mutate(immet_id = as.character(immet_id)),
    by = c("immet_id", "exam_order")
  ) |>
  filter(exam_order %in% c("m0", "m6")) |>
  distinct()

d14_bouts <- d14_join |>
  dplyr::mutate(
    t100_t400_b1_time  = t100_b1m80_time_min  - t400_b1m80_time_min,
    t40_t100_b1_time   = t40_b1m80_time_min   - t100_b1m80_time_min - t400_b1m80_time_min,
    t20_t40_b1_time    = t20_b1m80_time_min   - t40_b1m80_time_min  - t100_b1m80_time_min - t400_b1m80_time_min,
    
    t100_t400_b5_time  = t100_b5m80_time_min  - t400_b5m80_time_min,
    t40_t100_b5_time   = t40_b5m80_time_min   - t100_b5m80_time_min - t400_b5m80_time_min,
    t20_t40_b5_time    = t20_b5m80_time_min   - t40_b5m80_time_min  - t100_b5m80_time_min - t400_b5m80_time_min,
    
    t100_t400_b10_time = t100_b10m80_time_min - t400_b10m80_time_min,
    t40_t100_b10_time  = t40_b10m80_time_min  - t100_b10m80_time_min - t400_b10m80_time_min,
    t20_t40_b10_time   = t20_b10m80_time_min  - t40_b10m80_time_min  - t100_b10m80_time_min - t400_b10m80_time_min,
    
    t100_t400_time = t100_time_min - t400_time_min,
    t40_t100_time = t40_time_min - t100_time_min - t400_time_min,
    t20_t40_time = t20_time_min - t40_time_min - t100_time_min - t400_time_min
  ) |>
  dplyr::rename(
    t400_t8000_b1_time  = t400_b1m80_time_min,
    t400_t8000_b5_time  = t400_b5m80_time_min,
    t400_t8000_b10_time = t400_b10m80_time_min,
    t400_t8000_time = t400_time_min
  ) |>
  dplyr::group_by(immet_id) |>
  dplyr::filter(dplyr::n() == 2) |>
  dplyr::ungroup()

# bout proměnné, které chceme sumarizovat
bout_vars <- c(
  "t20_t40_b1_time", "t40_t100_b1_time", "t100_t400_b1_time", "t400_t8000_b1_time",
  "t20_t40_b5_time", "t40_t100_b5_time", "t100_t400_b5_time", "t400_t8000_b5_time",
  "t20_t40_b10_time","t40_t100_b10_time","t100_t400_b10_time","t400_t8000_b10_time"
)

# **********************************************************************
# 2. Tabulka: disease_subtype x exam_order (mean + sd) ----
# **********************************************************************
# tab_bout_subtype <- d14_bouts |> # previous version
tab_bout_subtype <- d14_join |>
  dplyr::group_by(disease_subtype, exam_order) |>
  dplyr::summarise(
    dplyr::across(
      dplyr::ends_with("time_min"), 
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd   = ~ sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# **********************************************************************
# 3. Tabulka: response_m0_m6 x exam_order (mean + sd) ----
# **********************************************************************
# tab_bout_response <- d14_join |> # previous version
tab_bout_response <- d14_join |>
  dplyr::group_by(response_m0_m6, exam_order) |>
  dplyr::summarise(
    dplyr::across(
      dplyr::ends_with("time_min"), 
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd   = ~ sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# **********************************************************************
# 4. Tabulka: pouze exam_order (mean + sd) ----
# **********************************************************************
# tab_bout_patient <- d14_bouts |># previous version
tab_bout_patient <- d14_join |>
  dplyr::group_by(exam_order) |>
  dplyr::summarise(
    dplyr::across(
      dplyr::ends_with("time_min"), 
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd   = ~ sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

export(list(disease_subtype = tab_bout_subtype, 
            response_m0_m6 = tab_bout_response,
            patient = tab_bout_patient),
       "output/tables/260129_bout_summary_01.xlsx")
