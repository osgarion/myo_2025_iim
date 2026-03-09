# **********************************************************************
# Ceiling/floor analysis wrapper: mvpa_t100_time_min
# Run this script as RStudio Background Job.
# **********************************************************************
local({
  options(
    myo_pacman_update = FALSE,
    myo_pacman_install = FALSE,
    myo_ncores = 4L,
    myo_run_renv_status = FALSE
  )
  run_repro_backup_block_ceil_01 <- FALSE
  sel_indep_var_ceil_01 <- c("mvpa_t100_time_min")
  output_root_ceil_01 <- file.path("output", "ceil_flo_mvpa_t100_time_min")
  source("scripts/Script_myo_2025_iim_ceilingEffect_01.R", local = TRUE)
})
