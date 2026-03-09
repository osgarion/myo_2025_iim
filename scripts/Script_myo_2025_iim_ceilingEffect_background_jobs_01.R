# **********************************************************************
# Batch launcher for ceiling/floor analyses in RStudio Background Jobs
# **********************************************************************

job_scripts_ceil_01 <- c(
  "scripts/Script_myo_2025_iim_ceilingEffect_p9931_enmo_01.R",
  "scripts/Script_myo_2025_iim_ceilingEffect_m5_enmo_01.R",
  "scripts/Script_myo_2025_iim_ceilingEffect_ig_intercept_enmo_01.R",
  "scripts/Script_myo_2025_iim_ceilingEffect_acc_day_spt_wei_01.R",
  "scripts/Script_myo_2025_iim_ceilingEffect_mvpa_t100_time_min_01.R"
)

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  for (sc_01 in job_scripts_ceil_01) {
    job_name_01 <- gsub("^scripts/|\\.R$", "", sc_01)
    rstudioapi::jobRunScript(
      path = sc_01,
      name = job_name_01,
      workingDir = getwd(),
      importEnv = FALSE
    )
  }
} else {
  # Fallback for non-interactive contexts (e.g., this launcher started as a background job).
  # Launch each wrapper as an independent Rscript process.
  rscript_bin_01 <- file.path(R.home("bin"), "Rscript.exe")
  if (!file.exists(rscript_bin_01)) {
    stop("Rscript executable not found at: ", rscript_bin_01)
  }

  for (sc_01 in job_scripts_ceil_01) {
    sc_abs_01 <- normalizePath(sc_01, winslash = "/", mustWork = TRUE)
    system2(
      command = rscript_bin_01,
      args = c("--vanilla", shQuote(sc_abs_01)),
      wait = FALSE
    )
  }
}

