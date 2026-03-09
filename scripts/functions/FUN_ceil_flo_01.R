# **********************************************************************
# Ceiling/floor: minimal libraries and conflict setup
# **********************************************************************

if (!require("pacman")) install.packages("pacman")

pacman_update_ceil_01 <- getOption("myo_pacman_update", FALSE)
pacman_install_ceil_01 <- getOption("myo_pacman_install", FALSE)

if ("package:conflicted" %in% search()) {
  detach("package:conflicted", unload = TRUE, character.only = TRUE)
}

pacman::p_load(
  update = pacman_update_ceil_01,
  install = pacman_install_ceil_01,
  tidyverse, purrr, renv,
  rio, janitor, GGally, corrplot, sjPlot,
  ggpubr, lme4, lmerTest, performance,
  broom.mixed, patchwork, future, furrr,
  openxlsx
)

if (requireNamespace("conflicted", quietly = TRUE)) {
  library(conflicted)
  conflicted::conflicts_prefer(
    dplyr::filter,
    dplyr::mutate,
    dplyr::rename,
    dplyr::summarize,
    dplyr::summarise,
    dplyr::select,
    dplyr::lag,
    purrr::map,
    tidyr::extract,
    janitor::clean_names
  )
}

# Conservative default for parallel sections called by sourced code.
if (is.null(getOption("myo_ncores"))) {
  options(myo_ncores = max(1L, min(4L, future::availableCores() - 1L)))
}

# Keep kable output behavior consistent with the rest of project scripts.
options(knitr.kable.NA = "", scipen = 999)
