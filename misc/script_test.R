if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, purrr,conflicted)

## libraries ----
pacman::p_load(update = T,  
               equatiomatic, beepr, tictoc, 
               tidyverse, purrr, furrr, easystats, rio, janitor, ggthemes, car,
               gtsummary, skimr, sjPlot, flextable, ggpubr, tidymodels,
               kableExtra, skimr, GGally, testthat, factoextra, gplots, uwot,
               lmerTest, dlookr, multilevelmod,furrr,ggforce, lazyWeave, paletteer,
               emmeans, openxlsx, GLMMadaptive
)

library(rstatix)

setwd("F:/Analysis/Vernerová Lucia/myo_2025_iim")

# Missing values, multivariate analyses
pacman::p_load(naniar, MVN) 

## function specification ----
conflicted::conflicts_prefer(
  janitor::remove_empty,
  dplyr::filter,
  dplyr::mutate,
  dplyr::rename,
  dplyr::summarize,
  dplyr::summarise,
  dplyr::select,
  purrr::map,
  tidyr::extract,
  janitor::clean_names,
  dplyr::relocate,
  lmerTest::lmer,
  dplyr::recode,
  GLMMadaptive::predict
)

pdf("output/figures/250820_podtyp_period_01.pdf", width = 14)
walk(eda_period_type_01$fig, print)
dev.off()