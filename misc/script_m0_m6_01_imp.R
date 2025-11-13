if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, purrr,conflicted)

## libraries ----
pacman::p_load(update = F,  
               equatiomatic, beepr, tictoc, 
               tidyverse, purrr, furrr, easystats, rio, janitor, ggthemes, car,
               gtsummary, skimr, sjPlot, flextable, ggpubr, rstatix, tidymodels,
               kableExtra, skimr, GGally, testthat, factoextra, gplots, uwot,
               lmerTest, dlookr, multilevelmod,furrr,ggforce, lazyWeave, paletteer,
               emmeans, openxlsx, GLMMadaptive, multidplyr
)

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
  dplyr::recode
)

# directories ----
setwd("F:/Analysis/Vernerová Lucia/myo_2025_iim")

# path_export <- "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Output/251002/m0_m6_imp/"
path_export <- "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Output/251017/m0_m6_imp/"


dir.create(path_export, recursive = TRUE)

# covariate models 2 ----
## only covariates ----
### models ----

cluster <- multidplyr::new_cluster(ncores)
multidplyr::cluster_library(cluster, c("dplyr","purrr","stringr","tibble", 
                                       "data.table", "tidyverse", "tidymodels",
                                       "multilevelmod","lme4","lmerTest",
                                       "broom.mixed", "easystats"))
multidplyr::cluster_copy(cluster, c("model_covar_2","model_min_2", "model_comp_int"))

res_mixMod_covar_02 <- d04_sel2 |>
  select(any_of(var_dep_01 |>  str_subset("odpoved_na_terapii_m0_vs_m6|kreatinin_umol_l", negate = TRUE)),
         any_of(var_indep_01),
         vek, denna_davka_gk_mg_kg_bw, kreatinin_umol_l, bcm, pohlavi, 
         podtyp_nemoci_zjednoduseny, jo_1, anti_hmgcr, projekt_id, poradie_vysetrenia) |> 
  mutate(across(.cols = c("ast", "alt", "ck", "crp", "haq", "ld", "mitax", "myoact", "myoglobin"), function(x) log(x+1)),
         jo_1 = as.factor(jo_1),
         anti_hmgcr = as.factor(anti_hmgcr)) |> 
  filter(poradie_vysetrenia %in% c("M0", "M6")) |> 
  mice::mice("pmm", m = 20,
       maxit = 15,
       printFlag = F) |> mice::complete() |>
  pivot_longer(cols = any_of(var_dep_01 |>  str_subset("odpoved_na_terapii_m0_vs_m6|kreatinin_umol_l", negate = TRUE)),
               names_to = "var_dep_name",
               values_to = "var_dep_value") |> 
  pivot_longer(cols = any_of(var_indep_01),
               names_to = "var_indep_name",
               values_to = "var_indep_value") |> 
  group_by(var_dep_name, var_indep_name) |> 
  multidplyr::partition(cluster) |> 
  do(
    data = data.table(.),
    mod_without = model_min_2(., filter_na = FALSE),
    mod_covar_gk = model_covar_2(., "denna_davka_gk_mg_kg_bw", filter_na = FALSE),
    mod_covar_age = model_covar_2(., "vek", filter_na = FALSE),
    mod_covar_bcm = model_covar_2(., "bcm", filter_na = FALSE),
    mod_covar_creatine = model_covar_2(., "kreatinin_umol_l", filter_na = FALSE)
  ) |> 
  collect()

walk(cluster, ~ .x$kill()); gc()


### tables ----
# without confounders
res_mixMod_covar_02_tab <- res_mixMod_covar_02 |> 
  select(var_dep_name, var_indep_name, mod_without) |> 
  mutate(mod_com_without = map(mod_without, ~.x$table)) |> 
  unnest(mod_com_without) |> 
  select(-mod_without) |> 
  arrange(p_val_covar) |> 
  rename_with(~ paste0(.x, "_without"),
              .cols = c("omega_sq_p",
                        "r2_mod",
                        "p_val_covar"))
# glucocorticoids dose per day and weight
res_mixMod_covar_02_tab_gk <- res_mixMod_covar_02 |> 
  select(var_dep_name, var_indep_name, mod_covar_gk) |> 
  mutate(mod_com_gk = map(mod_covar_gk, ~.x$table)) |> 
  unnest(mod_com_gk) |> 
  select(-mod_covar_gk) |> 
  arrange(p_val_covar) |> 
  rename_with(~ paste0(.x, "_gk"),
              .cols = c("omega_sq_p",
                        "r2_mod",
                        "p_val_covar"))
# age
res_mixMod_covar_02_tab_age <- res_mixMod_covar_02 |> 
  select(var_dep_name, var_indep_name, mod_covar_age) |> 
  mutate(mod_com_age = map(mod_covar_age, ~.x$table)) |> 
  unnest(mod_com_age) |> 
  select(-mod_covar_age) |> 
  arrange(p_val_covar) |> 
  rename_with(~ paste0(.x, "_age"),
              .cols = c("omega_sq_p",
                        "r2_mod",
                        "p_val_covar"))
# creatine
res_mixMod_covar_02_tab_creatine <- res_mixMod_covar_02 |> 
  select(var_dep_name, var_indep_name, mod_covar_creatine) |> 
  mutate(mod_com_creatine = map(mod_covar_creatine, ~.x$table)) |> 
  unnest(mod_com_creatine) |> 
  select(-mod_covar_creatine) |> 
  arrange(p_val_covar) |> 
  rename_with(~ paste0(.x, "_creatine"),
              .cols = c("omega_sq_p",
                        "r2_mod",
                        "p_val_covar"))
# body cell mass
res_mixMod_covar_02_tab_bcm <- res_mixMod_covar_02 |> 
  select(var_dep_name, var_indep_name, mod_covar_bcm) |> 
  mutate(mod_com_bcm = map(mod_covar_bcm, ~.x$table)) |> 
  unnest(mod_com_bcm) |> 
  select(-mod_covar_bcm) |> 
  arrange(p_val_covar) |> 
  rename_with(~ paste0(.x, "_bcm"),
              .cols = c("omega_sq_p",
                        "r2_mod",
                        "p_val_covar"))
# funal table
res_mixMod_covar_02_tab_confounders <- res_mixMod_covar_02_tab |> 
  full_join(res_mixMod_covar_02_tab_age, by = c("var_dep_name", "var_indep_name")) |> 
  full_join(res_mixMod_covar_02_tab_gk, by = c("var_dep_name", "var_indep_name")) |> 
  full_join(res_mixMod_covar_02_tab_bcm, by = c("var_dep_name", "var_indep_name")) |> 
  full_join(res_mixMod_covar_02_tab_creatine, by = c("var_dep_name", "var_indep_name")) |> 
  arrange(var_dep_name)

## subgrouping ----
### age ----
res_mixMod_covar_02_sub_age <- res_mixMod_covar_02 |> 
  select(var_dep_name, var_indep_name, data) |> 
  inner_join(res_mixMod_covar_02_tab_age |> 
               filter(if_any(last_col(), ~ . < 0.06)) |> 
               select(var_dep_name, var_indep_name),
             by = c("var_dep_name", "var_indep_name")) |> 
  mutate(fig_sex = pmap(list(data = data, 
                             ylab = var_dep_name, 
                             xlab = var_indep_name,
                             covar = "vek",
                             covar_spec = "pohlavi"),
                        plot_lmer_grid_basic_02),
         fig_sub = pmap(list(data = data, 
                             ylab = var_dep_name, 
                             xlab = var_indep_name,
                             covar = "vek",
                             covar_spec = "podtyp_nemoci_zjednoduseny"),
                        plot_lmer_grid_basic_02),
         fig_jo = pmap(list(data = data, 
                            ylab = var_dep_name, 
                            xlab = var_indep_name,
                            covar = "vek",
                            covar_spec = "jo_1"),
                       plot_lmer_grid_basic_02),
         fig_hmgcr = pmap(list(data = data, 
                               ylab = var_dep_name, 
                               xlab = var_indep_name,
                               covar = "vek",
                               covar_spec = "anti_hmgcr"),
                          plot_lmer_grid_basic_02)
  )

### gk ----
res_mixMod_covar_02_sub_gk <- res_mixMod_covar_02 |> 
  select(var_dep_name, var_indep_name, data) |> 
  inner_join(res_mixMod_covar_02_tab_gk |> 
               filter(if_any(last_col(), ~ . < 0.06)) |> 
               select(var_dep_name, var_indep_name),
             by = c("var_dep_name", "var_indep_name")) |> 
  mutate(fig_sex = pmap(list(data = data, 
                             ylab = var_dep_name, 
                             xlab = var_indep_name,
                             covar = "denna_davka_gk_mg_kg_bw",
                             covar_spec = "pohlavi"),
                        plot_lmer_grid_basic_02),
         fig_sub = pmap(list(data = data, 
                             ylab = var_dep_name, 
                             xlab = var_indep_name,
                             covar = "denna_davka_gk_mg_kg_bw",
                             covar_spec = "podtyp_nemoci_zjednoduseny"),
                        plot_lmer_grid_basic_02),
         fig_jo = pmap(list(data = data, 
                            ylab = var_dep_name, 
                            xlab = var_indep_name,
                            covar = "denna_davka_gk_mg_kg_bw",
                            covar_spec = "jo_1"),
                       plot_lmer_grid_basic_02),
         fig_hmgcr = pmap(list(data = data, 
                               ylab = var_dep_name, 
                               xlab = var_indep_name,
                               covar = "denna_davka_gk_mg_kg_bw",
                               covar_spec = "anti_hmgcr"),
                          plot_lmer_grid_basic_02)
  )

### bcm ----
res_mixMod_covar_02_sub_bcm <- res_mixMod_covar_02 |> 
  select(var_dep_name, var_indep_name, data) |> 
  inner_join(res_mixMod_covar_02_tab_bcm |> 
               filter(if_any(last_col(), ~ . < 0.06)) |> 
               select(var_dep_name, var_indep_name),
             by = c("var_dep_name", "var_indep_name")) |> 
  mutate(fig_sex = pmap(list(data = data, 
                             ylab = var_dep_name, 
                             xlab = var_indep_name,
                             covar = "bcm",
                             covar_spec = "pohlavi"),
                        plot_lmer_grid_basic_02),
         fig_sub = pmap(list(data = data, 
                             ylab = var_dep_name, 
                             xlab = var_indep_name,
                             covar = "bcm",
                             covar_spec = "podtyp_nemoci_zjednoduseny"),
                        plot_lmer_grid_basic_02),
         fig_jo = pmap(list(data = data, 
                            ylab = var_dep_name, 
                            xlab = var_indep_name,
                            covar = "bcm",
                            covar_spec = "jo_1"),
                       plot_lmer_grid_basic_02),
         fig_hmgcr = pmap(list(data = data, 
                               ylab = var_dep_name, 
                               xlab = var_indep_name,
                               covar = "bcm",
                               covar_spec = "anti_hmgcr"),
                          plot_lmer_grid_basic_02)
  )

### creatine ----
res_mixMod_covar_02_sub_creatine <- res_mixMod_covar_02 |> 
  select(var_dep_name, var_indep_name, data) |> 
  inner_join(res_mixMod_covar_02_tab_creatine |> 
               filter(if_any(last_col(), ~ . < 0.06)) |> 
               select(var_dep_name, var_indep_name),
             by = c("var_dep_name", "var_indep_name")) |> 
  mutate(fig_sex = pmap(list(data = data, 
                             ylab = var_dep_name, 
                             xlab = var_indep_name,
                             covar = "kreatinin_umol_l",
                             covar_spec = "pohlavi"),
                        plot_lmer_grid_basic_02),
         fig_sub = pmap(list(data = data, 
                             ylab = var_dep_name, 
                             xlab = var_indep_name,
                             covar = "kreatinin_umol_l",
                             covar_spec = "podtyp_nemoci_zjednoduseny"),
                        plot_lmer_grid_basic_02),
         fig_jo = pmap(list(data = data, 
                            ylab = var_dep_name, 
                            xlab = var_indep_name,
                            covar = "kreatinin_umol_l",
                            covar_spec = "jo_1"),
                       plot_lmer_grid_basic_02),
         fig_hmgcr = pmap(list(data = data, 
                               ylab = var_dep_name, 
                               xlab = var_indep_name,
                               covar = "kreatinin_umol_l",
                               covar_spec = "anti_hmgcr"),
                          plot_lmer_grid_basic_02)
  )

### gk + vek  ----
res_mixMod_covar_02_sub_gk_age <- res_mixMod_covar_02 |> 
  select(var_dep_name, var_indep_name, data) |> 
  inner_join(res_mixMod_covar_02_tab_gk|> 
               filter(if_any(last_col(), ~ . < 0.06)) |> 
               select(var_dep_name, var_indep_name),
             by = c("var_dep_name", "var_indep_name")) |> 
  mutate(fig_sex = pmap(list(data = data, 
                             ylab = var_dep_name, 
                             xlab = var_indep_name,
                             covar = "denna_davka_gk_mg_kg_bw + vek",
                             covar_spec = "pohlavi"),
                        plot_lmer_grid_basic_02),
         fig_sub = pmap(list(data = data, 
                             ylab = var_dep_name, 
                             xlab = var_indep_name,
                             covar = "denna_davka_gk_mg_kg_bw + vek",
                             covar_spec = "podtyp_nemoci_zjednoduseny"),
                        plot_lmer_grid_basic_02),
         fig_jo = pmap(list(data = data, 
                            ylab = var_dep_name, 
                            xlab = var_indep_name,
                            covar = "denna_davka_gk_mg_kg_bw + vek",
                            covar_spec = "jo_1"),
                       plot_lmer_grid_basic_02),
         fig_hmgcr = pmap(list(data = data, 
                               ylab = var_dep_name, 
                               xlab = var_indep_name,
                               covar = "denna_davka_gk_mg_kg_bw + vek",
                               covar_spec = "anti_hmgcr"),
                          plot_lmer_grid_basic_02)
  )



# export ----
## table ----
export(res_mixMod_covar_02_tab_confounders, 
       paste0(path_export,"covariates_01.xlsx"))

## age ----
path_export_age <- paste0(path_export,"age/")
dir.create(path_export_age, recursive = TRUE)

### sex ----
pdf(paste0(path_export_age, "sex.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_age$fig_sex,print)
dev.off()

### subdiag ----
pdf(paste0(path_export_age, "subdiag.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_age$fig_sub,print)
dev.off()

### jo_1 ----
pdf(paste0(path_export_age, "jo_1.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_age$fig_jo,print)
dev.off()

### anti_hgmcr ----
pdf(paste0(path_export_age, "anti_hgmcr.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_age$fig_hmgcr,print)
dev.off()

## bcm ----
path_export_bcm <- paste0(path_export,"bcm/")
dir.create(path_export_bcm, recursive = TRUE)

### sex ----
pdf(paste0(path_export_bcm, "sex.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_bcm$fig_sex,print)
dev.off()

### subdiag ----
pdf(paste0(path_export_bcm, "subdiag.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_bcm$fig_sub,print)
dev.off()

### jo_1 ----
pdf(paste0(path_export_bcm, "jo_1.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_bcm$fig_jo,print)
dev.off()

### anti_hgmcr ----
pdf(paste0(path_export_bcm, "anti_hgmcr.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_bcm$fig_hmgcr,print)
dev.off()

## creatine ----
path_export_creatine <- paste0(path_export,"creatine/")
dir.create(path_export_creatine, recursive = TRUE)

### sex ----
pdf(paste0(path_export_creatine, "sex.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_creatine$fig_sex,print)
dev.off()

### subdiag ----
pdf(paste0(path_export_creatine, "subdiag.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_creatine$fig_sub,print)
dev.off()

### jo_1 ----
pdf(paste0(path_export_creatine, "jo_1.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_creatine$fig_jo,print)
dev.off()

### anti_hgmcr ----
pdf(paste0(path_export_creatine, "anti_hgmcr.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_creatine$fig_hmgcr,print)
dev.off()

## gk ----
path_export_gk <- paste0(path_export,"gk/")
dir.create(path_export_gk, recursive = TRUE)

### sex ----
pdf(paste0(path_export_gk, "sex.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_gk$fig_sex,print)
dev.off()

### subdiag ----
pdf(paste0(path_export_gk, "subdiag.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_gk$fig_sub,print)
dev.off()

### jo_1 ----
pdf(paste0(path_export_gk, "jo_1.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_gk$fig_jo,print)
dev.off()

### anti_hgmcr ----
pdf(paste0(path_export_gk, "anti_hgmcr.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_gk$fig_hmgcr,print)
dev.off()

## gk + vek----
path_export_gk_age <- paste0(path_export,"gk_age/")
dir.create(path_export_gk_age, recursive = TRUE)

### sex ----
pdf(paste0(path_export_gk_age, "sex.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_gk_age$fig_sex,print)
dev.off()

### subdiag ----
pdf(paste0(path_export_gk_age, "subdiag.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_gk_age$fig_sub,print)
dev.off()

### jo_1 ----
pdf(paste0(path_export_gk_age, "jo_1.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_gk_age$fig_jo,print)
dev.off()

### anti_hgmcr ----
pdf(paste0(path_export_gk_age, "anti_hgmcr.pdf"), width = 15, height = 8)
walk(res_mixMod_covar_02_sub_gk_age$fig_hmgcr,print)
dev.off()
