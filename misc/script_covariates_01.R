## models interaction ----
tic()
Sys.time()

cluster <- multidplyr::new_cluster(ncores)
multidplyr::cluster_library(cluster, c("dplyr","purrr","stringr","tibble",
                                       "data.table", "tidyverse", "tidymodels",
                                       "multilevelmod","lme4","lmerTest",
                                       "broom.mixed", "easystats"))
multidplyr::cluster_copy(cluster, c("model_comp","model_min", "model_comp_int"))

res_mixMod_covar_01b <- d04_sel2 |>
  dplyr::mutate(across(.cols = c("ast", "alt", "ck", "crp", "haq", "ld", "mitax", "myoact", "myoglobin"), function(x) log(x + 1))) |>
  pivot_longer(
    cols = any_of(var_dep_01 |> str_subset("odpoved_na_terapii_m0_vs_m6|kreatinin_umol_l", negate = TRUE)),
    names_to = "var_dep_name",
    values_to = "var_dep_value"
  ) |>
  pivot_longer(
    cols = any_of(var_indep_01),
    names_to = "var_indep_name",
    values_to = "var_indep_value"
  ) |>
  group_by(var_dep_name, var_indep_name) |>
  # filter(!str_detect(var_dep_name, "_total")) |>
  multidplyr::partition(cluster) |>
  do(
    data = data.table(.),
    mod_int_vek = model_comp_int(data = ., "vek")
  )  |>
  collect()

walk(cluster, ~ .x$kill()); gc()

toc()
beep(4)

tic()
Sys.time()
res_mixMod_covar_01 <- d04_sel2 |> 
  mutate(across(.cols = c("ast", "alt", "ck", "crp", "haq", "ld", "mitax", "myoact", "myoglobin"), function(x) log(x+1))) |> 
  pivot_longer(cols = any_of(var_dep_01 |>  str_subset("odpoved_na_terapii_m0_vs_m6|kreatinin_umol_l", negate = TRUE)),
               names_to = "var_dep_name",
               values_to = "var_dep_value") |> 
  pivot_longer(cols = any_of(var_indep_01),
               names_to = "var_indep_name",
               values_to = "var_indep_value") |> 
  group_by(var_dep_name, var_indep_name) |> 
  nest() |> 
  ungroup() |> 
  # filter(!str_detect(var_dep_name, "_total")) |>           # drop of mmt-8 and mmt-10
  mutate(mod_vek = map(data, ~model_comp(.x,"vek")),
         mod_gk = map(data, ~model_comp(.x,"denna_davka_gk_mg_kg_bw")),
         mod_ck = map(data, ~model_comp(.x,"kreatinin_umol_l")),
         mod_bcm = map(data, ~model_comp(.x,"bcm")),
         mod_sex = map(data, ~model_comp(.x,"pohlavi")),
         mod_sub = map(data, ~model_comp(.x,"podtyp_nemoci_zjednoduseny")),
         mod_jo = map(data, ~model_comp(.x,"jo_1")),
         mod_hmgcr = map(data, ~model_comp(.x,"anti_hmgcr")),
         mod_int_sex = map(data, ~model_comp_int(data = ., "vek")),
         mod_covar_vek = map(mod_vek, ~.x$model),
         mod_covar_gk = map(mod_gk, ~.x$model),
         mod_covar_ck = map(mod_ck, ~.x$model),
         mod_covar_bcm = map(mod_bcm, ~.x$model),
         mod_covar_sex = map(mod_sex, ~.x$model),
         mod_covar_sub = map(mod_sub, ~.x$model),
         mod_covar_jo = map(mod_jo, ~.x$model),
         mod_covar_hmgcr = map(mod_hmgcr, ~.x$model),
         mod_covar_int_sex = map(mod_int_sex, ~.x$model),
         mod_NOcovar_vek = map(mod_vek, ~.x$model_noCovar),
         mod_NOcovar_gk = map(mod_gk, ~.x$model_noCovar),
         mod_NOcovar_ck = map(mod_ck, ~.x$model_noCovar),
         mod_NOcovar_bcm = map(mod_bcm, ~.x$model_noCovar),
         mod_NOcovar_sex = map(mod_sex, ~.x$model_noCovar),
         mod_NOcovar_sub = map(mod_sub, ~.x$model_noCovar),
         mod_NOcovar_jo = map(mod_jo, ~.x$model_noCovar),
         mod_NOcovar_hmgcr = map(mod_hmgcr, ~.x$model_noCovar),
         mod_NOcovar_int_sex = map(mod_int_sex, ~.x$model_noCovar),
         mod_without = map(data, ~model_min(.x))
  ) 

toc()
beep(4)
