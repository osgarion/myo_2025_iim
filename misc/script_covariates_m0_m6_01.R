# covariate models 2 ----
## functions ----
model_min_2 <- function(data, filter_na = TRUE) {
  
  if (filter_na){
    data <- data |> 
      select(var_dep_value, var_indep_value, poradie_vysetrenia, projekt_id) |> 
      na.omit()
  }
  
  ## tidymodels
  lmer_mod <- 
    linear_reg() |> 
    set_engine("lmer", REML = FALSE) |> 
    set_mode("regression")
  
  mod1 <- fit(
    lmer_mod,
    formula = var_dep_value ~ poradie_vysetrenia + var_indep_value + (1 | projekt_id),
    data = data
  )
  
  omega <- sjstats::anova_stats(mod1$fit)$partial.omegasq[[2]]
  r2_tbl <- performance::r2_nakagawa(mod1$fit)
  r2_mod <- r2_tbl$R2_conditional
  beta_covar <- model_parameters(mod1$fit)$Coefficient[[3]]
  lower_covar <- model_parameters(mod1$fit)$CI_low[[3]]
  upper_covar <- model_parameters(mod1$fit)$CI_high[[3]]
  p_val_covar <- model_parameters(mod1$fit)$p[[3]]
  std_beta_covar <- tryCatch(
    parameters::standardize_parameters(mod1$fit, method = "refit")$Std_Coefficient[[3]],
    error = function(e) parameters::standardize_parameters(mod1$fit, method = "posthoc")$Std_Coefficient[[3]]
  )

  return(list(table = tibble(omega_sq_p = omega,
                             r2_mod = r2_mod,
                             beta_covar = beta_covar,
                             std_beta_covar = std_beta_covar,
                             lower_covar = lower_covar,
                             upper_covar = upper_covar,
                             p_val_covar = p_val_covar),
              model = mod1)
  )


}

model_covar_2 <- function(data, covar, filter_na = TRUE) {
  
  if (filter_na){
    data <- data |> 
      select(var_dep_value, var_indep_value, poradie_vysetrenia, all_of(covar), projekt_id) |> 
      na.omit()
  }
  
  ## formula specification
  rhs_terms <- c("poradie_vysetrenia", "var_indep_value", covar, "(1 | projekt_id)")
  fml <- as.formula(paste("var_dep_value ~", paste(rhs_terms, collapse = " + ")))
  
  ## tidymodels
  lmer_mod <- 
    linear_reg() |> 
    set_engine("lmer", REML = FALSE) |> 
    set_mode("regression")
  
  mod1 <- fit(
    lmer_mod,
    formula = fml,
    data = data
  )
  
  
  omega <- sjstats::anova_stats(mod1$fit)$partial.omegasq[[2]]
  r2_tbl <- performance::r2_nakagawa(mod1$fit)
  r2_mod <- r2_tbl$R2_conditional
  beta_covar <- model_parameters(mod1$fit)$Coefficient[[3]]
  lower_covar <- model_parameters(mod1$fit)$CI_low[[3]]
  upper_covar <- model_parameters(mod1$fit)$CI_high[[3]]
  p_val_covar <- model_parameters(mod1$fit)$p[[3]]
  std_beta_covar <- tryCatch(
    parameters::standardize_parameters(mod1$fit, method = "refit")$Std_Coefficient[[3]],
    error = function(e) parameters::standardize_parameters(mod1$fit, method = "posthoc")$Std_Coefficient[[3]]
  )

  return(list(table = tibble(omega_sq_p = omega,
                             r2_mod = r2_mod,
                             beta_covar = beta_covar,
                             std_beta_covar = std_beta_covar,
                             lower_covar = lower_covar,
                             upper_covar = upper_covar,
                             p_val_covar = p_val_covar),
              model = mod1)
  )

}

## only covariates ----
### models ----
options(renv.config.auto.snapshot = FALSE)
tic()
Sys.time()

future::plan(multisession, workers = ncores)

res_mixMod_covar_02 <- d04_sel2 |>
  filter(poradie_vysetrenia %in% c("M0", "M6")) |>
  select(any_of(var_dep_01 |>  str_subset("odpoved_na_terapii_m0_vs_m6|kreatinin_umol_l", negate = TRUE)),
         any_of(var_indep_01),
         vek, denna_davka_gk_mg_kg_bw, kreatinin_umol_l, bcm, pohlavi,
         podtyp_nemoci_zjednoduseny, jo_1, anti_hmgcr, projekt_id, poradie_vysetrenia) |>
  mutate(across(.cols = c("ast", "alt", "ck", "crp", "haq", "ld", "mitax", "myoact", "myoglobin"), function(x) log(x+1)),
         jo_1 = as.factor(jo_1),
         anti_hmgcr = as.factor(anti_hmgcr)) |>
  # mice::mice("pmm", m = 20,
  #      maxit = 15,
  #      printFlag = F) |> complete() |>
  pivot_longer(cols = any_of(var_dep_01 |>  str_subset("odpoved_na_terapii_m0_vs_m6|kreatinin_umol_l", negate = TRUE)),
               names_to = "var_dep_name",
               values_to = "var_dep_value") |>
  pivot_longer(cols = any_of(var_indep_01),
               names_to = "var_indep_name",
               values_to = "var_indep_value") |>
  group_by(var_dep_name, var_indep_name) |>
  nest() |>
  mutate(
    mod_without       = furrr::future_map(data, ~model_min_2(.x,                              filter_na = FALSE), .options = furrr_options(seed = TRUE)),
    mod_covar_gk      = furrr::future_map(data, ~model_covar_2(.x, "denna_davka_gk_mg_kg_bw", filter_na = FALSE), .options = furrr_options(seed = TRUE)),
    mod_covar_age     = furrr::future_map(data, ~model_covar_2(.x, "vek",                     filter_na = FALSE), .options = furrr_options(seed = TRUE)),
    mod_covar_bcm     = furrr::future_map(data, ~model_covar_2(.x, "bcm",                     filter_na = FALSE), .options = furrr_options(seed = TRUE)),
    mod_covar_creatine = furrr::future_map(data, ~model_covar_2(.x, "kreatinin_umol_l",        filter_na = FALSE), .options = furrr_options(seed = TRUE))
  ) |>
  ungroup()

future::plan(sequential); gc()

toc()
beep(4)

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
                        "beta_covar",
                        "std_beta_covar",
                        "lower_covar",
                        "upper_covar",
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
                        "beta_covar",
                        "std_beta_covar",
                        "lower_covar",
                        "upper_covar",
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
                        "beta_covar",
                        "std_beta_covar",
                        "lower_covar",
                        "upper_covar",
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
                        "beta_covar",
                        "std_beta_covar",
                        "lower_covar",
                        "upper_covar",
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
                        "beta_covar",
                        "std_beta_covar",
                        "lower_covar",
                        "upper_covar",
                        "p_val_covar"))
# final table
res_mixMod_covar_02_tab_confounders <- res_mixMod_covar_02_tab |> 
  full_join(res_mixMod_covar_02_tab_age, by = c("var_dep_name", "var_indep_name")) |> 
  full_join(res_mixMod_covar_02_tab_gk, by = c("var_dep_name", "var_indep_name")) |> 
  full_join(res_mixMod_covar_02_tab_bcm, by = c("var_dep_name", "var_indep_name")) |> 
  full_join(res_mixMod_covar_02_tab_creatine, by = c("var_dep_name", "var_indep_name"))

### export ----
export(res_mixMod_covar_02_tab_confounders, 
       paste0("output/tables/", format(Sys.Date(), "%y%m%d"),"_covariates_0_m6_misc_01.xlsx"))