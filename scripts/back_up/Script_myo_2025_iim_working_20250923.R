# Back-up
back_up("scripts/functions/FUN_01.R") # the the destination subdirectory specify using 'path_dest'
back_up("scripts/functions/OBJ_01.R") # the the destination subdirectory specify using 'path_dest'
back_up("scripts/Script_myo_2025_iim_working.R") # the the destination subdirectory specify using 'path_dest'


# 250919 ----
# covariate models
# tic()
# Sys.time()
# 
# res_mixMod_covar_01 <- d04_sel2 |> 
#   mutate(across(.cols = c("ast", "alt", "ck", "crp", "haq", "ld", "mitax", "myoact", "myoglobin"), function(x) log(x+1))) |> 
#   pivot_longer(cols = any_of(var_dep_01 |>  str_subset("odpoved_na_terapii_m0_vs_m6|kreatinin_umol_l", negate = TRUE)),
#                names_to = "var_dep_name",
#                values_to = "var_dep_value") |> 
#   pivot_longer(cols = any_of(var_indep_01),
#                names_to = "var_indep_name",
#                values_to = "var_indep_value") |> 
#   group_by(var_dep_name, var_indep_name) |> 
#   nest() |> 
#   ungroup() |> 
#   filter(!str_detect(var_dep_name, "_total")) |>
#   mutate(mod_vek = map(data, ~model_comp(.x,"vek")),
#          mod_gk = map(data, ~model_comp(.x,"denna_davka_gk_mg_kg_bw")),
#          mod_ck = map(data, ~model_comp(.x,"kreatinin_umol_l")),
#          mod_bcm = map(data, ~model_comp(.x,"bcm")),
#          mod_sex = map(data, ~model_comp(.x,"pohlavi")),
#          mod_sub = map(data, ~model_comp(.x,"podtyp_nemoci_zjednoduseny")),
#          mod_jo = map(data, ~model_comp(.x,"jo_1")),
#          mod_hmgcr = map(data, ~model_comp(.x,"anti_hmgcr")),
#          mod_covar_vek = map(mod_vek, ~.x$model),
#          mod_covar_gk = map(mod_gk, ~.x$model),
#          mod_covar_ck = map(mod_ck, ~.x$model),
#          mod_covar_bcm = map(mod_bcm, ~.x$model),
#          mod_covar_sex = map(mod_sex, ~.x$model),
#          mod_covar_sub = map(mod_sub, ~.x$model),
#          mod_covar_jo = map(mod_jo, ~.x$model),
#          mod_covar_hmgcr = map(mod_hmgcr, ~.x$model),
#          mod_NOcovar_vek = map(mod_vek, ~.x$model_noCovar),
#          mod_NOcovar_gk = map(mod_gk, ~.x$model_noCovar),
#          mod_NOcovar_ck = map(mod_ck, ~.x$model_noCovar),
#          mod_NOcovar_bcm = map(mod_bcm, ~.x$model_noCovar),
#          mod_NOcovar_sex = map(mod_sex, ~.x$model_noCovar),
#          mod_NOcovar_sub = map(mod_sub, ~.x$model_noCovar),
#          mod_NOcovar_jo = map(mod_jo, ~.x$model_noCovar),
#          mod_NOcovar_hmgcr = map(mod_hmgcr, ~.x$model_noCovar),
#          mod_without = map(data, ~model_min(.x))
#   ) 
# 
# toc()
# beep(4)
# 
# # covariates overview
# res_mixMod_covar_01_tab <- res_mixMod_covar_01 |>
#   mutate(
#     mod_com_without = map(mod_without, ~.x$table),
#     mod_com_vek = map(mod_vek, ~.x$table),
#     mod_com_gk = map(mod_gk, ~.x$table),
#     mod_com_ck = map(mod_ck, ~.x$table),
#     mod_com_bcm = map(mod_bcm, ~.x$table),
#     mod_com_sex = map(mod_sex, ~.x$table),
#     mod_com_sub = map(mod_sub, ~.x$table),
#     mod_com_jo = map(mod_jo, ~.x$table),
#     mod_com_hmgcr = map(mod_bcm, ~.x$table)
#   ) |> 
#   select(starts_with("var_"), contains("mod_com")) |> 
#   unnest_wider(mod_com_without, names_sep = "_") |>
#   unnest_wider(mod_com_vek, names_sep = "_") |>
#   unnest_wider(mod_com_gk,  names_sep = "_")  |>
#   unnest_wider(mod_com_ck,  names_sep = "_") |> 
#   unnest_wider(mod_com_bcm,  names_sep = "_") |> 
#   unnest_wider(mod_com_sex,  names_sep = "_") |> 
#   unnest_wider(mod_com_sub,  names_sep = "_") |> 
#   unnest_wider(mod_com_jo,  names_sep = "_") |> 
#   unnest_wider(mod_com_hmgcr,  names_sep = "_") |> 
#   rename_with(~ str_remove(., "mod_com_"), everything()) |> 
#   mutate(across(where(is.numeric), ~round(.x, 3)))
# 
# # covariates overview kableExtra
# (res_mixMod_covar_01_tab_kable <- res_mixMod_covar_01_tab |> 
#   kable(escape = F, format = "html",
#         col.names = c("Dependent", "Independent", 
#                       "R2", "P-value - only independet", 
#                       "P-value - comparison", "R2", "P-Value - covariate",
#                       "P-value - comparison", "R2", "P-Value - covariate",
#                       "P-value - comparison", "R2", "P-Value - covariate",
#                       "P-value - comparison", "R2", "P-Value - covariate",                      "P-value - comparison", "R2", "P-Value - covariate",
#                       "P-value - comparison", "R2", "P-Value - covariate",
#                       "P-value - comparison", "R2", "P-Value - covariate",
#                       "P-value - comparison", "R2", "P-Value - covariate")) |> 
#   kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
#                 full_width = FALSE) |> 
#   add_header_above(c("Variables" = 2,  "Without covariate" = 2,
#                      "Age" = 3, "GC per day and weight" = 3,
#                      "Creatinine" = 3, "BCM" = 3, "Gender" = 3,
#                      "Subdiagnosis" = 3, "Anti-JO-1"= 3, "Anti-HMGCR" = 3)) |> 
#   collapse_rows(1, valign = "top") )
# 
# 
# 
# 
# res_mixMod_covar_02_vek <- res_mixMod_covar_01 |> 
#   select(starts_with("var_"), ends_with("_vek")) |> 
#   left_join(res_mixMod_covar_01_tab |> 
#               select(starts_with("var_"), vek_p_val_comp)) |> 
#   group_by(var_dep_name) |> 
#   mutate(
#     covariate = "age",
#     n_sig = sum(vek_p_val_comp < 0.06, na.rm = TRUE),
#     mod_final = if_else(n_sig > 2, mod_covar_vek, mod_NOcovar_vek),
#     adjustment  = if_else(n_sig > 2, "adjusted", "without"),
#     tidier = map(mod_final, ~model_parameters(
#       .x$fit,
#       ci = 0.95,
#       df_method = "satterthwaite"))
#   ) |>
#   ungroup() |>
#   unnest(tidier) |> 
#   filter(Parameter == "var_indep_value") |> 
#   select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
#          -CI, -Effects, -SE, -df_error) 
# 
# res_mixMod_covar_02_gk <- res_mixMod_covar_01 |> 
#   select(starts_with("var_"), ends_with("_gk")) |> 
#   left_join(res_mixMod_covar_01_tab |> 
#               select(starts_with("var_"), gk_p_val_comp)) |> 
#   group_by(var_dep_name) |> 
#   mutate(
#     covariate = "glucocoticoids",
#     n_sig = sum(gk_p_val_comp < 0.06, na.rm = TRUE),
#     mod_final = if_else(n_sig > 2, mod_covar_gk, mod_NOcovar_gk),
#     adjustment  = if_else(n_sig > 2, "adjusted", "without"),
#     tidier = map(mod_final, ~model_parameters(
#       .x$fit,
#       ci = 0.95,
#       df_method = "satterthwaite"))
#   ) |>
#   ungroup() |>
#   unnest(tidier) |> 
#   filter(Parameter == "var_indep_value") |> 
#   select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
#          -CI, -Effects, -SE, -df_error) 
# 
# res_mixMod_covar_02_ck <- res_mixMod_covar_01 |> 
#   select(starts_with("var_"), ends_with("_ck")) |> 
#   left_join(res_mixMod_covar_01_tab |> 
#               select(starts_with("var_"), ck_p_val_comp)) |> 
#   group_by(var_dep_name) |> 
#   mutate(
#     covariate = "creatinkinase",
#     n_sig = sum(ck_p_val_comp < 0.06, na.rm = TRUE),
#     mod_final = if_else(n_sig > 2, mod_covar_ck, mod_NOcovar_ck),
#     adjustment  = if_else(n_sig > 2, "adjusted", "without"),
#     tidier = map(mod_final, ~model_parameters(
#       .x$fit,
#       ci = 0.95,
#       df_method = "satterthwaite"))
#   ) |>
#   ungroup() |>
#   unnest(tidier) |> 
#   filter(Parameter == "var_indep_value") |> 
#   select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
#          -CI, -Effects, -SE, -df_error) 
# 
# 
# res_mixMod_covar_02_bmc <- res_mixMod_covar_01 |> 
#   select(starts_with("var_"), ends_with("_bcm")) |> 
#   left_join(res_mixMod_covar_01_tab |> 
#               select(starts_with("var_"), bcm_p_val_comp)) |> 
#   group_by(var_dep_name) |> 
#   mutate(
#     covariate = "body cell mass",
#     n_sig = sum(bcm_p_val_comp < 0.06, na.rm = TRUE),
#     mod_final = if_else(n_sig > 2, mod_covar_bcm, mod_NOcovar_bcm),
#     adjustment  = if_else(n_sig > 2, "adjusted", "without"),
#     tidier = map(mod_final, ~model_parameters(
#       .x$fit,
#       ci = 0.95,
#       df_method = "satterthwaite"))
#   ) |>
#   ungroup() |>
#   unnest(tidier) |> 
#   filter(Parameter == "var_indep_value") |> 
#   select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
#          -CI, -Effects, -SE, -df_error) 
# 
# res_mixMod_covar_02_sex <- res_mixMod_covar_01 |> 
#   select(starts_with("var_"), ends_with("_sex")) |> 
#   left_join(res_mixMod_covar_01_tab |> 
#               select(starts_with("var_"), sex_p_val_comp)) |> 
#   group_by(var_dep_name) |> 
#   mutate(
#     covariate = "gender",
#     n_sig = sum(sex_p_val_comp < 0.06, na.rm = TRUE),
#     mod_final = if_else(n_sig > 2, mod_covar_sex, mod_NOcovar_sex),
#     adjustment  = if_else(n_sig > 2, "adjusted", "without"),
#     tidier = map(mod_final, ~model_parameters(
#       .x$fit,
#       ci = 0.95,
#       df_method = "satterthwaite"))
#   ) |>
#   ungroup() |>
#   unnest(tidier) |> 
#   filter(Parameter == "var_indep_value") |> 
#   select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
#          -CI, -Effects, -SE, -df_error) 
# 
# res_mixMod_covar_02_sub <- res_mixMod_covar_01 |> 
#   select(starts_with("var_"), ends_with("_sub")) |> 
#   left_join(res_mixMod_covar_01_tab |> 
#               select(starts_with("var_"), sub_p_val_comp)) |> 
#   group_by(var_dep_name) |> 
#   mutate(
#     covariate = "subdiagnosis",
#     n_sig = sum(sub_p_val_comp < 0.06, na.rm = TRUE),
#     mod_final = if_else(n_sig > 2, mod_covar_sub, mod_NOcovar_sub),
#     adjustment  = if_else(n_sig > 2, "adjusted", "without"),
#     tidier = map(mod_final, ~model_parameters(
#       .x$fit,
#       ci = 0.95,
#       df_method = "satterthwaite"))
#   ) |>
#   ungroup() |>
#   unnest(tidier) |> 
#   filter(Parameter == "var_indep_value") |> 
#   select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
#          -CI, -Effects, -SE, -df_error) 
# 
# res_mixMod_covar_02_jo <- res_mixMod_covar_01 |> 
#   select(starts_with("var_"), ends_with("_jo")) |> 
#   left_join(res_mixMod_covar_01_tab |> 
#               select(starts_with("var_"), jo_p_val_comp)) |> 
#   group_by(var_dep_name) |> 
#   mutate(
#     covariate = "anti-jo-1",
#     n_sig = sum(jo_p_val_comp < 0.06, na.rm = TRUE),
#     mod_final = if_else(n_sig > 2, mod_covar_jo, mod_NOcovar_jo),
#     adjustment  = if_else(n_sig > 2, "adjusted", "without"),
#     tidier = map(mod_final, ~model_parameters(
#       .x$fit,
#       ci = 0.95,
#       df_method = "satterthwaite"))
#   ) |>
#   ungroup() |>
#   unnest(tidier) |> 
#   filter(Parameter == "var_indep_value") |> 
#   select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
#          -CI, -Effects, -SE, -df_error) 
# 
# 
# res_mixMod_covar_02_hmgcr <- res_mixMod_covar_01 |> 
#   select(starts_with("var_"), ends_with("_hmgcr")) |> 
#   left_join(res_mixMod_covar_01_tab |> 
#               select(starts_with("var_"), hmgcr_p_val_comp)) |> 
#   group_by(var_dep_name) |> 
#   mutate(
#     covariate = "anti-hmgcr",
#     n_sig = sum(hmgcr_p_val_comp < 0.06, na.rm = TRUE),
#     mod_final = if_else(n_sig > 2, mod_covar_hmgcr, mod_NOcovar_hmgcr),
#     adjustment  = if_else(n_sig > 2, "adjusted", "without"),
#     tidier = map(mod_final, ~model_parameters(
#       .x$fit,
#       ci = 0.95,
#       df_method = "satterthwaite"))
#   ) |>
#   ungroup() |>
#   unnest(tidier) |> 
#   filter(Parameter == "var_indep_value") |> 
#   select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
#          -CI, -Effects, -SE, -df_error) 
# 
# 
# export(list(overall = res_mixMod_covar_01_tab,
#             age = res_mixMod_covar_02_vek,
#             gk = res_mixMod_covar_02_gk,
#             ck = res_mixMod_covar_02_ck,
#             bmc = res_mixMod_covar_02_bmc,
#             gender = res_mixMod_covar_02_sex,
#             subdiagnosis = res_mixMod_covar_02_sub,
#             anti_jo_1 = res_mixMod_covar_02_jo,
#             anti_hmgcr = res_mixMod_covar_02_hmgcr),
#        "output/tables/250922_covariates_01.xlsx")
# 
# 
#  
#   
# test_data <- test_res$data[[1]]
# model_comp(test_data, "vek")
# 
# model_comp <- function(data, covar) {
#   
#   data <- data |> 
#     select(var_dep_value, poradie_vysetrenia, var_indep_value, all_of(covar), projekt_id) |> 
#     na.omit()
#   
#   rhs_terms <- c("poradie_vysetrenia", "var_indep_value", covar, "(1 | projekt_id)")
#   fml <- as.formula(paste("var_dep_value ~", paste(rhs_terms, collapse = " + ")))
#   
#   ## tidymodels
#   lmer_mod <- 
#     linear_reg() |> 
#     set_engine("lmer", REML = FALSE) |> 
#     set_mode("regression")
#   
#   mod1 <- fit(
#     lmer_mod,
#     formula = var_dep_value ~ poradie_vysetrenia + var_indep_value + (1 | projekt_id),
#     data = data
#   )
#   
#   mod2 <- fit(
#     lmer_mod,
#     formula = fml,
#     data = data
#   )
#   
#   p_val_comp <- stats::anova(mod1$fit, mod2$fit)$`Pr(>Chisq)`[[2]]
#   r2_mod <- performance(mod2, estimator = "ML")$R2_conditional
#   p_val_covar <- model_parameters(mod2$fit)$p[[6]]
#   
#   return(list(p_val_comp = p_val_comp,
#        r2_mod = r2_mod,
#        p_val_covar = p_val_covar))
#   
# }
# 
# 
# 
# 
# stats::anova(mod1$fit, mod2$fit)


## Veru requests ----
# sel_nam_01 <- sel_val_01 |> 
#   filter(!is.na(full_name)) |> 
#   pull(col)
# 
# data_m0_02 <- d03 |> 
#   group_by(projekt_id) |> 
#   mutate(across(mi_2:anti_hmgcr,
#                 ~ as.integer(any(.x == 1, na.rm = TRUE)))) |> 
#   ungroup() |> 
#   filter(poradie_vysetrenia == "M0" & !is.na(klastr_dle_exprese)) |> 
#   select(klastr_dle_exprese, all_of(sel_nam_01), -sae, -pl_7, -ej, -sas, -aba) |> 
#   relocate(odpoved_na_terapii_m0_vs_m6) |> 
#   # select(where(~mean(is.na(.x)) < 0.8)) |> 
#   mice::mice(method = "pmm", m = 5, printFlag = FALSE) |> 
#   mice::complete() |> 
#   select(where(~ !anyNA(.x)))
# 
# 
# res_pca <- prcomp(data_m0_02[, 3:ncol(data_m0_02)] |> 
#                     mutate(across(where(is.factor), as.numeric)),  
#                   scale = TRUE)
# fviz_pca_biplot(res_pca, label = "var", habillage=data_m0_02$klastr_dle_exprese,
#                 addEllipses=TRUE, ellipse.level=0.95,
#                 ggtheme = theme_minimal())
# fviz_pca_var(res_pca, col.var = "cos2",
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
#              repel = TRUE # Avoid text overlapping
# )
# 
# fviz_cos2(res_pca, choice = "var", axes = 1:2)
# 
# 
# data_m0_02 |> 
#   ggplot(aes(klastr_dle_exprese, ild)) +
#   geom_boxplot() +
#   geom_jitter() +
#   facet_wrap(~odpoved_na_terapii_m0_vs_m6)
# 
# 
# data_m0_02 |> 
#   ggplot(aes(klastr_dle_exprese, fill = odpoved_na_terapii_m0_vs_m6)) +
#   geom_bar(position = "fill")
# 
# 
# mosaicplot(data_m0_02 |> 
#              select(klastr_dle_exprese, odpoved_na_terapii_m0_vs_m6) |> 
#              table(),
#            main = "Mosaic plot",
#            color = TRUE
# )

# 250722 ----
## others - sex and age adjusted ----
### model ----
res_mixMod_age_sex <- d04_sel1_nested |> 
  filter(var_indep_name == "mstn" | !str_detect(var_dep_name, "mmt")) |>
  mutate(data = map(data, ~ .x |> 
                      mutate(var_indep_value_scl = log(var_indep_value + 1),
                             var_dep_value_scl = log(var_dep_value + 1),
                             gk = factor(gk)
                      )),
         mod_raw = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value*vek*pohlavi + (1 | projekt_id),
           data = .x
         )),
         mod_log_10 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value*vek*pohlavi + (1 | projekt_id),
           data = .x
         )),
         mod_log_01 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value_scl*vek*pohlavi + (1 | projekt_id),
           data = .x
         )),
         mod_log_11 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value_scl*vek*pohlavi + (1 | projekt_id),
           data = .x
         )),
         shapiro_p_raw = map_dbl(mod_raw, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_10 = map_dbl(mod_log_10, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_01 = map_dbl(mod_log_01, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_11 = map_dbl(mod_log_11, ~ shapiro.test(residuals(.x$fit))$p.value)
  )

res_mixMod_age_sex <- res_mixMod_age_sex |> 
  select(-contains("shap")) |> 
  pivot_longer(cols = contains("mod"),
               names_to = "mod_name",
               values_to = "mod_value") |> 
  mutate(mod_name = str_remove(mod_name, "mod_")) |> 
  left_join(res_mixMod_age_sex |> 
              select(-contains("mod_"), -data) |> 
              pivot_longer(cols = contains("shap"),
                           names_to = "shap_name",
                           values_to = "shap_value") |>
              mutate(shap_name = str_remove(shap_name, "shapiro_p_")),
            by = c("var_dep_name", "var_indep_name", "mod_name" = "shap_name")) |> 
  group_by(var_indep_name, var_dep_name) |> 
  slice_max(shap_value, n = 1) |> 
  ungroup() |> 
  mutate(tidier = map(mod_value, lmer_tidier),
         fig = pmap(list(mod_value, var_indep_name,var_dep_name), 
                    ~plot_mixed_terms_03(model = ..1, xlab = ..2, ylab = ..3, var_indep2 = "pohlavi"))) 


plot_mixed_terms_03(model = model, xlab = "test1", ylab = "test2")
plot_mixed_terms_04(model = model, xlab = "test1", ylab = "test2", data = data)

# export(res_mixMod_age_sex_tab, "output/tables/250722_others_age_sex_01.xlsx")

### figures ----
# pdf("output/figures/250722_others_age_sex_01.pdf", height = 6, width = 12)
walk(res_mixMod_age_sex$fig, print)
# dev.off()
# 

data <- res_mixMod_type$data[[1]]
model <- res_mixMod_type$mod_value[[1]]# davka_gk_mg_den

# podtyp nemoci + poradi vysetreni
plot_mixed_terms_03 <- function(model, 
                                var_pattern = "var_indep_",
                                var_indep2 = "podtyp_nemoci_zjednoduseny",
                                xlab        = NULL, 
                                ylab        = NULL,
                                data        = data) {
  # model        : a parsnip model_fit (with $fit = lmerMod) or a bare lmerMod
  # var_pattern  : regex/substring to match raw vs scaled predictor
  # xlab, ylab   : axis labels; NULL will default to the matched term / "Predicted response"
  
  library(ggeffects)
  library(ggplot2)
  library(broom.mixed)
  library(lmerTest)
  
  # 1) pull out the lmerMod
  engine <- if (inherits(model, "model_fit")) model$fit else model
  
  # 2) find the exact fixed‐effect term(s)
  fe      <- broom.mixed::tidy(as_lmerModLmerTest(engine), effects = "fixed")
  matches <- grep(var_pattern, fe$term, value = TRUE)
  if (length(matches) == 0) {
    stop("No fixed-effect term matches: ", var_pattern)
  }
  # prefer the scaled variant if present:
  predictor <- if ("var_indep_value_scl" %in% matches) {
    "var_indep_value_scl"
  } else matches[1]
  
  # 3) extract its p-value
  pval <- fe %>% 
    filter(term == predictor) %>% 
    pull(p.value) %>% 
    signif(2)
  
  
  
  # 4) get ggpredict data
  preds <- ggpredict(
    engine,
    terms = c(predictor, "podtyp_nemoci_zjednoduseny")
  )
  
  # 5) build the base plot
  fml <- as.formula(paste0("~", var_indep2))
  
  ct <- emtrends(engine,
                 fml,

              var =predictor) |>
    summary(infer = c(TRUE, TRUE)) |> 
    data.frame()
  
  preds2 <- ggpredict(
    engine,
    terms = c(predictor, "podtyp_nemoci_zjednoduseny", "poradie_vysetrenia")
  )
  
  raw_data <- attr(preds2, "rawdata")

  # 1) Build a lookup table of p‐values per facet
  pval_df <- ct %>%
    mutate(
      group   = as.character(.data[[var_indep2]]),         # vezmeme hodnoty ve sloupci "pohlavi"
      p.label = paste0("p = ", signif(p.value, 3))         # naformátujeme p‑hodnotu
    ) %>%
    select(group, p.label)
  
  data_new <- data
  dep_var <- deparse(extract_fit_engine(model)@call$formula)[1] |> 
    str_extract("^[^ ]+")

  colors_20 <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000",
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
    "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
    "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
  )
  
  p <- plot(preds) +
    geom_point(
      data = data_new,
      aes(x = !!sym(predictor), 
          y = !!sym(dep_var), 
          colour = odpoved_na_terapii_m0_vs_m6, 
          fill = odpoved_na_terapii_m0_vs_m6, 
          shape = odpoved_na_terapii_m0_vs_m6),
      size = 3,
      alpha       = 0.3,
      inherit.aes = FALSE
    ) +

    paletteer::scale_fill_paletteer_d("palettesForR::Cranes") +
    paletteer::scale_colour_paletteer_d("palettesForR::Cranes") +
    scale_shape_manual(values = c(21, 22, 23, 25,21, 22, 23, 25), name = "Examination") +
    # paletteer::scale_colour_paletteer_d("tvthemes::Steven") +
    # paletteer::scale_fill_paletteer_d("tvthemes::Steven") +
    facet_wrap(~ group, scales = "free_y") +
    labs(
      title  = paste0("p-value = ", pval),
      x      = xlab %||% predictor,
      y      = ylab %||% "Predicted response",
      colour = ""
    ) +
    theme_minimal(base_size = 18) +
    theme(
      strip.text       = element_text(face = "bold", size = 22),
      axis.title       = element_text(face = "bold", size = 22),
      axis.text        = element_text(size = 21),
      plot.title       = element_text(face = "bold", size = 25, hjust = 0.5),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
  
  # 3) Add a geom_text layer to annotate each facet
  p <- p + 
    geom_text(
      data = pval_df,
      aes(x = Inf, y = Inf, label = p.label),
      hjust = 1.1,    # nudge left from right edge
      vjust = 1.1,    # nudge down from top
      size  = 6, 
      inherit.aes = FALSE
    )
  
  return(p)
  
}




# 250722 ----
## others - sex and age adjusted ----
### model ----
# res_mixMod_age_sex <- d04_sel1_nested |> 
#   filter(var_indep_name == "mstn" | !str_detect(var_dep_name, "mmt")) |>
#   mutate(data = map(data, ~ .x |> 
#                       mutate(var_indep_value_scl = log(var_indep_value + 1),
#                              var_dep_value_scl = log(var_dep_value + 1),
#                              gk = factor(gk)
#                       )),
#          mod_raw = map(data, ~ fit(
#            lmer_mod,
#            formula = var_dep_value ~ poradie_vysetrenia + var_indep_value*vek*pohlavi + (1 | projekt_id),
#            data = .x
#          )),
#          mod_log_10 = map(data, ~ fit(
#            lmer_mod,
#            formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value*vek*pohlavi + (1 | projekt_id),
#            data = .x
#          )),
#          mod_log_01 = map(data, ~ fit(
#            lmer_mod,
#            formula = var_dep_value ~ poradie_vysetrenia + var_indep_value_scl*vek*pohlavi + (1 | projekt_id),
#            data = .x
#          )),
#          mod_log_11 = map(data, ~ fit(
#            lmer_mod,
#            formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value_scl*vek*pohlavi + (1 | projekt_id),
#            data = .x
#          )),
#          shapiro_p_raw = map_dbl(mod_raw, ~ shapiro.test(residuals(.x$fit))$p.value),
#          shapiro_p_log_10 = map_dbl(mod_log_10, ~ shapiro.test(residuals(.x$fit))$p.value),
#          shapiro_p_log_01 = map_dbl(mod_log_01, ~ shapiro.test(residuals(.x$fit))$p.value),
#          shapiro_p_log_11 = map_dbl(mod_log_11, ~ shapiro.test(residuals(.x$fit))$p.value)
#   )
# 
# res_mixMod_age_sex <- res_mixMod_age_sex |> 
#   select(-contains("shap")) |> 
#   pivot_longer(cols = contains("mod"),
#                names_to = "mod_name",
#                values_to = "mod_value") |> 
#   mutate(mod_name = str_remove(mod_name, "mod_")) |> 
#   left_join(res_mixMod_age_sex |> 
#               select(-contains("mod_"), -data) |> 
#               pivot_longer(cols = contains("shap"),
#                            names_to = "shap_name",
#                            values_to = "shap_value") |>
#               mutate(shap_name = str_remove(shap_name, "shapiro_p_")),
#             by = c("var_dep_name", "var_indep_name", "mod_name" = "shap_name")) |> 
#   group_by(var_indep_name, var_dep_name) |> 
#   slice_max(shap_value, n = 1) |> 
#   ungroup() |> 
#   mutate(tidier = map(mod_value, lmer_tidier),
#          fig = pmap(list(mod_value, var_indep_name,var_dep_name), 
#                     ~plot_mixed_terms_03(model = ..1, xlab = ..2, ylab = ..3, var_indep2 = "pohlavi"))) 
# 
# 
# ### table ----
# res_mixMod_age_sex_tab <- res_mixMod_age_sex |> 
#   unnest(tidier) |> 
#   filter(str_detect(term, "var_indep_")) |> 
#   select(var_indep_name, var_dep_name, mod_name, estimate, p.value,
#          conf.low, conf.high,  shap_value) |> 
#   mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
#   arrange(var_indep_name)
# 
# 
# # export(res_mixMod_age_sex_tab, "output/tables/250722_others_age_sex_01.xlsx")
# 
# ### figures ----
# # pdf("output/figures/250722_others_age_sex_01.pdf", height = 6, width = 12)
# walk(res_mixMod_age_sex$fig, print)
# # dev.off()
# # 
# 
# data <- res_mixMod_type$data[[1]]
# model <- res_mixMod_type$mod_value[[1]]# davka_gk_mg_den
# 
# # podtyp nemoci + poradi vysetreni
# plot_mixed_terms_03 <- function(model, 
#                                 var_pattern = "var_indep_",
#                                 var_indep2 = "podtyp_nemoci_zjednoduseny",
#                                 xlab        = NULL, 
#                                 ylab        = NULL) {
#   # model        : a parsnip model_fit (with $fit = lmerMod) or a bare lmerMod
#   # var_pattern  : regex/substring to match raw vs scaled predictor
#   # xlab, ylab   : axis labels; NULL will default to the matched term / "Predicted response"
#   
#   library(ggeffects)
#   library(ggplot2)
#   library(broom.mixed)
#   library(lmerTest)
#   
#   # 1) pull out the lmerMod
#   engine <- if (inherits(model, "model_fit")) model$fit else model
#   
#   # 2) find the exact fixed‐effect term(s)
#   fe      <- broom.mixed::tidy(as_lmerModLmerTest(engine), effects = "fixed")
#   matches <- grep(var_pattern, fe$term, value = TRUE)
#   if (length(matches) == 0) {
#     stop("No fixed-effect term matches: ", var_pattern)
#   }
#   # prefer the scaled variant if present:
#   predictor <- if ("var_indep_value_scl" %in% matches) {
#     "var_indep_value_scl"
#   } else matches[1]
#   
#   # 3) extract its p-value
#   pval <- fe %>% 
#     filter(term == predictor) %>% 
#     pull(p.value) %>% 
#     signif(2)
#   
#   
#   
#   # 4) get ggpredict data
#   preds <- ggpredict(
#     engine,
#     terms = c(predictor, var_indep2)
#   )
#   
#   # 5) build the base plot
#   fml <- as.formula(paste0("~", var_indep2))
#   
#   ct <- emtrends(engine,
#                  fml,
#               var =predictor) |>
#     summary(infer = c(TRUE, TRUE)) |> 
#     data.frame()
#   
#   preds2 <- ggpredict(
#     engine,
#     terms = c(predictor, var_indep2, "poradie_vysetrenia")
#   )
#   
#   raw_data <- attr(preds2, "rawdata")
#   
#   # 1) Build a lookup table of p‐values per facet
#   pval_df <- ct %>%
#     mutate(
#       group   = as.character(.data[[var_indep2]]),         # vezmeme hodnoty ve sloupci "pohlavi"
#       p.label = paste0("p = ", signif(p.value, 3))         # naformátujeme p‑hodnotu
#     ) %>%
#     select(group, p.label)
#   
#   
#   fml_new <- extract_fit_engine(model)@call$formula |> 
#     deparse() |> 
#     paste(collapse = " ") |> 
#     str_squish() |> 
#     str_replace("(\\+\\s*)podtyp_nemoci_zjednoduseny",
#                 "* podtyp_nemoci_zjednoduseny") |> 
#     str_replace("(\\+\\s*)var_indep_value",
#                 "* var_indep_value") |>
#     as.formula()
#   
#   data_new <- extract_fit_engine(model)@frame
#   
#   preds3 <- fit(lmer_mod, formula = fml_new, data = data_new)
#   
#   predictor2 <- paste0(predictor, " [all]")
#   
#   preds_int <- ggpredict(
#     preds3,
#     terms = c(predictor2, var_indep2)
#   )   
#   
#   
#   # 2) Your base plot
#   p <- plot(preds_int) +
#     geom_point(
#       data        = raw_data,
#       aes(x = x, y = response, colour = facet, fill = facet, shape = facet),
#       size = 3,
#       alpha       = 0.3,
#       inherit.aes = FALSE
#     ) +
#     paletteer::scale_fill_paletteer_d("palettesForR::Cranes") +
#     paletteer::scale_colour_paletteer_d("palettesForR::Cranes") +
#     scale_shape_manual(values = c(21, 22, 23, 25,21, 22, 23, 25), name = "Examination") +
#     # paletteer::scale_colour_paletteer_d("tvthemes::Steven") +
#     # paletteer::scale_fill_paletteer_d("tvthemes::Steven") +
#     facet_wrap(~ group, scales = "free_y") +
#     labs(
#       title  = paste0("p-value = ", pval),
#       x      = xlab %||% predictor,
#       y      = ylab %||% "Predicted response",
#       colour = ""
#     ) +
#     theme_minimal(base_size = 18) +
#     theme(
#       strip.text       = element_text(face = "bold", size = 22),
#       axis.title       = element_text(face = "bold", size = 22),
#       axis.text        = element_text(size = 21),
#       plot.title       = element_text(face = "bold", size = 25, hjust = 0.5),
#       panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
#       panel.grid.minor = element_blank()
#     )
#   
#   # 3) Add a geom_text layer to annotate each facet
#   p <- p + 
#     geom_text(
#       data = pval_df,
#       aes(x = Inf, y = Inf, label = p.label),
#       hjust = 1.1,    # nudge left from right edge
#       vjust = 1.1,    # nudge down from top
#       size  = 6, 
#       inherit.aes = FALSE
#     )
#   
#   return(p)
#   
# }
# 
# 
# 

# data <- res_mixMod_type$data[[1]]
# model <- data_test <- res_mixMod_type$mod_value[[1]]

# plot_mixed_terms_04 <- function(model, 
#                                 var_pattern = "var_indep_", 
#                                 xlab        = NULL, 
#                                 ylab        = NULL,
#                                 data        = data) {
#   # model        : a parsnip model_fit (with $fit = lmerMod) or a bare lmerMod
#   # var_pattern  : regex/substring to match raw vs scaled predictor
#   # xlab, ylab   : axis labels; NULL will default to the matched term / "Predicted response"
#   
#   library(ggeffects)
#   library(ggplot2)
#   library(broom.mixed)
#   library(lmerTest)
#   
#   # 1) pull out the lmerMod
#   engine <- if (inherits(model, "model_fit")) model$fit else model
#   
#   # 2) find the exact fixed‐effect term(s)
#   fe      <- broom.mixed::tidy(as_lmerModLmerTest(engine), effects = "fixed")
#   matches <- grep(var_pattern, fe$term, value = TRUE)
#   if (length(matches) == 0) {
#     stop("No fixed-effect term matches: ", var_pattern)
#   }
#   # prefer the scaled variant if present:
#   predictor <- if ("var_indep_value_scl" %in% matches) {
#     "var_indep_value_scl"
#   } else matches[1]
#   
#   # 3) extract its p-value
#   pval <- fe %>% 
#     filter(term == predictor) %>% 
#     pull(p.value) %>% 
#     signif(2)
#   
#   
#   
#   # 4) get ggpredict data
#   preds <- ggpredict(
#     engine,
#     terms = c(predictor, "podtyp_nemoci_zjednoduseny")
#   )
#   
#   # 5) build the base plot
#   ct <- emmeans(engine, "podtyp_nemoci_zjednoduseny") |> 
#     contrast(adjust="none") |> data.frame()
#   
#   preds2 <- ggpredict(
#     engine,
#     terms = c(predictor, "podtyp_nemoci_zjednoduseny", "poradie_vysetrenia")
#   )
#   
#   raw_data <- attr(preds2, "rawdata")
#   
#   # 1) Build a lookup table of p‐values per facet
#   pval_df <- ct %>% 
#     # strip off the common prefix so that it matches your 'group' levels
#     mutate(group = str_extract(contrast, "\\d+"),
#            p.label = paste0("p = ", signif(p.value, 3))) %>%
#     select(group, p.label)
#   
#   
#   fml_new <- extract_fit_engine(model)@call$formula |> 
#     deparse() |> 
#     paste(collapse = " ") |> 
#     str_squish() |> 
#     str_replace("(\\+\\s*)podtyp_nemoci_zjednoduseny",
#                 "* podtyp_nemoci_zjednoduseny") |> 
#     str_replace("(\\+\\s*)var_indep_value",
#                 "* var_indep_value") |>
#     as.formula()
#   
#   data_new <- data
#   
#   preds3 <- fit(lmer_mod, formula = fml_new, data = data_new)
#   
#   predictor2 <- paste0(predictor, " [all]")
#   
#   preds_int <- ggpredict(
#     preds3,
#     terms = c(predictor2, "podtyp_nemoci_zjednoduseny")
#   )   
#   
#   dep_var <- deparse(fml_new)[1] |> str_extract("^[^ ]+")
#   
#   # 2) Your base plot
#   p <- plot(preds_int) +
#     geom_point(
#       data = data_new,
#       aes(x = !!sym(predictor), 
#           y = !!sym(dep_var), 
#           colour = odpoved_na_terapii_m0_vs_m6, 
#           fill = odpoved_na_terapii_m0_vs_m6, 
#           shape = odpoved_na_terapii_m0_vs_m6),
#       size = 3,
#       alpha       = 0.3,
#       inherit.aes = FALSE
#     ) +
#     scale_fill_brewer(palette = "Dark2") +
#     scale_color_brewer(palette = "Dark2") +
#     scale_shape_manual(values = c(21, 22, 23, 24, 25,21, 22, 23, 24, 25), name = "Examination") +
#     # paletteer::scale_colour_paletteer_d("tvthemes::Steven") +
#     # paletteer::scale_fill_paletteer_d("tvthemes::Steven") +
#     facet_wrap(~ group, scales = "free_y") +
#     labs(
#       title  = paste0("p-value = ", pval),
#       x      = xlab %||% predictor,
#       y      = ylab %||% "Predicted response",
#       colour = ""
#     ) +
#     theme_minimal(base_size = 18) +
#     theme(
#       strip.text       = element_text(face = "bold", size = 22),
#       axis.title       = element_text(face = "bold", size = 22),
#       axis.text        = element_text(size = 21),
#       plot.title       = element_text(face = "bold", size = 25, hjust = 0.5),
#       panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
#       panel.grid.minor = element_blank()
#     )
#   
#   # 3) Add a geom_text layer to annotate each facet
#   p <- p + 
#     geom_text(
#       data = pval_df,
#       aes(x = Inf, y = Inf, label = p.label),
#       hjust = 1.1,    # nudge left from right edge
#       vjust = 1.1,    # nudge down from top
#       size  = 6, 
#       inherit.aes = FALSE
#     )
#   
#   return(p)
#   
# }







# 250716 ----
# 1) Re‑fit with an interaction
mod_int <- lmer(
  var_dep_value ~ var_indep_value * podtyp_nemoci_zjednoduseny
  + poradie_vysetrenia
  + (1 | projekt_id),
  data = data_test
)

# 2) Recompute the ggpredict grid
library(ggeffects)
preds_int <- ggpredict(
  mod_int,
  terms = c("var_indep_value [all]", "podtyp_nemoci_zjednoduseny")
)

# 3) Plot—now each facet gets its own slope
plot(preds_int, show_data = TRUE) +
  facet_wrap(~ group, scales = "free_y") +
  labs(
    x      = "var_indep_value",
    y      = "Predicted var_dep_value",
    colour = "Subtype"
  ) +
  theme_minimal()


model <- res_mixMod_type$mod_value[[50]]


plot_mixed_terms_03 <- function(model, 
                                var_pattern = "var_indep_", 
                                xlab        = NULL, 
                                ylab        = NULL) {
  # model        : a parsnip model_fit (with $fit = lmerMod) or a bare lmerMod
  # var_pattern  : regex/substring to match raw vs scaled predictor
  # xlab, ylab   : axis labels; NULL will default to the matched term / "Predicted response"
  
  library(ggeffects)
  library(ggplot2)
  library(broom.mixed)
  library(lmerTest)
  
  # 1) pull out the lmerMod
  engine <- if (inherits(model, "model_fit")) model$fit else model
  
  # 2) find the exact fixed‐effect term(s)
  fe      <- broom.mixed::tidy(as_lmerModLmerTest(engine), effects = "fixed")
  matches <- grep(var_pattern, fe$term, value = TRUE)
  if (length(matches) == 0) {
    stop("No fixed-effect term matches: ", var_pattern)
  }
  # prefer the scaled variant if present:
  predictor <- if ("var_indep_value_scl" %in% matches) {
    "var_indep_value_scl"
  } else matches[1]
  
  # 3) extract its p-value
  pval <- fe %>% 
    filter(term == predictor) %>% 
    pull(p.value) %>% 
    signif(2)
  
  
  
  # 4) get ggpredict data
  preds <- ggpredict(
    engine,
    terms = c(predictor, "podtyp_nemoci_zjednoduseny")
  )
  
  # 5) build the base plot
  ct <- emmeans(engine, "podtyp_nemoci_zjednoduseny") |> 
    contrast() |> data.frame()
  
  preds2 <- ggpredict(
    engine,
    terms = c(predictor, "podtyp_nemoci_zjednoduseny", "poradie_vysetrenia")
  )
  
  raw_data <- attr(preds2, "rawdata")
  
  # 1) Build a lookup table of p‐values per facet
  pval_df <- ct %>% 
    # strip off the common prefix so that it matches your 'group' levels
    mutate(group = str_extract(contrast, "\\d+"),
           p.label = paste0("p = ", signif(p.value, 3))) %>%
    select(group, p.label)
  
  
  fml_new <- extract_fit_engine(model)@call$formula |> 
    deparse() |> 
    paste(collapse = " ") |> 
    str_squish() |> 
    str_replace("(\\+\\s*)podtyp_nemoci_zjednoduseny",
                "* podtyp_nemoci_zjednoduseny") |> 
    as.formula()
  
  data_new <- extract_fit_engine(model)@frame
  
  preds3 <- fit(lmer_mod, formula = fml_new, data = data_new)
  
  predictor2 <- paste0(predictor, " [all]")
  
  preds_int <- ggpredict(
    preds3,
    terms = c(predictor2, "podtyp_nemoci_zjednoduseny")
  )   
  
  
  # 2) Your base plot
  p <- plot(preds_int) +
    geom_point(
      data        = raw_data,
      aes(x = x, y = response, colour = facet, fill = facet, shape = facet),
      size = 3,
      alpha       = 0.3,
      inherit.aes = FALSE
    ) +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(values = c(21, 22, 23, 25,21, 22, 23, 25), name = "Examination") +
    # paletteer::scale_colour_paletteer_d("tvthemes::Steven") +
    # paletteer::scale_fill_paletteer_d("tvthemes::Steven") +
    facet_wrap(~ group, scales = "free_y") +
    labs(
      title  = paste0("p-value = ", pval),
      x      = xlab %||% predictor,
      y      = ylab %||% "Predicted response",
      colour = ""
    ) +
    theme_minimal(base_size = 18) +
    theme(
      strip.text       = element_text(face = "bold", size = 22),
      axis.title       = element_text(face = "bold", size = 22),
      axis.text        = element_text(size = 21),
      plot.title       = element_text(face = "bold", size = 25, hjust = 0.5),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
  
  # 3) Add a geom_text layer to annotate each facet
  p <- p + 
    geom_text(
      data = pval_df,
      aes(x = Inf, y = Inf, label = p.label),
      hjust = 1.1,    # nudge left from right edge
      vjust = 1.1,    # nudge down from top
      size  = 5, 
      inherit.aes = FALSE
    )
  
  return(p)
  
}


# 250714 ----





mod_test <- res_mixMod_type$mod_value[[1]]







plot_model(mod_test, type = "emm", terms = "var_indep_value_scl", show.p = TRUE)

plot_mixed_terms_02 (mod_test, xlab="test1", ylab="test2")





# 1) grab fixed‐effect terms + p‐values
fe <- broom.mixed::tidy(as_lmerModLmerTest(mod_test$fit), effects = "fixed")

# 2) select those matching your pattern
sel_terms <- grep(var_pattern, fe$term, value = TRUE)
if (length(sel_terms) == 0) {
  stop("No fixed-effect terms matched pattern: ", var_pattern)
}
if (length(sel_terms) > 1) {
  warning("More than one term matched; plotting all. p-value will use the first.")
}

# 3) extract the p-value of the first matched term
pval <- fe %>% 
  filter(term == sel_terms[1]) %>% 
  pull(p.value)

# 4) build the emmeans plot
p <- sjPlot::plot_model(
  mod_test,
  type      = "emm",
  terms     = sel_terms,
  show.data = TRUE
)

# 5) apply axis labels if given
if (!is.null(xlab) || !is.null(ylab)) {
  p <- p + labs(
    x = xlab %||% waiver(),
    y = ylab %||% waiver()
  )
}

# 6) annotate p-value in top-right
p <- p +
  annotate(
    "text",
    x     = Inf,
    y     = Inf,
    label = paste0("p = ", signif(pval, 2)),
    hjust = 1.1,
    vjust = 1.5,
    size  = p_label_size
  ) +
  
  # 7) custom theme
  theme_sjplot2() +
  theme(
    plot.title       = element_blank(),
    plot.subtitle    = element_text(size = 20, color = "grey40", hjust = 0),
    axis.title       = element_text(face = "bold", size = 20),
    axis.text        = element_text(size = 18, color = "grey20"),
    panel.grid.major = element_line(linetype = "dotted", color = "grey80"),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )





         # mod_raw_tidier = map(mod_raw, lmer_tidier),
         # mod_log_10_tidier = map(mod_log_10, lmer_tidier),
         # mod_log_01_tidier = map(mod_log_01, lmer_tidier),
         # mod_log_11_tidier = map(mod_log_11, lmer_tidier)

res_mixMod_vek




















test <- res_mixMod_crp_ast_alt_ck |> 
  unnest(mod_raw_tidier) |> 
  filter(term == "var_indep_value") |> 
  select(var_dep_name, var_dep_name, estimate, p.value,
         conf.low, conf.high, shapiro_p_raw) 

data_test <- res_mixMod_crp_ast_alt_ck$data[[2]] |> na.omit()

mod_tmb <- glmmTMB(
  (var_dep_value+1) ~ var_indep_value + poradie_vysetrenia + (1 | projekt_id),
  data = data_test,
  family = Gamma(link = "log")
)

plot(ggeffects::ggpredict(mod_tmb, terms = "var_indep_value"))


data_test |> 
  ggplot(aes(var_indep_value, var_dep_value)) +
  geom_point()



#_________________________________________________________________
res_mixMod_crp_ast_alt_ck <- d04_sel1_nested |> 
  # filter(var_dep_name %in% c("ast", "alt", "ck", "crp")) |> 
  mutate(data = map(data, ~ .x |> 
                      mutate(var_indep_value_scl = scale(var_indep_value),
                             var_dep_value_scl = scale(var_dep_value)
                      )),
         mod_raw = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value + (1 | projekt_id),
           data = .x
         )),
         mod_log_10 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value + (1 | projekt_id),
           data = .x
         )),
         mod_log_01 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value_scl + (1 | projekt_id),
           data = .x
         )),
         mod_log_11 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value_scl + (1 | projekt_id),
           data = .x
         )),
         mod_raw_tidier = map(mod_raw, lmer_tidier),
         mod_log_10_tidier = map(mod_log_10, lmer_tidier),
         mod_log_01_tidier = map(mod_log_01, lmer_tidier),
         mod_log_11_tidier = map(mod_log_11, lmer_tidier),
         shapiro_p_raw = map_dbl(mod_raw, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_10 = map_dbl(mod_log_10, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_01 = map_dbl(mod_log_01, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_11 = map_dbl(mod_log_11, ~ shapiro.test(residuals(.x$fit))$p.value)
  )


test <- res_mixMod_crp_ast_alt_ck |> 
  unnest(mod_raw_tidier) |> 
  filter(term == "var_indep_value") |> 
  select(var_dep_name, var_dep_name, estimate, p.value,
         conf.low, conf.high, shapiro_p_raw) 



#______________________________________________________________________
res_mixMod_crp_ast_alt_ck <- d04_sel1_nested |> 
  # filter(var_dep_name %in% c("ast", "alt", "ck", "crp")) |> 
  mutate(data = map(data, ~ .x |> 
                      mutate(var_indep_value_scl = log(var_indep_value + 1),
                             var_dep_value_scl = log(var_dep_value + 1)
                             )),
         mod_raw = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value + (1 | projekt_id),
           data = .x
         )),
         mod_log_10 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value + (1 | projekt_id),
           data = .x
         )),
         mod_log_01 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value_scl + (1 | projekt_id),
           data = .x
         )),
         mod_log_11 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value_scl + (1 | projekt_id),
           data = .x
         )),
         mod_raw_tidier = map(mod_raw, lmer_tidier),
         mod_log_10_tidier = map(mod_log_10, lmer_tidier),
         mod_log_01_tidier = map(mod_log_01, lmer_tidier),
         mod_log_11_tidier = map(mod_log_11, lmer_tidier),
         shapiro_p_raw = map_dbl(mod_raw, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_10 = map_dbl(mod_log_10, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_01 = map_dbl(mod_log_01, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_11 = map_dbl(mod_log_11, ~ shapiro.test(residuals(.x$fit))$p.value)
         )


test <- res_mixMod_crp_ast_alt_ck |> 
  unnest(mod_raw_tidier) |> 
  filter(term == "var_indep_value") |> 
  select(var_dep_name, var_dep_name, estimate, p.value,
         conf.low, conf.high, shapiro_p_raw) 



res_mixMod_crp_ast_alt_ck$data[[10]] |> 
  # filter(var_indep_value <4000) |> 
  ggplot(aes(var_indep_value, var_dep_value, col = factor(podtyp_nemoci))) + 
  geom_point() +
  stat_smooth(method = "lm", se = FALSE) +
  facet_wrap(~poradie_vysetrenia)





recipe_01 <- recipe(var_dep_value ~ poradie_vysetrenia + var_indep_value + (1|projekt_id), 
                 test) |> 
  step_log(var_indep_value, var_dep_name, mod_log_10, shapiro_p_log_10) |> 
  

lmer_wflow <- 
  workflow() |> 
  add_recipe(recipe(var_dep_value ~ poradie_vysetrenia + var_indep_value , data = test)) %>%
  # or use
  # add_variables(outcomes = Reaction, predictors = c(Days, Subject))
  add_model(lmer_mod, formula = var_dep_value ~ poradie_vysetrenia + var_indep_value  + (1 | projekt_id))




# This skips workflow, and works directly:
mod_test <- fit(
  lmer_mod,
  formula = var_dep_value ~ poradie_vysetrenia + var_indep_value + (1|projekt_id),
  data = test |> mutate(var_dep_value = log(var_dep_value+1))
)

mod_test |> broom.mixed::tidy()

check_model(mod_test)














plot_dbscan_pca_umap <- function(data_numeric, groups, eps = 0.5, minPts = 5) {
  # Checks
  stopifnot(nrow(data_numeric) == length(groups))
  
  # Packages
  require(ggplot2)
  require(dbscan)
  require(umap)
  require(gridExtra)
  
  # Run DBSCAN
  db <- dbscan::dbscan(data_numeric, eps = eps, minPts = minPts)
  cluster <- as.factor(db$cluster)  # 0 means noise
  
  group_factor <- as.factor(groups)
  
  # PCA
  pca <- prcomp(data_numeric, scale. = TRUE)
  pca_df <- as.data.frame(pca$x[, 1:2])
  colnames(pca_df) <- c("PC1", "PC2")
  pca_df$Cluster <- cluster
  pca_df$Group <- group_factor
  
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, shape = Group)) +
    geom_point(size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "PCA + DBSCAN", color = "DBSCAN Cluster", shape = "Group")
  
  # UMAP
  umap_result <- umap::umap(data_numeric)
  umap_df <- as.data.frame(umap_result$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$Cluster <- cluster
  umap_df$Group <- group_factor
  
  p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster, shape = Group)) +
    geom_point(size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "UMAP + DBSCAN", color = "DBSCAN Cluster", shape = "Group")
  
  # Combine
  gridExtra::grid.arrange(p_pca, p_umap, nrow = 1)
}



plot_dbscan_pca_umap(data_m0, group_m0, eps = 0.1, minPts = 2)




pacman::p_load(dbscan, factoextra)

hdbscan(data_m0, minPts = 5)


cl <- hdbscan(data_m0, minPts = 5)
plot(data_m0, col=cl$cluster+1, pch=20)

pacman::p_load(fpc, factoextra)

db <- fpc::dbscan(data_m0, eps = 0.15, MinPts = 5)
plot(db, data_m0, main = "DBSCAN", frame = FALSE)


db <- dbscan(data_m0, eps = 0.5, minPts = 5)
df$cluster <- as.factor(db$cluster)

pca_data <- as.data.frame(pca_m0$x[,1:2])
pca_data$cluster <- df$cluster
pca_data$group <- df$group




###############################################################################


# 
# 
# mod_mixMODmmt(d04_sel1_nested$data[[5]], 100)
# 
# data <- d04_sel1_nested$data[[5]]
# 
# 
# mod_mixMODmmt(d04_sel1_nested$data[[1]], 80)
# 
# 
# 
# data <- d04_sel1_nested$data[[5]]
# upper = 100






# fig_mixMODmmt <- function(data, mod, var_indep_name, var_dep_name, upper) {
#   data_set <- data |> 
#     mutate(
#       var_indep_value_scl = log(var_indep_value),
#       y_capped = pmin(var_dep_value, upper),
#       ind = as.integer(var_dep_value >= upper)
#     )
#   
#   pred_df <- GLMMadaptive::effectPlotData(mod, newdata = na.omit(data_set))
#   
#   fig <- ggplot(pred_df, aes(x = var_indep_value, y = var_dep_value)) +
#     geom_point(aes(color = poradie_vysetrenia), alpha = 0.6, show.legend = FALSE) +
#     geom_line(aes(y = pred, color = poradie_vysetrenia), size = 1, show.legend = FALSE) +
#     geom_ribbon(aes(ymax = upp, ymin = low, fill = poradie_vysetrenia), 
#                 alpha = 0.3, linetype = 0) +
#     facet_wrap(~ poradie_vysetrenia, scales = "free_y") +
#     scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
#     labs(
#       x = var_indep_name,
#       y = var_dep_name,
#       fill = "Time point"
#     ) +
#     theme_sjplot2() +
#     theme(
#       axis.title = element_text(size = 14, face = "bold"),
#       strip.text = element_text(face = "bold"),
#       strip.background = element_rect(fill = "gray90", color = NA)
#     )
#   
#   return(fig)
# }
# 
