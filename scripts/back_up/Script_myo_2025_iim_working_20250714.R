# Back-up
back_up("scripts/functions/FUN_01.R") # the the destination subdirectory specify using 'path_dest'
back_up("scripts/functions/OBJ_01.R") # the the destination subdirectory specify using 'path_dest'
back_up("scripts/Script_myo_2025_iim_working.R") # the the destination subdirectory specify using 'path_dest'


# 250714 ----


res_mixMod_vek <- d04_sel1_nested |> 
  # filter(var_dep_name %in% c("ast", "alt", "ck", "crp")) |> 
  filter(!str_detect(var_dep_name, "mmt")) |>
  mutate(data = map(data, ~ .x |> 
                      mutate(var_indep_value_scl = log(var_indep_value + 1),
                             var_dep_value_scl = log(var_dep_value + 1)
                      )),
         mod_raw = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value + vek + (1 | projekt_id),
           data = .x
         )),
         mod_log_10 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value + vek + (1 | projekt_id),
           data = .x
         )),
         mod_log_01 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value_scl +  vek + (1 | projekt_id),
           data = .x
         )),
         mod_log_11 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value_scl +  vek + (1 | projekt_id),
           data = .x
         )),
         shapiro_p_raw = map_dbl(mod_raw, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_10 = map_dbl(mod_log_10, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_01 = map_dbl(mod_log_01, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_11 = map_dbl(mod_log_11, ~ shapiro.test(residuals(.x$fit))$p.value)
  )

res_mixMod_vek <- res_mixMod_vek |> 
  select(-contains("shap")) |> 
  pivot_longer(cols = contains("mod"),
               names_to = "mod_name",
               values_to = "mod_value") |> 
  mutate(mod_name = str_remove(mod_name, "mod_")) |> 
  left_join(res_mixMod_vek |> 
              select(-contains("mod_"), -data) |> 
              pivot_longer(cols = contains("shap"),
                           names_to = "shap_name",
                           values_to = "shap_value") |>
              mutate(shap_name = str_remove(shap_name, "shapiro_p_")),
            by = c("var_dep_name", "var_indep_name", "mod_name" = "shap_name")) |> 
  group_by(var_dep_name, var_indep_name) |> 
  slice_max(shap_value, n = 1) |> 
  ungroup() |> 
  mutate(tidier = map(mod_value, lmer_tidier),
         fig = pmap(list(mod_value, var_indep_name,var_dep_name), 
                    ~plot_mixed_terms_01(model = ..1, xlab = ..2, ylab = ..3)))

# table
res_mixMod_vek_tab <- res_mixMod_vek |> 
  unnest(tidier) |> 
  filter(term == "var_indep_value") |> 
  select(var_dep_name, var_indep_name, estimate, p.value,
         conf.low, conf.high,  shap_value) 

# figures 
walk(res_mixMod_vek$fig)



mod_test <- res_mixMod_vek$mod_value[[1]]


plot_model(mod_test, type = "emm", terms = "var_indep_value_scl", show.p = TRUE)

plot_mixed_terms(mod, xlab="test1", ylab="test2")











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
