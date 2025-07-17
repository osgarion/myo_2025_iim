# Back-up
back_up("scripts/functions/FUN_01.R") # the the destination subdirectory specify using 'path_dest'
back_up("scripts/functions/OBJ_01.R") # the the destination subdirectory specify using 'path_dest'
back_up("scripts/Script_myo_2025_iim_working.R") # the the destination subdirectory specify using 'path_dest'

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
