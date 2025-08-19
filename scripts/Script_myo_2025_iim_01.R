# Header 1 ----
## Header 2 ----
### Header 3 ----

# References ----
# model diagnosis
## https://easystats.github.io/see/articles/performance.html
## https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
# statistical data and model report
## report(model) # report(model)
## https://github.com/easystats/report
# Multivariate analyses
## https://cran.r-project.org/web/packages/MVN/vignettes/MVN.pdf
# Regression model summary
## https://strengejacke.github.io/sjPlot/articles/tab_model_estimates.html
## https://cran.r-project.org/web/packages/jtools/vignettes/summ.html
# Regression model plots
## https://cran.r-project.org/web/packages/sjPlot/vignettes/plot_model_estimates.html
# Markdown session info
## https://stackoverflow.com/questions/33156946/how-can-i-format-sessioninfo-in-rmarkdown

# Directories ----
dir.create("data/raw", recursive = TRUE)
dir.create("data/processed", recursive = TRUE)
dir.create("scripts/functions", recursive = TRUE)
dir.create("output/figures", recursive = TRUE)
dir.create("output/tables", recursive = TRUE)
dir.create("misc")
dir.create("reports")



# Libraries & functions ----
# required
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, purrr,conflicted)

# Functions and libraries uploading
list.files("scripts/functions/", pattern="*.*", full.names=TRUE) |> map(~source(.))

# Back-up
back_up("scripts/functions/FUN_01.R") # the the destination subdirectory specify using 'path_dest'
back_up("scripts/functions/OBJ_01.R") # the the destination subdirectory specify using 'path_dest'
back_up("scripts/Script_myo_2025_iim_01.R") # the the destination subdirectory specify using 'path_dest'

# Save to .RData
save(xx1,
     file = "reports/markD_01.RData")
save.image("reports/markD_01.RData")
# Save the tables into data/tables.RData by listing them individually
cgwtools::resave(xx2,xx3, file = "reports/markD_01.RData") # resave a list of tables that I'll use in the .Rmd file.
# Save the tables into data/tables.RData using "patterns" 
cgwtools::resave(list=ls(pattern="tbl"), file = "reports/markD_01.RData")

# copy-to-clipboard botton
# ```{r xaringanExtra-clipboard, echo=FALSE}
# xaringanExtra::use_clipboard()
# ```

# EDA ----
## Duplicates ----
eval_dup_01 <- d03 |> get_dupes()  |>  capture_messages()

## Missing values ----
### all dataset ----
vis_miss(d03)
gg_miss_var(d03)
gg_miss_var(d03 ,facet=poradie_vysetrenia)
gg_miss_case(d03)
gg_miss_case(d03 ,facet = poradie_vysetrenia)
(miss_tab01 <- miss_var_summary(d03))
(miss_tab02 <- miss_case_summary(d03))

gg_miss_upset(d03) # give an overall pattern of missingness
gg_miss_fct(x = d03, fct = poradie_vysetrenia)

## selection ----
vis_miss(d04_sel1)
gg_miss_var(d04_sel1)
gg_miss_var(d04_sel1 ,facet=poradie_vysetrenia)
gg_miss_case(d04_sel1)
gg_miss_case(d04_sel1 ,facet = poradie_vysetrenia)
(miss_tab01 <- miss_var_summary(d04_sel1))
(miss_tab02 <- miss_case_summary(d04_sel1))

gg_miss_upset(d04_sel1) # give an overall pattern of missingness
gg_miss_fct(x = d04_sel1, fct = poradie_vysetrenia)


## overview ----
d03 |>
  select(-projekt_id) |> 
  my_skim()

## correlation ----
### ggally ----
lower = 0.6; upper <- 0.7; examination = "M0"
cor_data <- d03 |> 
  filter(poradie_vysetrenia == examination) |> 
  select(where(is.numeric)) |> 
  select(where(~mean(is.na(.x)) < 0.8))

cor_mat <- cor(cor_data, method = "kendall", use = "pairwise.complete.obs")
cor_mat[is.na(cor_mat)] <- 0 # nahradí všechny NA za 0

cor_df <- reshape2::melt(cor_mat)
names(cor_df) <- c("Var1", "Var2", "value")

param_name_sig <- cor_df |>
  filter(
    as.numeric(Var1) > as.numeric(Var2),
    abs(value) > lower,
    abs(value) < upper
  )|> 
  select(-value) |> 
  pivot_longer(cols = everything(),
               values_to = "param_name_sig") |> 
  distinct(param_name_sig) |> 
  pull(param_name_sig) |> 
  droplevels()

ggpairs(cor_data,
        columns = cor_data |> select(all_of(param_name_sig)) |> names(), 
        lower = list(continuous = my_fn, method = "kendall")) +
  theme_sjplot2()

### selected variable ----
eda_06 <- d04_sel1 |> 
  select(!contains("mmt"), -projekt_id) |> 
  group_by(poradie_vysetrenia) |> 
  correlate(all_of(var_indep_01)) |> 
  filter(!var2 %in% var_indep_01) |> 
  plot()

## heatmap ----
# m0
heatmap.2(scale(data_m0), scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")
# m3
heatmap.2(scale(data_m3), scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")
# m6
heatmap.2(scale(data_m6), scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")
# m18
heatmap.2(scale(data_m18), scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")


## PCA ----
# m0
pca_m0 <- prcomp(data_m0, 
                 scale = TRUE)
fviz_eig(pca_m0)
fviz_pca_ind(pca_m0, addEllipses=TRUE, ellipse.level=0.95, habillage = group_m0)

# m3
pca_m3 <- prcomp(data_m3, 
                 scale = TRUE)
fviz_eig(pca_m3)
fviz_pca_ind(pca_m3, addEllipses=TRUE, ellipse.level=0.95, habillage = group_m3)

# m6
pca_m6 <- prcomp(data_m6, 
                 scale = TRUE)
fviz_eig(pca_m6)
fviz_pca_ind(pca_m6, addEllipses=TRUE, ellipse.level=0.95, habillage = group_m6)

# m18
pca_m18 <- prcomp(data_m18, 
                 scale = TRUE)
fviz_eig(pca_m18)
fviz_pca_ind(pca_m18, addEllipses=TRUE, ellipse.level=0.95, habillage = group_m18)

## UMAP ----
# m0
# uwot::umap(pca_m0$x) |> plot(main = "UMAP - M0")
plot_umap(data_m0, col = group_m0, pca = pca_m0$x, main = "UMAP - M0")
# m3
# uwot::umap(pca_m3$x) |> plot(main = "UMAP - M3")
plot_umap(data_m3, col = group_m3, pca = pca_m3$x, main = "UMAP - M3")
# m6
# uwot::umap(pca_m6$x) |> plot(main = "UMAP - M6")
plot_umap(data_m6, col = group_m6, pca = pca_m6$x, main = "UMAP - M6")
# m18
# uwot::umap(pca_m18$x) |> plot(main = "UMAP - M18")
plot_umap(data_m18, col = group_m18, pca = pca_m18$x, main = "UMAP - M18")

## pca + umap + dbscan ----

## densitogram ----
eda_01 <- d04_sel1 |> 
  select(!contains("mmt"), -projekt_id, -odpoved_na_terapii_m0_vs_m6) |> 
  pivot_longer(cols = fi_2:last_col(),
               names_to = "param_name",
               values_to = "param_value") |> 
  ggplot(aes(param_value, fill=poradie_vysetrenia)) +
  geom_density(alpha = 0.1) +
  facet_wrap(~param_name, scales = "free")

eda_02 <- d04_sel1 |> 
  select(!contains("mmt"), -projekt_id, -poradie_vysetrenia) |> 
  relocate(odpoved_na_terapii_m0_vs_m6, .before = fi_2) |>
  pivot_longer(cols = fi_2:last_col(),
               names_to = "param_name",
               values_to = "param_value") |> 
  ggplot(aes(param_value, fill=odpoved_na_terapii_m0_vs_m6)) +
  geom_density(alpha = 0.1) +
  facet_wrap(~param_name, scales = "free")


## scatter plot ----
eda_03 <- d04_sel1 |> 
  select(!contains("mmt"), -projekt_id) |> 
  mutate(poradie_vysetrenia = str_remove(poradie_vysetrenia, "M") |> as.numeric()) |> 
  relocate(odpoved_na_terapii_m0_vs_m6, poradie_vysetrenia, .before = fi_2) |> 
  pivot_longer(cols = fi_2:last_col(),
               names_to = "param_name",
               values_to = "param_value") |> 
  ggplot(aes(poradie_vysetrenia, param_value, col=odpoved_na_terapii_m0_vs_m6)) +
  geom_jitter() +
  stat_smooth(method="loess", se=FALSE) +
  facet_wrap(~param_name, scales = "free")

eda_04 <- d04_sel1 |> 
  select(!contains("mmt"), -projekt_id, -all_of(var_indep_01[2:length(var_indep_01)])) |> 
  relocate(odpoved_na_terapii_m0_vs_m6, poradie_vysetrenia, mstn,.before = fi_2) |> 
  pivot_longer(cols = fi_2:last_col(),
               names_to = "param_name",
               values_to = "param_value") |> 
  ggplot(aes(log(mstn), param_value, col=poradie_vysetrenia)) +
  geom_point() +
  stat_smooth(method="loess", se=FALSE) +
  theme_sjplot2() +
  facet_wrap(~param_name, scales = "free")


## boxplot ----
eda_05 <- d04_sel1 |> 
  select(!contains("mmt"), -projekt_id, -odpoved_na_terapii_m0_vs_m6) |> 
  pivot_longer(cols = fi_2:last_col(),
               names_to = "param_name",
               values_to = "param_value") |> 
  ggplot(aes(poradie_vysetrenia, param_value, col=poradie_vysetrenia)) +
  geom_boxplot() +
  facet_wrap(~param_name, scales = "free")




# Statistics ----
## mmt -----
### model ----
# model
res_mixMod_mmt <-  d04_sel1_nested |> 
  filter(str_detect(var_dep_name, "mmt8")) |> 
  mutate(mod = map(data, ~possmixMODmmt(.x, upper = 80)),
         tidier = map(mod, broom.mixed::tidy),
         fig = pmap(list(data, mod, var_indep_name, var_dep_name, upper = 80),fig_mixMODmmt)
  ) |> 
  bind_rows(d04_sel1_nested |> 
              filter(str_detect(var_dep_name, "mmt10")) |> 
              mutate(mod = map(data, ~possmixMODmmt(.x, upper = 100)),
                     tidier = map(mod, broom.mixed::tidy),
                     fig = pmap(list(data, mod, var_indep_name, var_dep_name, upper = 100),fig_mixMODmmt)
              ))


data_set <- prepare_for_censored_model(d04_sel1_nested$data[[5]],
                                       
                                       dep_var = "var_dep_value",
                                       indep_var = "var_indep_value",
                                       group_var = "projekt_id",
                                       upper = 100
)

data_mod <- GLMMadaptive::mixed_model(
  fixed  = cbind(y_capped, ind) ~ poradie_vysetrenia + var_indep_value_scl,
  random = ~ 1 | projekt_id,
  family = GLMMadaptive::censored.normal(),
  data   = data_set
)
# data for mmt_10 and mstn analyses has to be specific for this analysis 
data_pred <- GLMMadaptive::effectPlotData(data_mod, newdata = na.omit(data_set))
data_fig <- ggplot(pred_df, aes(x = var_indep_value, y = var_dep_value)) +
  geom_point(aes(color = poradie_vysetrenia), alpha = 0.6, show.legend = FALSE) +
  geom_line(aes(y = pred, color = poradie_vysetrenia), size = 1, show.legend = FALSE) +
  geom_ribbon(aes(ymax = upp, ymin = low, fill = poradie_vysetrenia), 
              alpha = 0.3, linetype = 0) +
  facet_wrap(~ poradie_vysetrenia, scales = "free_y") +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(
    x = "mstn",
    y = "mmt10_total",
    fill = "Time point"
  ) +
  theme_sjplot2() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90", color = NA)
  )

rx_mmt_mstn <- which(d04_sel1_nested$var_dep_name == "mmt10_total" &
                       d04_sel1_nested$var_indep_name == "mstn")
res_mixMod_mmt$data[[rx_mmt_mstn]] <- data_set
res_mixMod_mmt$mod[[rx_mmt_mstn]] <- data_mod
res_mixMod_mmt$tidier [[rx_mmt_mstn]] <- data_mod |> tidy()
res_mixMod_mmt$fig[[rx_mmt_mstn]] <- data_fig

### figures ----
# pdf("output/figures/250704_mmt_01.pdf", width = 10, height = 6.5)
walk(res_mixMod_mmt$fig, print)
# dev.off()

### table ----
res_mixMod_mmt_tab <- res_mixMod_mmt |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep")) |> 
  select(-data, -mod, -term, -std.error, -statistic, -fig)

# export(res_mixMod_mmt_tab, "output/tables/250704_mmt_01.xlsx")

# kable
res_mixMod_mmt_tab |> 
  mutate(p.value = ifelse(p.value < 0.05, 
                          ifelse(p.value < 0.001,
                                 "<b><0.001</b>", 
                                 sprintf("<b>%.3f</b>", p.value)), 
                          sprintf("%.3f", p.value))) |> 
  kable(escape = F, format = "html",
        digits = 3,
        col.names = c("Dependent variable", "Independent variable", "Estimate",
                      "p-value", "2.5% CI", "97.5% CI")) |> 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                full_width = FALSE) |> 
  collapse_rows(1, valign = "top") 

## others - without adjustment ----
### model ----
res_mixMod_raw <- d04_sel1_nested |> 
  # filter(var_dep_name %in% c("ast", "alt", "ck", "crp")) |> 
  filter(!str_detect(var_dep_name, "mmt")) |>
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
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value_scl +  (1 | projekt_id),
           data = .x
         )),
         mod_log_11 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value_scl + (1 | projekt_id),
           data = .x
         )),
         shapiro_p_raw = map_dbl(mod_raw, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_10 = map_dbl(mod_log_10, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_01 = map_dbl(mod_log_01, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_11 = map_dbl(mod_log_11, ~ shapiro.test(residuals(.x$fit))$p.value)
  )

res_mixMod_raw <- res_mixMod_raw |> 
  select(-contains("shap")) |> 
  pivot_longer(cols = contains("mod"),
               names_to = "mod_name",
               values_to = "mod_value") |> 
  mutate(mod_name = str_remove(mod_name, "mod_")) |> 
  left_join(res_mixMod_raw |> 
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
                    ~plot_mixed_terms_01(model = ..1, xlab = ..2, ylab = ..3))) 


### table ----
res_mixMod_raw_tab <- res_mixMod_raw |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_")) |> 
  select(var_indep_name, var_dep_name, mod_name, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name )


# export(res_mixMod_raw_tab, "output/tables/250714_others_raw_01.xlsx")

### figures ----
pdf("output/figures/250714_others_raw_01.pdf")
walk(res_mixMod_raw$fig, print)
dev.off()


## others - age adjusted ----
### model ----
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
  group_by(var_indep_name, var_dep_name) |> 
  slice_max(shap_value, n = 1) |> 
  ungroup() |> 
  mutate(tidier = map(mod_value, lmer_tidier),
         fig = pmap(list(mod_value, var_indep_name,var_dep_name), 
                    ~plot_mixed_terms_01(model = ..1, xlab = ..2, ylab = ..3))) 


### table ----
res_mixMod_vek_tab <- res_mixMod_vek |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_")) |> 
  select(var_indep_name, var_dep_name, mod_name, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name )


# export(res_mixMod_vek_tab, "output/tables/250714_others_vek_01.xlsx")

### figures ----
# pdf("output/figures/250714_others_vek_01.pdf")
walk(res_mixMod_vek$fig, print)
# dev.off()


## others - podnemoc adjusted ----
### model ----
res_mixMod_type <- d04_sel1_nested |> 
  # filter(var_dep_name %in% c("ast", "alt", "ck", "crp")) |> 
  filter(!str_detect(var_dep_name, "mmt")) |>
  mutate(data = map(data, ~ .x |> 
                      mutate(var_indep_value_scl = log(var_indep_value + 1),
                             var_dep_value_scl = log(var_dep_value + 1),
                             podtyp_nemoci_zjednoduseny = factor(podtyp_nemoci_zjednoduseny)
                      )),
         mod_raw = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value*podtyp_nemoci_zjednoduseny + (1 | projekt_id),
           data = .x
         )),
         mod_log_10 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value*podtyp_nemoci_zjednoduseny + (1 | projekt_id),
           data = .x
         )),
         mod_log_01 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value_scl*podtyp_nemoci_zjednoduseny +  (1 | projekt_id),
           data = .x
         )),
         mod_log_11 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value_scl*podtyp_nemoci_zjednoduseny + (1 | projekt_id),
           data = .x
         )),
         shapiro_p_raw = map_dbl(mod_raw, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_10 = map_dbl(mod_log_10, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_01 = map_dbl(mod_log_01, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_11 = map_dbl(mod_log_11, ~ shapiro.test(residuals(.x$fit))$p.value)
  )

res_mixMod_type <- res_mixMod_type |> 
  select(-contains("shap")) |> 
  pivot_longer(cols = contains("mod"),
               names_to = "mod_name",
               values_to = "mod_value") |> 
  mutate(mod_name = str_remove(mod_name, "mod_")) |> 
  left_join(res_mixMod_type |> 
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
                    ~plot_mixed_terms_03(model = ..1, xlab = ..2, ylab = ..3)),
         fig2 = pmap(list(mod_value, var_indep_name,var_dep_name, data), 
                    ~plot_mixed_terms_04(model = ..1, xlab = ..2, ylab = ..3, data = ..4))) 


### table ----
res_mixMod_type_tab <- res_mixMod_type |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_") & !str_detect(term, ":")) |>
  select(var_indep_name, var_dep_name, mod_name, term, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name )


# export(res_mixMod_type_tab, "output/tables/250722_others_type_01.xlsx")

### figures ----
# podtyp nemoci + poradi vysetreni
pdf("output/figures/250722_others_type_01.pdf", height = 8, width = 16)
walk(res_mixMod_type$fig, print)
dev.off()

#podtyp nemoci + odpoved po 6 mesicic
pdf("output/figures/250722_others_type_02.pdf", height = 8, width = 16)
safe_print <- possibly(print, otherwise = NULL)
walk(res_mixMod_type$fig2, ~ safe_print(.x))
dev.off()

## others - glucocorticoids adjusted ----
### model ----
res_mixMod_gk <- d04_sel1_nested |> 
  # filter(var_dep_name %in% c("ast", "alt", "ck", "crp")) |> 
  filter(!str_detect(var_dep_name, "mmt")) |>
  mutate(data = map(data, ~ .x |> 
                      mutate(var_indep_value_scl = log(var_indep_value + 1),
                             var_dep_value_scl = log(var_dep_value + 1),
                             gk = factor(gk)
                      )),
         mod_raw = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value + gk + (1 | projekt_id),
           data = .x
         )),
         mod_log_10 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value + gk + (1 | projekt_id),
           data = .x
         )),
         mod_log_01 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value_scl +  gk + (1 | projekt_id),
           data = .x
         )),
         mod_log_11 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value_scl +  gk + (1 | projekt_id),
           data = .x
         )),
         shapiro_p_raw = map_dbl(mod_raw, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_10 = map_dbl(mod_log_10, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_01 = map_dbl(mod_log_01, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_11 = map_dbl(mod_log_11, ~ shapiro.test(residuals(.x$fit))$p.value)
  )

res_mixMod_gk <- res_mixMod_gk |> 
  select(-contains("shap")) |> 
  pivot_longer(cols = contains("mod"),
               names_to = "mod_name",
               values_to = "mod_value") |> 
  mutate(mod_name = str_remove(mod_name, "mod_")) |> 
  left_join(res_mixMod_gk |> 
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
                    ~plot_mixed_terms_01(model = ..1, xlab = ..2, ylab = ..3))) 


### table ----
res_mixMod_gk_tab <- res_mixMod_gk_dose |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_")) |> 
  select(var_indep_name, var_dep_name, mod_name, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name )


# export(res_mixMod_gk_tab, "output/tables/250714_others_gk_01.xlsx")

### figures ----
# pdf("output/figures/250714_others_gk_01.pdf")
walk(res_mixMod_gk$fig, print)
# dev.off()

## others - glucocorticoids dose adjusted ----
### model ----
res_mixMod_gk_dose <- d04_sel1_nested |> 
  # filter(var_dep_name %in% c("ast", "alt", "ck", "crp")) |> 
  filter(!str_detect(var_dep_name, "mmt")) |>
  mutate(data = map(data, ~ .x |> 
                      mutate(var_indep_value_scl = log(var_indep_value + 1),
                             var_dep_value_scl = log(var_dep_value + 1),
                             gk = factor(gk)
                      )),
         mod_raw = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value + davka_gk_mg_den + (1 | projekt_id),
           data = .x
         )),
         mod_log_10 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value + davka_gk_mg_den + (1 | projekt_id),
           data = .x
         )),
         mod_log_01 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value_scl +  davka_gk_mg_den + (1 | projekt_id),
           data = .x
         )),
         mod_log_11 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value_scl +  davka_gk_mg_den + (1 | projekt_id),
           data = .x
         )),
         shapiro_p_raw = map_dbl(mod_raw, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_10 = map_dbl(mod_log_10, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_01 = map_dbl(mod_log_01, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_11 = map_dbl(mod_log_11, ~ shapiro.test(residuals(.x$fit))$p.value)
  )

res_mixMod_gk_dose <- res_mixMod_gk_dose |> 
  select(-contains("shap")) |> 
  pivot_longer(cols = contains("mod"),
               names_to = "mod_name",
               values_to = "mod_value") |> 
  mutate(mod_name = str_remove(mod_name, "mod_")) |> 
  left_join(res_mixMod_gk_dose |> 
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
                    ~plot_mixed_terms_01(model = ..1, xlab = ..2, ylab = ..3))) 


### table ----
res_mixMod_gk_dose_tab <- res_mixMod_gk_dose |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_")) |> 
  select(var_indep_name, var_dep_name, mod_name, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name )


# export(res_mixMod_gk_dose_tab, "output/tables/250722_others_gk_dose_01.xlsx")

### figures ----
# pdf("output/figures/250722_others_gk_dose_01.pdf")
walk(res_mixMod_gk_dose$fig, print)
# dev.off()

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


### table ----
res_mixMod_age_sex_tab <- res_mixMod_age_sex |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_")) |> 
  select(var_indep_name, var_dep_name, mod_name, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name)



# export(res_mixMod_age_sex_tab, "output/tables/250722_others_age_sex_01.xlsx")

### figures ----
# pdf("output/figures/250722_others_age_sex_01.pdf", height = 6, width = 12)
walk(res_mixMod_age_sex$fig, print)
# dev.off()


