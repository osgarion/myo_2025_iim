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
pacman::p_load(tidyverse, purrr, conflicted, renv)

# Functions and libraries uploading
list.files("scripts/functions/", pattern="*.*", full.names=TRUE) |> map(~source(.))

# Environment reproducibility ----
# renv::init()        # once
# renv::status()      # sync?
# renv::snapshot()    # after changing deps
# renv::restore()     # on a new machine

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
gg_miss_var(d04_sel2 ,facet=poradie_vysetrenia)
gg_miss_case(d03)
gg_miss_case(d03 ,facet = poradie_vysetrenia)
gg_miss_case(d04_sel2 ,facet = poradie_vysetrenia)
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

eda_07 <- d03 |>
  select(all_of(var_indep_01), podtyp_nemoci_zjednoduseny) |>
  filter(!is.na(podtyp_nemoci_zjednoduseny)) |>
  ggpairs(ggplot2::aes(colour=podtyp_nemoci_zjednoduseny),
          lower = list(continuous = wrap(lowerFn, method = "lm")))

pdf("output/figures/250820_corr_indep_01.pdf")
print(eda_07)
dev.off()


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

## podtyp across periods ----
eda_period_type_01 <- d04_sel1_nested |> 
  ungroup() |> 
  select(-var_indep_name) |>
  distinct(var_dep_name, .keep_all = TRUE) |> 
  mutate(fig = pmap(list(data, var_dep_name),
                    ~plot_period_type_01(data=..1, ylab=..2)))

pdf("output/figures/250820_podtyp_period_01.pdf", width = 14)
walk(eda_period_type_01$fig, print)
dev.off()

## lmer overview ----
eda_lmer_overview_01 <- d04_sel1_nested |> 
  ungroup() |> 
  mutate(fig = pmap(list(data, var_indep_name, var_dep_name),
                    ~plot_lmer_grid_basic_01(data=..1, xlab=..2, ylab=..3)))

# pdf("output/figures/250820_lmer_overview_01.pdf", width = 14)
walk(eda_lmer_overview_01$fig, print)
# dev.off()

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

## features without adjustment ----
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
export(res_mixMod_raw_tab, "output/tables/250820_others_raw_01.xlsx")

### figures ----
# pdf("output/figures/250714_others_raw_01.pdf")
walk(res_mixMod_raw$fig, print)
# dev.off()


## feature adjustment ----
### age ----
#### model ----
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


#### table ----
res_mixMod_vek_tab <- res_mixMod_vek |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_")) |> 
  select(var_indep_name, var_dep_name, mod_name, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name )


# export(res_mixMod_vek_tab, "output/tables/250714_others_vek_01.xlsx")
export(res_mixMod_vek_tab, "output/tables/250820_others_vek_01.xlsx")

#### figures ----
# pdf("output/figures/250714_others_vek_01.pdf")
pdf("output/figures/250820_others_vek_01.pdf")
walk(res_mixMod_vek$fig, print)
dev.off()


### podnemoc ----
#### model ----
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


#### table ----
res_mixMod_type_tab <- res_mixMod_type |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_") & !str_detect(term, ":")) |>
  select(var_indep_name, var_dep_name, mod_name, term, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name )

# export(res_mixMod_type_tab, "output/tables/250722_others_type_01.xlsx")
export(res_mixMod_type_tab, "output/tables/250820_others_type_01.xlsx")

#### figures ----
# podtyp nemoci + poradi vysetreni
# pdf("output/figures/250722_others_type_01.pdf", height = 8, width = 16)
pdf("output/figures/250820_others_type_01.pdf", height = 8, width = 16)
walk(res_mixMod_type$fig, print)
dev.off()

#podtyp nemoci + odpoved po 6 mesicic
# pdf("output/figures/250722_others_type_02.pdf", height = 8, width = 16)
pdf("output/figures/250820_others_type_02.pdf", height = 8, width = 16)
safe_print <- possibly(print, otherwise = NULL)
walk(res_mixMod_type$fig2, ~ safe_print(.x))
dev.off()

### podnemoc - simple ----
res_mixMod_type_simple <- d04_sel1_nested |> 
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

res_mixMod_type_simple <- res_mixMod_type_simple |> 
  select(-contains("shap")) |> 
  pivot_longer(cols = contains("mod"),
               names_to = "mod_name",
               values_to = "mod_value") |> 
  mutate(mod_name = str_remove(mod_name, "mod_")) |> 
  left_join(res_mixMod_type_simple |> 
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


#### table ----
res_mixMod_type_simple_tab <- res_mixMod_type_simple |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_") & !str_detect(term, ":")) |>
  select(var_indep_name, var_dep_name, mod_name, term, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name )

export(res_mixMod_type_simple_tab, "output/tables/250819_others_type_simple_01.xlsx")


#### figures ----
# podtyp nemoci + poradi vysetreni
pdf("output/figures/250819_others_type_simple_01.pdf", height = 8, width = 16)
walk(res_mixMod_type_simple$fig, print)
dev.off()

#podtyp nemoci + odpoved po 6 mesicic
pdf("output/figures/250819_others_type_simple_02.pdf", height = 8, width = 16)
safe_print <- possibly(print, otherwise = NULL)
walk(res_mixMod_type_simple$fig2, ~ safe_print(.x))
dev.off()

### mstn-nomr podnemoc - simple ----
#### model ----
res_mixMod_mstnNorm_type_simple <- d06_norm_merge_nest |> 
  # filter(var_dep_name %in% c("ast", "alt", "ck", "crp")) |> 
  filter(!str_detect(var_dep_name, "mmt")) |>
  mutate(data = map(data, ~ .x |> 
                      mutate(var_indep_value_scl = log(var_indep_value + 1),
                             var_dep_value_scl = log(var_dep_value + 1),
                             gk = factor(gk)
                      )),
         mod_raw = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ var_indep_value + poradie_vysetrenia + podtyp_nemoci_zjednoduseny + (1 | projekt_id),
           data = .x
         )),
         mod_log_10 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ var_indep_value + poradie_vysetrenia + podtyp_nemoci_zjednoduseny + (1 | projekt_id),
           data = .x
         )),
         mod_log_01 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ var_indep_value_scl + poradie_vysetrenia + podtyp_nemoci_zjednoduseny + (1 | projekt_id),
           data = .x
         )),
         mod_log_11 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ var_indep_value_scl + poradie_vysetrenia + podtyp_nemoci_zjednoduseny + (1 | projekt_id),
           data = .x
         )),
         shapiro_p_raw = map_dbl(mod_raw, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_10 = map_dbl(mod_log_10, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_01 = map_dbl(mod_log_01, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_11 = map_dbl(mod_log_11, ~ shapiro.test(residuals(.x$fit))$p.value)
  )

res_mixMod_mstnNorm_type_simple <- res_mixMod_mstnNorm_type_simple |>
  select(-contains("shap")) |> 
  pivot_longer(cols = contains("mod"),
               names_to = "mod_name",
               values_to = "mod_value") |> 
  mutate(mod_name = str_remove(mod_name, "mod_")) |> 
  left_join(res_mixMod_mstnNorm_type_simple |> 
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
                    ~plot_mixed_terms_05(model = ..1, xlab = ..2, ylab = ..3)), 
         fig2 = pmap(list(mod_value, var_indep_name,var_dep_name, data), 
                     ~plot_mixed_terms_04(model = ..1, xlab = ..2, ylab = ..3, data = ..4)))


#### table ----
res_mixMod_mstnNorm_type_simple_tab <- res_mixMod_mstnNorm_type_simple |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_")) |> 
  select(var_indep_name, var_dep_name, mod_name, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name )

# export(res_mixMod_mstnNorm_type_simple_tab, "output/tables/250722_others_mstnNorm_type_simple_01.xlsx")
export(res_mixMod_mstnNorm_type_simple_tab, "output/tables/250820_others_mstnNorm_type_simple_01.xlsx")



#### figures ----
# podtyp nemoci + poradi vysetreni
pdf("output/figures/250819_others_mstnNorm_type_simple_01.pdf", height = 8, width = 16)
walk(res_mixMod_mstnNorm_type_simple$fig, print)
dev.off()

#podtyp nemoci + odpoved po 6 mesicic
pdf("output/figures/250819_others_mstnNorm_type_simple_02.pdf", height = 8, width = 16)
safe_print <- possibly(print, otherwise = NULL)
walk(res_mixMod_mstnNorm_type_simple$fig2, ~ safe_print(.x))
dev.off()

### glucocorticoids ----
#### model ----
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


#### table ----
res_mixMod_gk_tab <- res_mixMod_gk |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_")) |> 
  select(var_indep_name, var_dep_name, mod_name, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name )


# export(res_mixMod_gk_tab, "output/tables/250714_others_gk_01.xlsx")
export(res_mixMod_gk_tab, "output/tables/250820_others_gk_01.xlsx")

#### figures ----
# pdf("output/figures/250714_others_gk_01.pdf")
pdf("output/figures/250820_others_gk_01.pdf")
walk(res_mixMod_gk$fig, print)
dev.off()

### glucocorticoids dose ----
#### model ----
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


#### table ----
res_mixMod_gk_dose_tab <- res_mixMod_gk_dose |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_")) |> 
  select(var_indep_name, var_dep_name, mod_name, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name )

# export(res_mixMod_gk_dose_tab, "output/tables/250722_others_gk_dose_01.xlsx")
export(res_mixMod_gk_dose_tab, "output/tables/250820_others_gk_dose_01.xlsx")



#### figures ----
# pdf("output/figures/250722_others_gk_dose_01.pdf")
pdf("output/figures/250820_others_gk_dose_01.pdf")
walk(res_mixMod_gk_dose$fig, print)
dev.off()

###  glucocorticoids to bw dose ----
#### model ----
res_mixMod_gk_bw_dose <- d04_sel1_nested |> 
  # filter(var_dep_name %in% c("ast", "alt", "ck", "crp")) |> 
  filter(!str_detect(var_dep_name, "mmt")) |>
  mutate(data = map(data, ~ .x |> 
                      mutate(var_indep_value_scl = log(var_indep_value + 1),
                             var_dep_value_scl = log(var_dep_value + 1),
                             gk = factor(gk)
                      )),
         mod_raw = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value + denna_davka_gk_mg_kg_bw + (1 | projekt_id),
           data = .x
         )),
         mod_log_10 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value + denna_davka_gk_mg_kg_bw + (1 | projekt_id),
           data = .x
         )),
         mod_log_01 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ poradie_vysetrenia + var_indep_value_scl +  denna_davka_gk_mg_kg_bw + (1 | projekt_id),
           data = .x
         )),
         mod_log_11 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ poradie_vysetrenia + var_indep_value_scl +  denna_davka_gk_mg_kg_bw + (1 | projekt_id),
           data = .x
         )),
         shapiro_p_raw = map_dbl(mod_raw, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_10 = map_dbl(mod_log_10, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_01 = map_dbl(mod_log_01, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_11 = map_dbl(mod_log_11, ~ shapiro.test(residuals(.x$fit))$p.value)
  )

res_mixMod_gk_bw_dose <- res_mixMod_gk_bw_dose |> 
  select(-contains("shap")) |> 
  pivot_longer(cols = contains("mod"),
               names_to = "mod_name",
               values_to = "mod_value") |> 
  mutate(mod_name = str_remove(mod_name, "mod_")) |> 
  left_join(res_mixMod_gk_bw_dose |> 
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


#### table ----
res_mixMod_gk_bw_dose_tab <- res_mixMod_gk_bw_dose |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_|denna_davka")) |>
  # select(var_indep_name, var_dep_name, mod_name, estimate, p.value,
  #        conf.low, conf.high,  shap_value) |>
  select(var_indep_name, var_dep_name, term, mod_name, estimate, p.value,
         conf.low, conf.high,  shap_value) |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name )

export(res_mixMod_gk_bw_dose_tab, "output/tables/250819_others_gk_bw_dose_01.xlsx")

#### figures ----
pdf("output/figures/250819_others_gk_bw_dose_01.pdf")
walk(res_mixMod_gk_bw_dose$fig, print)
dev.off()

### mstn-nomr glucocorticoids dose ----
#### model ----
res_mixMod_mstnNorm_gk_dose <- d06_norm_merge_nest |> 
  # filter(var_dep_name %in% c("ast", "alt", "ck", "crp")) |> 
  filter(!str_detect(var_dep_name, "mmt")) |>
  mutate(data = map(data, ~ .x |> 
                      mutate(var_indep_value_scl = log(var_indep_value + 1),
                             var_dep_value_scl = log(var_dep_value + 1),
                             gk = factor(gk)
                      )),
         mod_raw = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ var_indep_value + poradie_vysetrenia + davka_gk_mg_den + (1 | projekt_id),
           data = .x
         )),
         mod_log_10 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ var_indep_value + poradie_vysetrenia + davka_gk_mg_den + (1 | projekt_id),
           data = .x
         )),
         mod_log_01 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ var_indep_value_scl + poradie_vysetrenia + davka_gk_mg_den + (1 | projekt_id),
           data = .x
         )),
         mod_log_11 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ var_indep_value_scl + poradie_vysetrenia + davka_gk_mg_den + (1 | projekt_id),
           data = .x
         )),
         shapiro_p_raw = map_dbl(mod_raw, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_10 = map_dbl(mod_log_10, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_01 = map_dbl(mod_log_01, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_11 = map_dbl(mod_log_11, ~ shapiro.test(residuals(.x$fit))$p.value)
  )

res_mixMod_mstnNorm_gk_dose <- res_mixMod_mstnNorm_gk_dose |>
  select(-contains("shap")) |> 
  pivot_longer(cols = contains("mod"),
               names_to = "mod_name",
               values_to = "mod_value") |> 
  mutate(mod_name = str_remove(mod_name, "mod_")) |> 
  left_join(res_mixMod_mstnNorm_gk_dose |> 
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
                    ~plot_mixed_terms_05(model = ..1, xlab = ..2, ylab = ..3))) 


#### table ----
res_mixMod_mstnNorm_gk_dose_tab <- res_mixMod_mstnNorm_gk_dose |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_")) |> 
  select(var_indep_name, var_dep_name, mod_name, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name )

# export(res_mixMod_mstnNorm_gk_dose_tab, "output/tables/250722_others_mstnNorm_gk_dose_01.xlsx")
export(res_mixMod_mstnNorm_gk_dose_tab, "output/tables/250820_others_mstnNorm_gk_dose_01.xlsx")



#### figures ----
pdf("output/figures/250820_others_mstnNorm_gk_dose_01.pdf")
walk(res_mixMod_mstnNorm_gk_dose$fig, print)
dev.off()

### creatine dose ----
#### model ----
res_mixMod_creatine_dose <- d03 |>
  mutate(
    odpoved_na_terapii_m0_vs_m6 =
      str_replace(odpoved_na_terapii_m0_vs_m6, "No Improvement", "No improvement")
  ) |>
  select(projekt_id, pohlavi, vek, gk, davka_gk_mg_den,denna_davka_gk_mg_kg_bw, podtyp_nemoci_zjednoduseny, poradie_vysetrenia, all_of(var_dep_01), all_of(var_indep_01)) |>
  filter(poradie_vysetrenia %in% c("M0", "M3", "M6", "M18")) |>
  pivot_longer(cols = any_of(var_dep_01 |>  str_subset("^(odpoved_na_terapii_m0_vs_m6|kreatinin_umol_l)$", negate = TRUE)),
               names_to = "var_dep_name",
               values_to = "var_dep_value") |> 
  pivot_longer(cols = any_of(var_indep_01),
               names_to = "var_indep_name",
               values_to = "var_indep_value") |> 
  group_by(var_dep_name, var_indep_name) |> 
  nest() |> 
  # filter(var_dep_name %in% c("ast", "alt", "ck", "crp")) |> 
  filter(!str_detect(var_dep_name, "mmt")) |>
  mutate(data = map(data, ~ .x |> 
                      mutate(var_indep_value_scl = log(var_indep_value + 1),
                             var_dep_value_scl = log(var_dep_value + 1),
                             gk = factor(gk)
                      )),
         mod_raw = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ var_indep_value + poradie_vysetrenia + kreatinin_umol_l + (1 | projekt_id),
           data = .x
         )),
         mod_log_10 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ var_indep_value + poradie_vysetrenia + kreatinin_umol_l + (1 | projekt_id),
           data = .x
         )),
         mod_log_01 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value ~ var_indep_value_scl + poradie_vysetrenia +  kreatinin_umol_l + (1 | projekt_id),
           data = .x
         )),
         mod_log_11 = map(data, ~ fit(
           lmer_mod,
           formula = var_dep_value_scl ~ var_indep_value_scl + poradie_vysetrenia +  kreatinin_umol_l + (1 | projekt_id),
           data = .x
         )),
         shapiro_p_raw = map_dbl(mod_raw, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_10 = map_dbl(mod_log_10, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_01 = map_dbl(mod_log_01, ~ shapiro.test(residuals(.x$fit))$p.value),
         shapiro_p_log_11 = map_dbl(mod_log_11, ~ shapiro.test(residuals(.x$fit))$p.value)
  )

res_mixMod_creatine_dose <- res_mixMod_creatine_dose |> 
  select(-contains("shap")) |> 
  pivot_longer(cols = contains("mod"),
               names_to = "mod_name",
               values_to = "mod_value") |> 
  mutate(mod_name = str_remove(mod_name, "mod_")) |> 
  left_join(res_mixMod_creatine_dose |> 
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
                    ~plot_mixed_terms_05(model = ..1, xlab = ..2, ylab = ..3))) 


#### table ----
res_mixMod_creatine_dose_tab <- res_mixMod_creatine_dose |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_|kreatin")) |> 
  select(var_indep_name, var_dep_name, term, mod_name, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name )

# export(res_mixMod_creatine_dose_tab, "output/tables/250722_others_creatine_dose_01.xlsx")
export(res_mixMod_creatine_dose_tab, "output/tables/250820_others_creatine_dose_01.xlsx")



#### figures ----
# pdf("output/figures/250722_others_creatine_dose_01.pdf")
pdf("output/figures/250820_others_creatine_dose_01.pdf")
walk(res_mixMod_creatine_dose$fig, print)
dev.off()

### sex and age ----
#### model ----
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


#### table ----
res_mixMod_age_sex_tab <- res_mixMod_age_sex |> 
  unnest(tidier) |> 
  filter(str_detect(term, "var_indep_")) |> 
  select(var_indep_name, var_dep_name, mod_name, estimate, p.value,
         conf.low, conf.high,  shap_value) |> 
  mutate(across(where(is.numeric), ~ round(.x, 3))) |> 
  arrange(var_indep_name)

# export(res_mixMod_age_sex_tab, "output/tables/250722_others_age_sex_01.xlsx")
export(res_mixMod_age_sex_tab, "output/tables/250820_others_age_sex_01.xlsx")


#### figures ----
# pdf("output/figures/250722_others_age_sex_01.pdf", height = 6, width = 12)
pdf("output/figures/250820_others_age_sex_01.pdf", height = 6, width = 12)
walk(res_mixMod_age_sex$fig, print)
dev.off()

## covariate models ----
### models ----
tic()
Sys.time()

cluster <- multidplyr::new_cluster(ncores)
multidplyr::cluster_library(cluster, c("dplyr","purrr","stringr","tibble", 
                                       "data.table", "tidyverse", "tidymodels",
                                       "multilevelmod","lme4","lmerTest",
                                       "broom.mixed", "easystats"))
multidplyr::cluster_copy(cluster, c("model_comp","model_min", "model_comp_int"))


res_mixMod_covar_01 <- d04_sel2 |> 
  mutate(across(.cols = c("ast", "alt", "ck", "crp", "haq", "ld", "mitax", "myoact", "myoglobin"), function(x) log(x+1))) |> 
  pivot_longer(cols = any_of(var_dep_01 |>  str_subset("odpoved_na_terapii_m0_vs_m6|kreatinin_umol_l", negate = TRUE)),
               names_to = "var_dep_name",
               values_to = "var_dep_value") |> 
  pivot_longer(cols = any_of(var_indep_01),
               names_to = "var_indep_name",
               values_to = "var_indep_value") |> 
  group_by(var_dep_name, var_indep_name) |> 
  # filter(!str_detect(var_dep_name, "_total")) |>           # drop of mmt-8 and mmt-10
  multidplyr::partition(cluster) |> 
  do(
    # data = data.table(.),
    mod_vek =model_comp(.,"vek"),
    mod_gk = model_comp(.,"denna_davka_gk_mg_kg_bw"),
    mod_ck = model_comp(.,"kreatinin_umol_l"),
    mod_bcm = model_comp(.,"bcm"),
    mod_sex = model_comp(.,"pohlavi"),
    mod_sub = model_comp(.,"podtyp_nemoci_zjednoduseny"),
    mod_jo = model_comp(.,"jo_1"),
    mod_hmgcr = model_comp(.,"anti_hmgcr"),
    mod_int_sex = model_comp_int(data = ., "vek"),
    mod_without = model_min(.)
  ) |> 
  collect() |> 
  mutate(
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
    mod_NOcovar_int_sex = map(mod_int_sex, ~.x$model_noCovar)
  ) 

walk(cluster, ~ .x$kill()); gc()

toc()
beep(4)

### covariates overview ----
res_mixMod_covar_01_tab <- res_mixMod_covar_01 |>
  mutate(
    mod_com_without = map(mod_without, ~.x$table),
    mod_com_vek = map(mod_vek, ~.x$table),
    mod_com_gk = map(mod_gk, ~.x$table),
    mod_com_ck = map(mod_ck, ~.x$table),
    mod_com_bcm = map(mod_bcm, ~.x$table),
    mod_com_sex = map(mod_sex, ~.x$table),
    mod_com_sub = map(mod_sub, ~.x$table),
    mod_com_jo = map(mod_jo, ~.x$table),
    mod_com_hmgcr = map(mod_hmgcr, ~.x$table),
    mod_com_int_sex = map(mod_int_sex, ~.x$table)
  ) |> 
  select(starts_with("var_"), contains("mod_com")) |> 
  unnest_wider(mod_com_without, names_sep = "_") |>
  unnest_wider(mod_com_vek, names_sep = "_") |>
  unnest_wider(mod_com_gk,  names_sep = "_")  |>
  unnest_wider(mod_com_ck,  names_sep = "_") |> 
  unnest_wider(mod_com_bcm,  names_sep = "_") |> 
  unnest_wider(mod_com_sex,  names_sep = "_") |> 
  unnest_wider(mod_com_sub,  names_sep = "_") |> 
  unnest_wider(mod_com_jo,  names_sep = "_") |> 
  unnest_wider(mod_com_hmgcr,  names_sep = "_") |> 
  unnest_wider(mod_com_int_sex,  names_sep = "_") |> 
  rename_with(~ str_remove(., "mod_com_"), everything()) |> 
  mutate(across(where(is.numeric), ~round(.x, 3))) |> 
  arrange(var_dep_name, var_indep_name) 

#### kableExtra ----
(res_mixMod_covar_01_tab_kable <- res_mixMod_covar_01_tab |> 
    kable(escape = F, format = "html",
          col.names = c("Dependent", "Independent", 
                        "R2", "P-value - only independet", 
                        "P-value - comparison", "R2", "P-Value - covariate",
                        "P-value - comparison", "R2", "P-Value - covariate",
                        "P-value - comparison", "R2", "P-Value - covariate",
                        "P-value - comparison", "R2", "P-Value - covariate",                      "P-value - comparison", "R2", "P-Value - covariate",
                        "P-value - comparison", "R2", "P-Value - covariate",
                        "P-value - comparison", "R2", "P-Value - covariate",
                        "P-value - comparison", "R2", "P-Value - covariate",
                        "P-value - comparison", "R2", "P-Value - covariate")) |> 
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                  full_width = FALSE) |> 
    add_header_above(c("Variables" = 2,  "Without covariate" = 2,
                       "Age" = 3, "GC per day and weight" = 3,
                       "Creatinine" = 3, "BCM" = 3, "Gender" = 3,
                       "Subdiagnosis" = 3, "Anti-JO-1"= 3, "Anti-HMGCR" = 3,
                       "With Gender Interaction" = 3)) |> 
    collapse_rows(1, valign = "top") )

### model comparisons ----
res_mixMod_covar_02_vek <- res_mixMod_covar_01 |> 
  select(starts_with("var_"), ends_with("_vek")) |> 
  left_join(res_mixMod_covar_01_tab |> 
              select(starts_with("var_"), vek_p_val_comp)) |> 
  group_by(var_dep_name) |> 
  mutate(
    covariate = "age",
    n_sig = sum(vek_p_val_comp < 0.06, na.rm = TRUE),
    mod_final = if_else(n_sig > 2, mod_covar_vek, mod_NOcovar_vek),
    adjustment  = if_else(n_sig > 2, "adjusted", "without"),
    tidier = map(mod_final, ~model_parameters(
      .x$fit,
      ci = 0.95,
      df_method = "satterthwaite"))
  ) |>
  ungroup() |>
  unnest(tidier) |> 
  filter(Parameter == "var_indep_value") |> 
  select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
         -CI, -Effects, -SE, -df_error) |> 
  arrange(var_dep_name, var_indep_name)

res_mixMod_covar_02_gk <- res_mixMod_covar_01 |> 
  select(starts_with("var_"), ends_with("_gk")) |> 
  left_join(res_mixMod_covar_01_tab |> 
              select(starts_with("var_"), gk_p_val_comp)) |> 
  group_by(var_dep_name) |> 
  mutate(
    covariate = "glucocoticoids",
    n_sig = sum(gk_p_val_comp < 0.06, na.rm = TRUE),
    mod_final = if_else(n_sig > 2, mod_covar_gk, mod_NOcovar_gk),
    adjustment  = if_else(n_sig > 2, "adjusted", "without"),
    tidier = map(mod_final, ~model_parameters(
      .x$fit,
      ci = 0.95,
      df_method = "satterthwaite"))
  ) |>
  ungroup() |>
  unnest(tidier) |> 
  filter(Parameter == "var_indep_value") |> 
  select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
         -CI, -Effects, -SE, -df_error) |> 
  arrange(var_dep_name, var_indep_name)

res_mixMod_covar_02_ck <- res_mixMod_covar_01 |> 
  select(starts_with("var_"), ends_with("_ck")) |> 
  left_join(res_mixMod_covar_01_tab |> 
              select(starts_with("var_"), ck_p_val_comp)) |> 
  group_by(var_dep_name) |> 
  mutate(
    covariate = "creatinkinase",
    n_sig = sum(ck_p_val_comp < 0.06, na.rm = TRUE),
    mod_final = if_else(n_sig > 2, mod_covar_ck, mod_NOcovar_ck),
    adjustment  = if_else(n_sig > 2, "adjusted", "without"),
    tidier = map(mod_final, ~model_parameters(
      .x$fit,
      ci = 0.95,
      df_method = "satterthwaite"))
  ) |>
  ungroup() |>
  unnest(tidier) |> 
  filter(Parameter == "var_indep_value") |> 
  select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
         -CI, -Effects, -SE, -df_error) |> 
  arrange(var_dep_name, var_indep_name)


res_mixMod_covar_02_bmc <- res_mixMod_covar_01 |> 
  select(starts_with("var_"), ends_with("_bcm")) |> 
  left_join(res_mixMod_covar_01_tab |> 
              select(starts_with("var_"), bcm_p_val_comp)) |> 
  group_by(var_dep_name) |> 
  mutate(
    covariate = "body cell mass",
    n_sig = sum(bcm_p_val_comp < 0.06, na.rm = TRUE),
    mod_final = if_else(n_sig > 2, mod_covar_bcm, mod_NOcovar_bcm),
    adjustment  = if_else(n_sig > 2, "adjusted", "without"),
    tidier = map(mod_final, ~model_parameters(
      .x$fit,
      ci = 0.95,
      df_method = "satterthwaite"))
  ) |>
  ungroup() |>
  unnest(tidier) |> 
  filter(Parameter == "var_indep_value") |> 
  select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
         -CI, -Effects, -SE, -df_error) |> 
  arrange(var_dep_name, var_indep_name)

res_mixMod_covar_02_sex <- res_mixMod_covar_01 |> 
  select(starts_with("var_"), ends_with("_sex")) |> 
  left_join(res_mixMod_covar_01_tab |> 
              select(starts_with("var_"), sex_p_val_comp)) |> 
  group_by(var_dep_name) |> 
  mutate(
    covariate = "gender",
    n_sig = sum(sex_p_val_comp < 0.06, na.rm = TRUE),
    mod_final = if_else(n_sig > 2, mod_covar_sex, mod_NOcovar_sex),
    adjustment  = if_else(n_sig > 2, "adjusted", "without"),
    tidier = map(mod_final, ~model_parameters(
      .x$fit,
      ci = 0.95,
      df_method = "satterthwaite"))
  ) |>
  ungroup() |>
  unnest(tidier) |> 
  filter(Parameter == "var_indep_value") |> 
  select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
         -CI, -Effects, -SE, -df_error) |> 
  arrange(var_dep_name, var_indep_name)

res_mixMod_covar_02_sub <- res_mixMod_covar_01 |> 
  select(starts_with("var_"), ends_with("_sub")) |> 
  left_join(res_mixMod_covar_01_tab |> 
              select(starts_with("var_"), sub_p_val_comp)) |> 
  group_by(var_dep_name) |> 
  mutate(
    covariate = "subdiagnosis",
    n_sig = sum(sub_p_val_comp < 0.06, na.rm = TRUE),
    mod_final = if_else(n_sig > 2, mod_covar_sub, mod_NOcovar_sub),
    adjustment  = if_else(n_sig > 2, "adjusted", "without"),
    tidier = map(mod_final, ~model_parameters(
      .x$fit,
      ci = 0.95,
      df_method = "satterthwaite"))
  ) |>
  ungroup() |>
  unnest(tidier) |> 
  filter(Parameter == "var_indep_value") |> 
  select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
         -CI, -Effects, -SE, -df_error) |> 
  arrange(var_dep_name, var_indep_name)

res_mixMod_covar_02_jo <- res_mixMod_covar_01 |> 
  select(starts_with("var_"), ends_with("_jo")) |> 
  left_join(res_mixMod_covar_01_tab |> 
              select(starts_with("var_"), jo_p_val_comp)) |> 
  group_by(var_dep_name) |> 
  mutate(
    covariate = "anti-jo-1",
    n_sig = sum(jo_p_val_comp < 0.06, na.rm = TRUE),
    mod_final = if_else(n_sig > 2, mod_covar_jo, mod_NOcovar_jo),
    adjustment  = if_else(n_sig > 2, "adjusted", "without"),
    tidier = map(mod_final, ~model_parameters(
      .x$fit,
      ci = 0.95,
      df_method = "satterthwaite"))
  ) |>
  ungroup() |>
  unnest(tidier) |> 
  filter(Parameter == "var_indep_value") |> 
  select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
         -CI, -Effects, -SE, -df_error) |> 
  arrange(var_dep_name, var_indep_name)


res_mixMod_covar_02_hmgcr <- res_mixMod_covar_01 |> 
  select(starts_with("var_"), ends_with("_hmgcr")) |> 
  left_join(res_mixMod_covar_01_tab |> 
              select(starts_with("var_"), hmgcr_p_val_comp)) |> 
  group_by(var_dep_name) |> 
  mutate(
    covariate = "anti-hmgcr",
    n_sig = sum(hmgcr_p_val_comp < 0.06, na.rm = TRUE),
    mod_final = if_else(n_sig > 2, mod_covar_hmgcr, mod_NOcovar_hmgcr),
    adjustment  = if_else(n_sig > 2, "adjusted", "without"),
    tidier = map(mod_final, ~model_parameters(
      .x$fit,
      ci = 0.95,
      df_method = "satterthwaite"))
  ) |>
  ungroup() |>
  unnest(tidier) |> 
  filter(Parameter == "var_indep_value") |> 
  select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
         -CI, -Effects, -SE, -df_error) |> 
  arrange(var_dep_name, var_indep_name)

res_mixMod_covar_02_int_sex <- res_mixMod_covar_01 |> 
  select(starts_with("var_"), ends_with("_int_sex")) |> 
  left_join(res_mixMod_covar_01_tab |> 
              select(starts_with("var_"), int_sex_p_val_comp)) |> 
  group_by(var_dep_name) |> 
  mutate(
    covariate = "sex-interaction",
    n_sig = sum(int_sex_p_val_comp < 0.06, na.rm = TRUE),
    mod_final = if_else(n_sig > 2, mod_covar_int_sex, mod_NOcovar_int_sex),
    adjustment  = if_else(n_sig > 2, "adjusted", "without"),
    tidier = map(mod_final, ~model_parameters(
      .x$fit,
      ci = 0.95,
      df_method = "satterthwaite"))
  ) |>
  ungroup() |>
  unnest(tidier) |> 
  filter(Parameter == "var_indep_value") |> 
  select(-n_sig, -starts_with("mod_"), -Parameter, - Group, -t, 
         -CI, -Effects, -SE, -df_error) |> 
  arrange(var_dep_name, var_indep_name)

### export ----
export(list(overall = res_mixMod_covar_01_tab,
            age = res_mixMod_covar_02_vek,
            gk = res_mixMod_covar_02_gk,
            creatine = res_mixMod_covar_02_ck,
            bmc = res_mixMod_covar_02_bmc,
            gender = res_mixMod_covar_02_sex,
            subdiagnosis = res_mixMod_covar_02_sub,
            anti_jo_1 = res_mixMod_covar_02_jo,
            anti_hmgcr = res_mixMod_covar_02_hmgcr,
            sex_interaction = res_mixMod_covar_02_int_sex),
       "output/tables/250922_covariates_01.xlsx")

## covariate models 2 ----
### only covariates ----
#### models ----
tic()
Sys.time()

cluster <- multidplyr::new_cluster(ncores)
multidplyr::cluster_library(cluster, c("dplyr","purrr","stringr","tibble", 
                                       "data.table", "tidyverse", "tidymodels",
                                       "multilevelmod","lme4","lmerTest",
                                       "broom.mixed", "easystats"))
multidplyr::cluster_copy(cluster, c("model_covar","model_min", "model_comp_int"))

res_mixMod_covar_02 <- d04_sel2 |>
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
  multidplyr::partition(cluster) |> 
  do(
    data = data.table(.),
    mod_without = model_min(., filter_na = FALSE),
    mod_covar_gk = model_covar(., "denna_davka_gk_mg_kg_bw", filter_na = FALSE),
    mod_covar_age = model_covar(., "vek", filter_na = FALSE),
    mod_covar_bcm = model_covar(., "bcm", filter_na = FALSE),
    mod_covar_creatine = model_covar(., "kreatinin_umol_l", filter_na = FALSE)
  ) |> 
  collect()

walk(cluster, ~ .x$kill()); gc()

toc()
beep(4)

#### tables ----
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
# final table
res_mixMod_covar_02_tab_confounders <- res_mixMod_covar_02_tab |> 
  full_join(res_mixMod_covar_02_tab_age, by = c("var_dep_name", "var_indep_name")) |> 
  full_join(res_mixMod_covar_02_tab_gk, by = c("var_dep_name", "var_indep_name")) |> 
  full_join(res_mixMod_covar_02_tab_bcm, by = c("var_dep_name", "var_indep_name")) |> 
  full_join(res_mixMod_covar_02_tab_creatine, by = c("var_dep_name", "var_indep_name"))

#### export ----
export(res_mixMod_covar_02_tab_confounders, "output/tables/250930_covariates_01.xlsx")

### subgrouping ----
#### age ----
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

#### gk ----
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

#### bcm ----
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

#### creatine ----
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

#### gk + bcm ----
res_mixMod_covar_02_sub_gk_bcm <- res_mixMod_covar_02 |> 
  select(var_dep_name, var_indep_name, data) |> 
  inner_join(res_mixMod_covar_02_tab_creatine |> 
               filter(if_any(last_col(), ~ . < 0.06)) |> 
               select(var_dep_name, var_indep_name),
             by = c("var_dep_name", "var_indep_name")) |> 
  mutate(fig_sex = pmap(list(data = data, 
                             ylab = var_dep_name, 
                             xlab = var_indep_name,
                             covar = "denna_davka_gk_mg_kg_bw + bcm",
                             covar_spec = "pohlavi"),
                        plot_lmer_grid_basic_02),
         fig_sub = pmap(list(data = data, 
                             ylab = var_dep_name, 
                             xlab = var_indep_name,
                             covar = "denna_davka_gk_mg_kg_bw + bcm",
                             covar_spec = "podtyp_nemoci_zjednoduseny"),
                        plot_lmer_grid_basic_02),
         fig_jo = pmap(list(data = data, 
                            ylab = var_dep_name, 
                            xlab = var_indep_name,
                            covar = "denna_davka_gk_mg_kg_bw + bcm",
                            covar_spec = "jo_1"),
                       plot_lmer_grid_basic_02),
         fig_hmgcr = pmap(list(data = data, 
                               ylab = var_dep_name, 
                               xlab = var_indep_name,
                               covar = "denna_davka_gk_mg_kg_bw + bcm",
                               covar_spec = "anti_hmgcr"),
                          plot_lmer_grid_basic_02)
  )



# Cluster analyses ----
## M0 ----
### FAMD ----
# https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/115-famd-factor-analysis-of-mixed-data-in-r-essentials
res_famd_m0 <- FAMD(
  # d08_m0 |> select( -odpoved_na_terapii_m0_vs_m6),
  d08_m0 |>filter(!is.na(odpoved_na_terapii_m0_vs_m6)) |> 
    select(-ends_with("_total")), # do tabulky byly dodány mda a mmt8 kategorie
  graph = TRUE
)
# figures
tiff("output/figures/cluster analyses/251113_m0_fmda_overview_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_famd_var(res_famd_m0,labelsize = 8,repel = TRUE) +
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))
dev.off()
## feature contribution
tiff("output/figures/cluster analyses/251113_m0_fmda_feat_contrib_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_contrib(res_famd_m0, "var", axes = 1,labelsize = 8) +
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, angle = 45, hjust = 1),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14))
dev.off()
## quantitative variables - continuous
tiff("output/figures/cluster analyses/251113_m0_fmda_continuous_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_famd_var(res_famd_m0, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              labelsize = 8,
              repel = TRUE) +
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))
dev.off()
## qualitative variables - factors
tiff("output/figures/cluster analyses/251113_m0_fmda_factor_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_famd_var(res_famd_m0, "quali.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              labelsize = 8,
              repel = TRUE) +
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))
dev.off()
  
## podtyp_nemoci_zjednoduseny
tiff("output/figures/cluster analyses/251113_m0_fmda_subdisease_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_mfa_ind(
  res_famd_m0,
  habillage = "podtyp_nemoci_zjednoduseny",   # color points by your factor
  palette = c("#00AFBB", "#E7B800", "#FC4E07", "#8B4513"),  # optional extra color for 4 levels
  addEllipses = TRUE,
  ellipse.type = "confidence",
  repel = FALSE,   # do not show text labels
  geom = "point"   # only draw points, no text
)+
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))
dev.off()
  
## odpoved_na_terapii_m0_vs_m6
tiff("output/figures/cluster analyses/251113_m0_fmda_respose_trt_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_mfa_ind(
  res_famd_m0,
  habillage = "odpoved_na_terapii_m0_vs_m6",   # color points by your factor
  palette = c("#00AFBB", "#E7B800", "#FC4E07", "#8B4513"),  # optional extra color for 4 levels
  addEllipses = TRUE,
  ellipse.type = "confidence",
  repel = FALSE,   # do not show text labels
  geom = "point"   # only draw points, no text
) +
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))
dev.off()

## mda_categories
tiff("output/figures/cluster analyses/251113_m0_fmda_mda_cat_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_mfa_ind(
  res_famd_m0,
  habillage = "mda_categories",   # color points by your factor
  palette = c("#00AFBB", "#E7B800", "#FC4E07", "#8B4513"),  # optional extra color for 4 levels
  addEllipses = TRUE,
  ellipse.type = "confidence",
  repel = FALSE,   # do not show text labels
  geom = "point"   # only draw points, no text
) +
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))
dev.off()

## mmt8_categories
tiff("output/figures/cluster analyses/251113_m0_fmda_mmt8_cat_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_mfa_ind(
  res_famd_m0,
  habillage = "mmt8_categories",   # color points by your factor
  palette = c("#00AFBB", "#E7B800", "#FC4E07", "#8B4513"),  # optional extra color for 4 levels
  addEllipses = TRUE,
  ellipse.type = "confidence",
  repel = FALSE,   # do not show text labels
  geom = "point"   # only draw points, no text
) +
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))
dev.off()

### UMAP ----
famd_scores_m0 <- res_famd_m0$ind$coord   # numeric low-dim representation

umap_res_m0 <- umap(famd_scores_m0, n_neighbors = 15, min_dist = 0.1)


# tiff("output/figures/cluster analyses/251113_m0_umap_01.tiff", 
#      compression = "lzw",
#      res = "300",
#      width = 1600,
#      height = 1200)
# plot(umap_res_m0, col = as.numeric(d08_m0$podtyp_nemoci_zjednoduseny), pch = 19)
# dev.off()

d08_m0_2 <- d08_m0 |>filter(!is.na(odpoved_na_terapii_m0_vs_m6))

umap_df <- as.data.frame(umap_res_m0) |>
  setNames(c("UMAP1", "UMAP2")) |>
  mutate(
    podtyp = d08_m0_2$podtyp_nemoci_zjednoduseny,
    mmt8 = d08_m0_2$mmt8_categories,
    mda = d08_m0_2$mda_categories,
    odpoved = d08_m0_2$odpoved_na_terapii_m0_vs_m6
  )


vars_to_plot <- c("podtyp", "mmt8", "mda", "odpoved")

walk(vars_to_plot, function(v) {
  
  p <- plot_umap(umap_df, v)
  
  tiff(
    filename = paste0("output/figures/cluster analyses/251113_umap_", v, ".tiff"),
    compression = "lzw",
    res = 300,
    width = 1600,
    height = 1200
  )
  print(p)
  dev.off()
})

### PLS-DA ----
# # disease subtype
# d08_m0_pls_X <- model.matrix(
#   ~ . - 1,
#   data = d08_m0_pls |> 
#     filter(!is.na(odpoved_na_terapii_m0_vs_m6)) |> 
#     select(-podtyp_nemoci_zjednoduseny)) 
# d08_m0_pls_Y <- d08_m0_pls |> 
#   filter(!is.na(odpoved_na_terapii_m0_vs_m6)) |> 
#   select(podtyp_nemoci_zjednoduseny) |> 
#   pull()
# 
# pls_m0 <- plsda(d08_m0_pls_X, d08_m0_pls_Y, 
#                 ncomp = 2, scale = TRUE)
# 
# ## visualization
# tiff("output/figures/cluster analyses/251113_m0_pls_subdisease_group_01.tiff", 
#      compression = "lzw",
#      res = "300",
#      width = 1600,
#      height = 1200)
# plotIndiv(pls_m0, comp = c(1,2), 
#           group = d08_m0_pls_Y,ellipse = TRUE, 
#           legend=T) 
# dev.off()
# tiff("output/figures/cluster analyses/251113_m0_pls_subdisease_contrib_01.tiff", 
#      compression = "lzw",
#      res = "300",
#      width = 1600,
#      height = 1200)
# plotLoadings(pls_m0, contrib = 'max', method = 'median')
# dev.off()
# 
# # response to treatmetn
# d08_m0_pls_X_trt <- model.matrix(
#   ~ . - 1,
#   data = d08_m0_pls |> 
#     filter(!is.na(odpoved_na_terapii_m0_vs_m6)) |> 
#     select(-odpoved_na_terapii_m0_vs_m6)) 
# d08_m0_pls_Y_trt <- d08_m0_pls |> 
#   filter(!is.na(odpoved_na_terapii_m0_vs_m6)) |> 
#   select(odpoved_na_terapii_m0_vs_m6) |> 
#   pull()
# 
# pls_m0_trt <- plsda(d08_m0_pls_X_trt, d08_m0_pls_Y_trt, 
#                     ncomp = 2, scale = TRUE)
# 
# ## visualization
# tiff("output/figures/cluster analyses/251113_m0_pls_response_group_01.tiff", 
#      compression = "lzw",
#      res = "300",
#      width = 1600,
#      height = 1200)
# plotIndiv(pls_m0_trt, comp = c(1,2), group = d08_m0_pls_Y_trt, ellipse = TRUE, legend=T)
# dev.off()
# tiff("output/figures/cluster analyses/251113_m0_pls_response_contrib_01.tiff", 
#      compression = "lzw",
#      res = "300",
#      width = 1600,
#      height = 1200)
# plotLoadings(pls_m0_trt, contrib = 'max', method = 'median')
# dev.off()

source("scripts/cluster_plsda_pipeline.R")

pls_results <- run_all_plsda(
  df = d08_m0_pls,
  prefix = "251113_m0_pls"
)

## delta ----
### FAMD ----
# https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/115-famd-factor-analysis-of-mixed-data-in-r-essentials
res_famd_delta <- FAMD(
  # d08_delta |> select( -odpoved_na_terapii_m0_vs_m6),
  d08_delta |>filter(!is.na(odpoved_na_terapii_m0_vs_m6)) |> 
    select(-contains("_total")), # do tabulky byly dodány mda a mmt8 kategorie
  graph = TRUE
)
# figures
tiff("output/figures/cluster analyses/251113_delta_fmda_overview_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_famd_var(res_famd_delta,labelsize = 8,repel = TRUE) +
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))
dev.off()
## feature contribution
tiff("output/figures/cluster analyses/251113_delta_fmda_feat_contrib_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_contrib(res_famd_delta, "var", axes = 1,labelsize = 8) +
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, angle = 45, hjust = 1),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14))
dev.off()
## quantitative variables - continuous
tiff("output/figures/cluster analyses/251113_delta_fmda_continuous_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_famd_var(res_famd_delta, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              labelsize = 8,
              repel = TRUE) +
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))
dev.off()
## qualitative variables - factors
tiff("output/figures/cluster analyses/251113_delta_fmda_factor_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_famd_var(res_famd_delta, "quali.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              labelsize = 8,
              repel = TRUE) +
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))
dev.off()

## podtyp_nemoci_zjednoduseny
tiff("output/figures/cluster analyses/251113_delta_fmda_subdisease_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_mfa_ind(
  res_famd_delta,
  habillage = "podtyp_nemoci_zjednoduseny",   # color points by your factor
  palette = c("#00AFBB", "#E7B800", "#FC4E07", "#8B4513"),  # optional extra color for 4 levels
  addEllipses = TRUE,
  ellipse.type = "confidence",
  repel = FALSE,   # do not show text labels
  geom = "point"   # only draw points, no text
)+
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))
dev.off()

## odpoved_na_terapii_m0_vs_m6
tiff("output/figures/cluster analyses/251113_delta_fmda_respose_trt_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_mfa_ind(
  res_famd_delta,
  habillage = "odpoved_na_terapii_m0_vs_m6",   # color points by your factor
  palette = c("#00AFBB", "#E7B800", "#FC4E07", "#8B4513"),  # optional extra color for 4 levels
  addEllipses = TRUE,
  ellipse.type = "confidence",
  repel = FALSE,   # do not show text labels
  geom = "point"   # only draw points, no text
) +
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))
dev.off()

## mda_categories
tiff("output/figures/cluster analyses/251113_delta_fmda_mda_cat_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_mfa_ind(
  res_famd_delta,
  habillage = "mda_categories",   # color points by your factor
  palette = c("#00AFBB", "#E7B800", "#FC4E07", "#8B4513"),  # optional extra color for 4 levels
  addEllipses = TRUE,
  ellipse.type = "confidence",
  repel = FALSE,   # do not show text labels
  geom = "point"   # only draw points, no text
) +
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))
dev.off()

## mmt8_categories
tiff("output/figures/cluster analyses/251113_delta_fmda_mmt8_cat_01.tiff", 
     compression = "lzw",
     res = "300",
     width = 1600,
     height = 1200)
fviz_mfa_ind(
  res_famd_delta,
  habillage = "mmt8_categories",   # color points by your factor
  palette = c("#00AFBB", "#E7B800", "#FC4E07", "#8B4513"),  # optional extra color for 4 levels
  addEllipses = TRUE,
  ellipse.type = "confidence",
  repel = FALSE,   # do not show text labels
  geom = "point"   # only draw points, no text
) +
  theme_minimal(base_size = 16) +   # increases all text size (titles, axis, legend)
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))
dev.off()

### UMAP ----
famd_scores_delta <- res_famd_delta$ind$coord   # numeric low-dim representation
umap_res_delta <- umap(famd_scores_delta, n_neighbors = 15, min_dist = 0.1)

# tiff("output/figures/cluster analyses/251113_delta_umap_01.tiff", 
#      compression = "lzw",
#      res = "300",
#      width = 1600,
#      height = 1200)
# plot(umap_res_delta, col = as.numeric(d08_delta$podtyp_nemoci_zjednoduseny), pch = 19)
# dev.off()

d08_delta_2 <- d08_delta |>filter(!is.na(odpoved_na_terapii_m0_vs_m6))

umap_df <- as.data.frame(umap_res_delta) |>
  setNames(c("UMAP1", "UMAP2")) |>
  mutate(
    podtyp = d08_delta_2$podtyp_nemoci_zjednoduseny,
    mmt8 = d08_delta_2$mmt8_categories,
    mda = d08_delta_2$mda_categories,
    odpoved = d08_delta_2$odpoved_na_terapii_m0_vs_m6
  )


vars_to_plot <- c("podtyp", "mmt8", "mda", "odpoved")

walk(vars_to_plot, function(v) {
  
  p <- plot_umap(umap_df, v)
  
  tiff(
    filename = paste0("output/figures/cluster analyses/251113_umap_delta_", v, ".tiff"),
    compression = "lzw",
    res = 300,
    width = 1600,
    height = 1200
  )
  print(p)
  dev.off()
})

### PLS-DA ----
# # disease subtype
# d08_delta_pls_X <- model.matrix(
#   ~ . - 1,
#   data = d08_delta_pls |> 
#     filter(!is.na(odpoved_na_terapii_m0_vs_m6)) |> 
#     select(-podtyp_nemoci_zjednoduseny)) 
# d08_delta_pls_Y <- d08_delta_pls |> 
#   filter(!is.na(odpoved_na_terapii_m0_vs_m6)) |> 
#   select(podtyp_nemoci_zjednoduseny) |> 
#   pull()
# 
# pls_delta <- plsda(d08_delta_pls_X, d08_delta_pls_Y, 
#                    ncomp = 2, scale = TRUE)
# 
# ## visualization
# tiff("output/figures/cluster analyses/251113_delta_pls_subdisease_group_01.tiff", 
#      compression = "lzw",
#      res = "300",
#      width = 1600,
#      height = 1200)
# plotIndiv(pls_delta, comp = c(1,2), 
#           group = d08_delta_pls_Y,ellipse = TRUE, 
#           legend=T) 
# dev.off()
# tiff("output/figures/cluster analyses/251113_delta_pls_subdisease_contrib_01.tiff", 
#      compression = "lzw",
#      res = "300",
#      width = 1600,
#      height = 1200)
# plotLoadings(pls_delta, contrib = 'max', method = 'median')
# dev.off()
# 
# # response to treatmetn
# d08_delta_pls_X_trt <- model.matrix(
#   ~ . - 1,
#   data = d08_delta_pls |> 
#     filter(!is.na(odpoved_na_terapii_m0_vs_m6)) |> 
#     select(-odpoved_na_terapii_m0_vs_m6)) 
# d08_delta_pls_Y_trt <- d08_delta_pls |> 
#   filter(!is.na(odpoved_na_terapii_m0_vs_m6)) |> 
#   select(odpoved_na_terapii_m0_vs_m6) |> 
#   pull()
# 
# pls_delta_trt <- plsda(d08_delta_pls_X_trt, d08_delta_pls_Y_trt, 
#                        ncomp = 2, scale = TRUE)
# 
# ## visualization
# tiff("output/figures/cluster analyses/251113_delta_pls_response_group_01.tiff", 
#      compression = "lzw",
#      res = "300",
#      width = 1600,
#      height = 1200)
# plotIndiv(pls_delta_trt, comp = c(1,2), group = d08_delta_pls_Y_trt, ellipse = TRUE, legend=T)
# dev.off()
# tiff("output/figures/cluster analyses/251113_delta_pls_response_contrib_01.tiff", 
#      compression = "lzw",
#      res = "300",
#      width = 1600,
#      height = 1200)
# plotLoadings(pls_delta_trt, contrib = 'max', method = 'median')
# dev.off()

source("scripts/cluster_plsda_pipeline.R")

pls_delta_results <- run_all_plsda(
  df = d08_delta_pls,
  prefix = "251113_delta_pls"
)



