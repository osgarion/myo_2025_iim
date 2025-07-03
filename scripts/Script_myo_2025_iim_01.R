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

### selected variables ----
# selection ----
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
lower = 0.6; upper <- 0.7; examination == "M0"
cor_data <- d03 |> 
  filter(poradie_vysetrenia == examination) |> 
  select(where(is.numeric)) |> 
  select(where(~mean(is.na(.x)) < 0.8))

cor_mat <- cor(cor_data, method = "kendall", use = "pairwise.complete.obs")
cor_mat[is.na(cor_mat)] <- 0 # nahradí všechny NA za 0

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

ggpairs(data |> 
          filter(poradie_vysetrenia == examination),
        columns = data |> select(all_of(param_name_sig)) |> names(), 
        lower = list(continuous = my_fn, method = "kendall")) +
  theme_sjplot2()

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
fviz_eig(pca_m0, col)
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


# Statistics ----
# 
# 


