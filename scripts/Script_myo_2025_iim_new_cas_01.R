# Script_myo_2025_iim_new_cas_01.R
# Auto-generated: d10_cas -> d10_cas_b, date-stamped outputs via out_date
out_date <- format(Sys.Date(), "%y%m%d")

local({

# 1. Libraries & functions ----
# *******************************************************
# required
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, purrr, conflicted, renv, rio)
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(recipes::update)
}
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}
# Functions and libraries uploading
list.files("scripts/functions/", pattern="*.*", full.names=TRUE) |>
  map(~source(.))

# -------------------------------------------------------
# d13_join/d13_mofa override: use d10_cas_b
# -------------------------------------------------------
d13_join <- d10_cas_b |>
  transmute(
    immet_id   = as.character(sample),
    exam_order = as.character(exam),
    across(where(is.numeric), identity)
  ) |>
  inner_join(
    d11_clinics_imp |>
      transmute(
        immet_id   = as.character(immet_id),
        exam_order = as.character(exam_order),
        across(-c(immet_id, exam_order), identity)
      ),
    by = c("immet_id", "exam_order")
  ) |>
  distinct() |>
  na.omit() |> 
  group_by(immet_id) |>
  filter(n() == 2) |>
  ungroup()

d13_mofa <- d13_join |>
  mutate(
    sample  = paste0(immet_id, "_", exam_order),
    subtype = as.character(disease_subtype) # convenience col for color_by
  ) |>
  distinct(sample, .keep_all = TRUE)
# -------------------------------------------------------
renv::status()

# 2. Environment reproducibility ----
# *******************************************************
# RENV
# renv::init()        # once
# renv::status()      # sync?
# renv::snapshot()    # after changing deps
# renv::restore()     # on a new machine

# Back-up
back_up("scripts/functions/FUN_01.R") # the the destination subdirectory specify using 'path_dest'
back_up("scripts/functions/OBJ_01.R") # the the destination subdirectory specify using 'path_dest'
back_up("scripts/Script_myo_2025_iim_02_ggir.R") # the the destination subdirectory specify using 'path_dest'

# Save to .RData
save(xx1,
     file = "reports/markD_01.RData")
save.image("reports/markD_01.RData")
save.image(file.path("data", paste0(out_date, "_backup.RData")))
# Save the tables into data/tables.RData by listing them individually
cgwtools::resave(xx2,xx3, file = "reports/markD_01.RData") # resave a list of tables that I'll use in the .Rmd file.
# Save the tables into data/tables.RData using "patterns" 
cgwtools::resave(list=ls(pattern="tbl"), file = "reports/markD_01.RData")

# *******************************************************
# 3. EDA ----
# *******************************************************

## figures ----
### ostatni ----
d10_ostatni |>
  ggplot(aes(ig_intercept_enmo, ig_gradient_enmo)) +
  geom_point(shape = 15, col = "darkblue") +
  geom_smooth(
    method = "loess", se = FALSE,
    linetype = "dashed", col = "darkblue"
  ) +
  geom_point(
    data = d10_ostatni |> filter(r2_enmo > 0.9),
    shape = 17, col = "darkgreen"
  ) +
  geom_smooth(
    data = d10_ostatni |> filter(r2_enmo > 0.9),
    method = "loess", se = FALSE,
    linetype = "dashed", col = "darkgreen"
  ) +
  labs(
    title = "Intercept vs. gradient (ENMO 0–24 h)",
    subtitle = "Dashed LOESS fit for all samples (blue) and high-fit subset R² > 0.9 (green)",
    x = "Intercept (ENMO, 0–24 h)",
    y = "Gradient (ENMO, 0–24 h)"
  ) +
  sjPlot::theme_sjplot2() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(size = 15, face = "italic"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, face = "bold")
  )

## q-q plots ----
# *******************************************************
### čas ----
# *******************************************************
d10_cas_b |> 
  filter(exam %in% c("M0", "M6")) |>
  pivot_longer(cols = where(is.numeric),
               names_to = "var",
               values_to = "values") |> 
  ggqqplot(x = "values",
           color = "exam",
           palette = c("#00AFBB", "#E7B800"),
           facet.by = "var", scales = "free")





## ggally ----
# *******************************************************
### čas ----
# *******************************************************
d10_cas_b |> 
  filter(exam %in% c("M0", "M6")) |>
  (\(df)
   ggpairs(
     df,
     columns = 3:ncol(df),
     mapping = ggplot2::aes(color = exam),
     lower = list(continuous = my_fn),
     upper = list(continuous = GGally::wrap("cor", method = "pearson"))
   )
  )()

### aktivita ----
# *******************************************************
d10_aktivita |> 
  filter(exam %in% c("M0", "M6")) |>
  (\(df)
   ggpairs(
     df,
     columns = 3:ncol(df),
     mapping = ggplot2::aes(color = exam),
     lower = list(continuous = my_fn),
     upper = list(continuous = GGally::wrap("cor", method = "pearson"))
   )
  )()

# *******************************************************
## missingness ----
# *******************************************************
vis_dat(d11_clinics)
vis_miss(d11_clinics |> 
           filter(exam_order %in% c("M0", "M6")) |> 
           select(
             where(~ mean(is.na(.x)) <= 0.20)
           ))

# **********************************************************************
## correlation ----
# **********************************************************************
### AKTIVITA ----
# **********************************************************************
#### Data preparation ----
# **********************************************************************
x_cols_corr_activ <- d10_aktivita |>
  dplyr::select(where(is.numeric)) |>
  names()

y_cols_corr_activ <- d11_clinics_imp |>
  dplyr::select(where(is.numeric)) |>
  names()

d12_join_all_activ <- d12_join

d12_join_m0_activ <- d12_join |>
  dplyr::filter(tolower(as.character(exam_order)) == "m0")

d12_join_m6_activ <- d12_join |>
  dplyr::filter(tolower(as.character(exam_order)) == "m6")

delta_df_activ <- d12_join |>
  dplyr::mutate(exam_order = tolower(as.character(exam_order))) |>
  dplyr::filter(exam_order %in% c("m0", "m6")) |>
  dplyr::select(immet_id, exam_order, all_of(x_cols_corr_activ), all_of(y_cols_corr_activ)) |>
  tidyr::pivot_wider(
    names_from  = exam_order,
    values_from = c(all_of(x_cols_corr_activ), all_of(y_cols_corr_activ))
  )

# keep only complete pairs (both m0 and m6 exist) - "soft" filter; NAs handled later via pairwise.complete.obs
delta_df_activ <- delta_df_activ |>
  dplyr::filter(
    if_all(dplyr::all_of(c(paste0(x_cols_corr_activ, "_m0"), paste0(x_cols_corr_activ, "_m6"))), ~ !is.na(.x)) |
      if_all(dplyr::all_of(c(paste0(y_cols_corr_activ, "_m0"), paste0(y_cols_corr_activ, "_m6"))), ~ !is.na(.x))
  )

# compute deltas into a new compact table with SAME column names as x_cols/y_cols
delta_X_activ <- purrr::map_dfc(
  x_cols_corr_activ,
  \(nm) delta_df_activ[[paste0(nm, "_m6")]] - delta_df_activ[[paste0(nm, "_m0")]]
) |>
  stats::setNames(x_cols_corr_activ)

delta_Y_activ <- purrr::map_dfc(
  y_cols_corr_activ,
  \(nm) delta_df_activ[[paste0(nm, "_m6")]] - delta_df_activ[[paste0(nm, "_m0")]]
) |>
  stats::setNames(y_cols_corr_activ)

d12_join_delta_activ <- dplyr::bind_cols(
  tibble::tibble(immet_id = delta_df_activ$immet_id),
  delta_X_activ,
  delta_Y_activ
) |> select(-age) |> 
  left_join(d12_join |>
              dplyr::filter(tolower(as.character(exam_order)) == "m0") |>
              dplyr::select(immet_id, age) |>
              dplyr::distinct())

# **********************************************************************
#### Analyses ----
# **********************************************************************
corr_res_all_activ <- make_corr_heatmap_kendall(
  d_join = d12_join_all_activ,
  x_cols = x_cols_corr_activ,
  y_cols = y_cols_corr_activ,
  title_suffix = "AKTIVITA: ALL"
)
cor_df_all_activ  <- corr_res_all_activ$cor_df
(fig_cor_all_activ <- corr_res_all_activ$fig)

corr_res_m0_activ <- make_corr_heatmap_kendall(
  d_join = d12_join_m0_activ,
  x_cols = x_cols_corr_activ,
  y_cols = y_cols_corr_activ,
  title_suffix = "AKTIVITA: M0"
)
cor_df_m0_activ  <- corr_res_m0_activ$cor_df
(fig_cor_m0_activ <- corr_res_m0_activ$fig)

corr_res_m6_activ <- make_corr_heatmap_kendall(
  d_join = d12_join_m6_activ,
  x_cols = x_cols_corr_activ,
  y_cols = y_cols_corr_activ,
  title_suffix = "AKTIVITA: M6"
)
cor_df_m6_activ  <- corr_res_m6_activ$cor_df
(fig_cor_m6_activ <- corr_res_m6_activ$fig)

corr_res_delta_activ <- make_corr_heatmap_kendall(
  d_join = d12_join_delta_activ,
  x_cols = x_cols_corr_activ,
  y_cols = y_cols_corr_activ,
  title_suffix = "AKTIVITA: DELTA (M6 - M0)"
)
cor_df_delta_activ  <- corr_res_delta_activ$cor_df
(fig_cor_delta_activ <- corr_res_delta_activ$fig)

# **********************************************************************
#### Export ----
# **********************************************************************
# figures
tiff(
  filename = file.path("output/figures/final", paste0(format(Sys.Date(), "%y%m%d"), "_corr_kendall_activ_all.tiff")),
  compression = "lzw", res = "150", width = 1200, height = 1200
)
print(fig_cor_all_activ)
dev.off()

tiff(
  filename = file.path("output/figures/final", paste0(format(Sys.Date(), "%y%m%d"), "_corr_kendall_activ_m0.tiff")),
  compression = "lzw", res = "150", width = 1200, height = 1200
)
print(fig_cor_m0_activ)
dev.off()

tiff(
  filename = file.path("output/figures/final", paste0(format(Sys.Date(), "%y%m%d"), "_corr_kendall_activ_m6.tiff")),
  compression = "lzw", res = "150", width = 1200, height = 1200
)
print(fig_cor_m6_activ)
dev.off()

tiff(
  filename = file.path("output/figures/final", paste0(format(Sys.Date(), "%y%m%d"), "_corr_kendall_activ_delta_m6_m0.tiff")),
  compression = "lzw", res = "150", width = 1200, height = 1200
)
print(fig_cor_delta_activ)
dev.off()

# # tables (xlsx)
# corr_table_all_activ <- format_corr_table(cor_df_all_activ)
# rio::export(
#   corr_table_all_activ,
#   file.path("output/tables/final", paste0(format(Sys.Date(), "%y%m%d"), "_kendall_aktivita_vs_clinical_ALL.xlsx"))
# )
# 
# corr_table_m0_activ <- format_corr_table(cor_df_m0_activ)
# rio::export(
#   corr_table_m0_activ,
#   file.path("output/tables/final", paste0(format(Sys.Date(), "%y%m%d"), "_kendall_aktivita_vs_clinical_M0.xlsx"))
# )
# 
# corr_table_m6_activ <- format_corr_table(cor_df_m6_activ)
# rio::export(
#   corr_table_m6_activ,
#   file.path("output/tables/final", paste0(format(Sys.Date(), "%y%m%d"), "_kendall_aktivita_vs_clinical_M6.xlsx"))
# )
# 
# corr_table_delta_activ <- format_corr_table(cor_df_delta_activ)
# rio::export(
#   corr_table_delta_activ,
#   file.path("output/tables/final", paste0(format(Sys.Date(), "%y%m%d"), "_kendall_aktivita_vs_clinical_DELTA_M6_M0.xlsx"))
# )


# **********************************************************************
### ČAS ----
# **********************************************************************
#### Data preparation ----
# **********************************************************************
x_cols_corr_time <- d10_cas_b |>
  dplyr::select(where(is.numeric)) |>
  names()

y_cols_corr_time <- d11_clinics_imp |>
  dplyr::select(where(is.numeric)) |>
  names()

d13_join_all_time <- d13_join

d13_join_m0_time <- d13_join |>
  dplyr::filter(tolower(as.character(exam_order)) == "m0")

d13_join_m6_time <- d13_join |>
  dplyr::filter(tolower(as.character(exam_order)) == "m6")

delta_df_time <- d13_join |>
  dplyr::mutate(exam_order = tolower(as.character(exam_order))) |>
  dplyr::filter(exam_order %in% c("m0", "m6")) |>
  dplyr::select(immet_id, exam_order, all_of(x_cols_corr_time), all_of(y_cols_corr_time)) |>
  tidyr::pivot_wider(
    names_from  = exam_order,
    values_from = c(all_of(x_cols_corr_time), all_of(y_cols_corr_time))
  )

# keep only complete pairs (both m0 and m6 exist) - "soft" filter
delta_df_time <- delta_df_time |>
  dplyr::filter(
    if_all(dplyr::all_of(c(paste0(x_cols_corr_time, "_m0"), paste0(x_cols_corr_time, "_m6"))), ~ !is.na(.x)) |
      if_all(dplyr::all_of(c(paste0(y_cols_corr_time, "_m0"), paste0(y_cols_corr_time, "_m6"))), ~ !is.na(.x))
  )

# compute deltas into a new compact table with SAME column names as x_cols/y_cols
delta_X_time <- purrr::map_dfc(
  x_cols_corr_time,
  \(nm) delta_df_time[[paste0(nm, "_m6")]] - delta_df_time[[paste0(nm, "_m0")]]
) |>
  stats::setNames(x_cols_corr_time)

delta_Y_time <- purrr::map_dfc(
  y_cols_corr_time,
  \(nm) delta_df_time[[paste0(nm, "_m6")]] - delta_df_time[[paste0(nm, "_m0")]]
) |>
  stats::setNames(y_cols_corr_time)

d13_join_delta_time <- dplyr::bind_cols(
  tibble::tibble(immet_id = delta_df_time$immet_id),
  delta_X_time,
  delta_Y_time
) |> select(-age) |> 
  left_join(d13_join |>
              dplyr::filter(tolower(as.character(exam_order)) == "m0") |>
              dplyr::select(immet_id, age) |>
              dplyr::distinct())

# **********************************************************************
#### Analyses ----
# **********************************************************************
corr_res_all_time <- make_corr_heatmap_kendall(
  d_join = d13_join_all_time,
  x_cols = x_cols_corr_time,
  y_cols = y_cols_corr_time,
  title_suffix = "ČAS: ALL"
)
cor_df_all_time  <- corr_res_all_time$cor_df
(fig_cor_all_time <- corr_res_all_time$fig)

corr_res_m0_time <- make_corr_heatmap_kendall(
  d_join = d13_join_m0_time,
  x_cols = x_cols_corr_time,
  y_cols = y_cols_corr_time,
  title_suffix = "ČAS: M0"
)
cor_df_m0_time  <- corr_res_m0_time$cor_df
(fig_cor_m0_time <- corr_res_m0_time$fig)

corr_res_m6_time <- make_corr_heatmap_kendall(
  d_join = d13_join_m6_time,
  x_cols = x_cols_corr_time,
  y_cols = y_cols_corr_time,
  title_suffix = "ČAS: M6"
)
cor_df_m6_time  <- corr_res_m6_time$cor_df
(fig_cor_m6_time <- corr_res_m6_time$fig)

corr_res_delta_time <- make_corr_heatmap_kendall(
  d_join = d13_join_delta_time,
  x_cols = x_cols_corr_time,
  y_cols = y_cols_corr_time,
  title_suffix = "ČAS: DELTA (M6 - M0)"
)
cor_df_delta_time  <- corr_res_delta_time$cor_df
(fig_cor_delta_time <- corr_res_delta_time$fig)

# **********************************************************************
#### Export ----
# **********************************************************************
# figures
tiff(
  filename = file.path("output/figures/final", paste0(format(Sys.Date(), "%y%m%d"), "_corr_kendall_time_all.tiff")),
  compression = "lzw", res = "150", width = 1200, height = 1200
)
print(fig_cor_all_time)
dev.off()

tiff(
  filename = file.path("output/figures/final", paste0(format(Sys.Date(), "%y%m%d"), "_corr_kendall_time_m0.tiff")),
  compression = "lzw", res = "150", width = 1200, height = 1200
)
print(fig_cor_m0_time)
dev.off()

tiff(
  filename = file.path("output/figures/final", paste0(format(Sys.Date(), "%y%m%d"), "_corr_kendall_time_m6.tiff")),
  compression = "lzw", res = "150", width = 1200, height = 1200
)
print(fig_cor_m6_time)
dev.off()

tiff(
  filename = file.path("output/figures/final", paste0(format(Sys.Date(), "%y%m%d"), "_corr_kendall_time_delta_m6_m0.tiff")),
  compression = "lzw", res = "150", width = 1200, height = 1200
)
print(fig_cor_delta_time)
dev.off()

# # tables (xlsx)
# corr_table_all_time <- format_corr_table(cor_df_all_time)
# rio::export(
#   corr_table_all_time,
#   file.path("output/tables/final", paste0(format(Sys.Date(), "%y%m%d"), "_kendall_cas_vs_clinical_ALL.xlsx"))
# )
# 
# corr_table_m0_time <- format_corr_table(cor_df_m0_time)
# rio::export(
#   corr_table_m0_time,
#   file.path("output/tables/final", paste0(format(Sys.Date(), "%y%m%d"), "_kendall_cas_vs_clinical_M0.xlsx"))
# )
# 
# corr_table_m6_time <- format_corr_table(cor_df_m6_time)
# rio::export(
#   corr_table_m6_time,
#   file.path("output/tables/final", paste0(format(Sys.Date(), "%y%m%d"), "_kendall_cas_vs_clinical_M6.xlsx"))
# )
# 
# corr_table_delta_time <- format_corr_table(cor_df_delta_time)
# rio::export(
#   corr_table_delta_time,
#   file.path("output/tables/final", paste0(format(Sys.Date(), "%y%m%d"), "_kendall_cas_vs_clinical_DELTA_M6_M0.xlsx"))
# )


# **********************************************************************
## radar charts ----
# **********************************************************************
# list of covariates
fig_sel_01 <- c("enmo_full_recording_mean", "inactivity", "light",
                "moderate", "t400_time_min")

fig_sel_02 <- c("ig_gradient_enmo", "ig_intercept_enmo", "m5_enmo")

fig_sel_03 <- c("t100_b1m80_time_min", "t100_b5m80_time_min", "t100_b10m80_time_min")

fig_sel_04 <- c("p9931_enmo", "p9861_enmo", "p9792_enmo", "p9583_enmo", "p9167_enmo")

# models
res_adjust_01 <- d15_scores |> 
  select(immet_id, exam_order, age_m0, bmi, bcm, ck, physician_vas, b_m_rate,
         all_of(c(fig_sel_01, fig_sel_02, fig_sel_03))) |> 
  pivot_longer(cols = all_of(c(fig_sel_01, fig_sel_02, fig_sel_03)),
               names_to = "covar_name",
               values_to = "covar_value") |> 
  group_by(covar_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod_0 = map(data, ~possLMER(data = .x, covar_value ~ exam_order + age_m0 + (1|immet_id))),
         mod_1 = map(data, ~possLMER(data = .x, covar_value ~ exam_order + age_m0 + bmi + (1|immet_id))),
         mod_2 = map(data, ~possLMER(data = .x, covar_value ~ exam_order + age_m0 + ck + (1|immet_id))),
         mod_3 = map(data, ~possLMER(data = .x, covar_value ~ exam_order + age_m0 + bcm + (1|immet_id))),
         mod_4 = map(data, ~possLMER(data = .x, covar_value ~ exam_order + age_m0 + b_m_rate + (1|immet_id))),
         mod_5 = map(data, ~possLMER(data = .x, covar_value ~ exam_order + age_m0 + physician_vas + (1|immet_id))),
         fig = pmap(list(mod_0, mod_1, mod_2, mod_3, mod_4, mod_5, covar_name), 
                    ~plot(compare_performance(..1, ..2, ..3, ..4, ..5, ..6, rank = TRUE, verbose = FALSE))+
                      ggtitle(..7)),
         tab = pmap(list(mod_0, mod_1, mod_2, mod_3, mod_4, mod_5, covar_name), 
                    ~compare_performance(..1, ..2, ..3, ..4, ..5, ..6, rank = TRUE, verbose = FALSE))
  )

# export - figures
pdf(file.path("output/figures", 
              paste0(format(Sys.Date(), "%y%m%d"),
                     "_", 
                     "adjust_01.pdf")), 
    width = 10, height = 7)
walk(res_adjust_01$fig, print)
dev.off()

# export - table
res_adjust_01_tab <- res_adjust_01 |> 
  unnest(tab) |> 
  select(-(data:fig))

export(res_adjust_01_tab, file.path("output/tables", paste0(out_date, "_adjust_model_performance.xlsx")))

# **********************************************************************
## boat stacked chart ----
# **********************************************************************
# 2) long jen pro total/b1/b5/b10/rest u definovaných ranges
bouts_long <- d14_bouts |>
  mutate(
    exam_order      = factor(tolower(as.character(exam_order)), levels = c("m0", "m6")),
    disease_subtype = as.factor(disease_subtype)
  ) |>
  add_rest_cols(ranges = ranges, bouts = bouts) |>
  pivot_longer(
    cols = matches("^(t20_t40|t40_t100|t100_t400|t400_t8000)_(b1|b5|b10|rest)?_?time$"),
    names_to  = "var",
    values_to = "minutes"
  ) |>
  mutate(
    range = str_extract(var, "^(t20_t40|t40_t100|t100_t400|t400_t8000)"),
    component = case_when(
      str_detect(var, "_b1_time$")   ~ "b1",
      str_detect(var, "_b5_time$")   ~ "b5",
      str_detect(var, "_b10_time$")  ~ "b10",
      str_detect(var, "_rest_time$") ~ "rest",
      str_detect(var, "_time$")      ~ "total",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(component))

# 3) vyhodit "total" a (volitelně) agregovat minuty po skupinách
bouts_sum <- bouts_long |>
  filter(component != "total") |>
  mutate(component = factor(component, levels = c("b1", "b5", "b10", "rest")),
         range = factor(range, 
                        levels = c("t20_t40", "t40_t100", "t100_t400", "t400_t8000")
         )) |>
  group_by(exam_order, disease_subtype, range, component) |>
  summarise(minutes = sum(minutes, na.rm = TRUE), .groups = "drop")

# 4) plot: vždy 100 % v každém sloupci (position="fill")
fig_stack_subtype <- ggplot(bouts_sum, aes(x = exam_order, y = minutes, fill = component)) +
  geom_col(position = "fill", width = 0.8) +
  facet_grid(range ~ disease_subtype) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Exam",
    y = "Podíl z total času v pásmu (100 % = tXX_YY_time)",
    fill = "bout"
  ) +
  sjPlot::theme_sjplot2() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
    axis.title = ggplot2::element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 20, face = "bold"), 
    strip.background = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(size = 15),
    axis.text = ggplot2::element_text(size = 15, face = "bold")
  ) +
  scale_fill_manual(values = list(`b1` = "#c6dbef", `b5` = "#6baed6", 
                                  `b10` = "#2171b5",`rest` = "#bdbdbd"))


tiff(
  filename = file.path("output/figures/final",
                       paste0(format(Sys.Date(), "%y%m%d"),
                              "_", 
                              "stacked_subtype_01.tiff")),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_stack_subtype)
dev.off()

# 5) plot: vždy 100 % v každém sloupci (position="fill")
bouts_sum_2 <- bouts_long |>
  filter(component != "total") |>
  mutate(component = factor(component, levels = c("b1", "b5", "b10", "rest")),
         range = factor(range, 
                        levels = c("t20_t40", "t40_t100", "t100_t400", "t400_t8000")
         )) |>  group_by(exam_order, response_m0_m6, range, component) |>
  summarise(minutes = sum(minutes, na.rm = TRUE), .groups = "drop")

fig_stack_response <- ggplot(bouts_sum_2, aes(x = exam_order, y = minutes, fill = component)) +
  geom_col(position = "fill", width = 0.8) +
  facet_grid(range ~ response_m0_m6) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Exam",
    y = "Podíl z total času v pásmu (100 % = tXX_YY_time)",
    fill = "bout"
  ) +
  sjPlot::theme_sjplot2() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
    axis.title = ggplot2::element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 20, face = "bold"), 
    strip.background = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(size = 15),
    axis.text = ggplot2::element_text(size = 15, face = "bold")
  ) +
  scale_fill_manual(values = list(`b1` = "#c6dbef", `b5` = "#6baed6", 
                                  `b10` = "#2171b5",`rest` = "#bdbdbd"))

# export - figure
tiff(
  filename = file.path("output/figures/final",
                       paste0(format(Sys.Date(), "%y%m%d"),
                              "_", 
                              "stacked_response_01.tiff")),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_stack_response)
dev.off()

# export - table
export(list(subtype = bouts_sum, response = bouts_sum_2), 
       file.path("output/tables", paste0(out_date, "_stacked_bout_data_01.xlsx")))

# *******************************************************
# 4a. rCCA ----
# **********************************************************************
## Data preparation ----
# **********************************************************************
# ACTIVITA: X = numeric columns from d10_aktivita only
x_cols_cca_activ <- d10_aktivita |>
  dplyr::select(where(is.numeric)) |>
  names()

X_cca_activ <- d12_join |>
  dplyr::select(all_of(x_cols_cca_activ)) |>
  as.matrix()

Y_df_cca_activ <- d12_join |>
  dplyr::select(-c(immet_id, exam_order), -all_of(x_cols_cca_activ))

Y_cca_activ <- model.matrix(~ . - 1, data = Y_df_cca_activ)

stopifnot(nrow(X_cca_activ) == nrow(Y_cca_activ))

X_cca_activ <- drop_const_cols(X_cca_activ)
Y_cca_activ <- drop_const_cols(Y_cca_activ)

X_cca_activ <- scale(X_cca_activ)
Y_cca_activ <- scale(Y_cca_activ)

# ČAS: X = numeric columns from d10_cas_b only
x_cols_cca_time <- d10_cas_b |>
  dplyr::select(where(is.numeric)) |>
  names()

X_cca_time <- d13_join |>
  dplyr::select(all_of(x_cols_cca_time)) |>
  as.matrix()

Y_df_cca_time <- d13_join |>
  dplyr::select(-c(immet_id, exam_order), -all_of(x_cols_cca_time))

Y_cca_time <- model.matrix(~ . - 1, data = Y_df_cca_time)

stopifnot(nrow(X_cca_time) == nrow(Y_cca_time))

X_cca_time <- drop_const_cols(X_cca_time)
Y_cca_time <- drop_const_cols(Y_cca_time)

X_cca_time <- scale(X_cca_time)
Y_cca_time <- scale(Y_cca_time)

# MYOKINY: myokines (already delta-expression, measured only at M0)
# 1) X block: 
x_cols_cca_myok <- d17_myokines |>
  dplyr::select(-immet_id) |>
  names()

X_cca_myok <- d17_myokines |>
  dplyr::select(all_of(x_cols_cca_myok)) |>
  as.matrix()

# 2) Y block: clinical at M0 only
d12_m0 <- d12_join |>
  dplyr::filter(tolower(as.character(exam_order)) == "m0") |>
  dplyr::distinct(immet_id, .keep_all = TRUE)

# IMPORTANT: align rows between X and Y via immet_id
common_ids_myok <- intersect(d17_myokines$immet_id, d12_m0$immet_id)

X_cca_myok <- d17_myokines |>
  dplyr::filter(immet_id %in% common_ids_myok) |>
  dplyr::arrange(match(immet_id, common_ids_myok)) |>
  dplyr::select(all_of(x_cols_cca_myok)) |>
  as.matrix()

Y_df_cca_myok <- d12_m0 |>
  dplyr::filter(immet_id %in% common_ids_myok) |>
  dplyr::arrange(match(immet_id, common_ids_myok)) |>
  dplyr::select(-c(immet_id, exam_order))

Y_cca_myok <- model.matrix(~ . - 1, data = Y_df_cca_myok)

stopifnot(nrow(X_cca_myok) == nrow(Y_cca_myok))

X_cca_myok <- drop_const_cols(X_cca_myok)
Y_cca_myok <- drop_const_cols(Y_cca_myok)

X_cca_myok <- scale(X_cca_myok)
Y_cca_myok <- scale(Y_cca_myok)

# **********************************************************************
## Analyses (rCCA fits) ----
# **********************************************************************

# activity vs clinical (activity-specific)
rcc_fit_activ <- mixOmics::rcc(
  X_cca_activ, Y_cca_activ,
  ncomp  = 2,
  method = "shrinkage" # or "ridge"
)

# time vs clinical (time-specific)
rcc_fit_time <- mixOmics::rcc(
  X_cca_time, Y_cca_time,
  ncomp  = 2,
  method = "shrinkage" # or "ridge"
)

# Optional quick look at top loadings
sort(rcc_fit_activ$loadings$X[, 1], decreasing = TRUE)[1:10]
sort(rcc_fit_activ$loadings$Y[, 1], decreasing = TRUE)[1:10]

sort(rcc_fit_time$loadings$X[, 1], decreasing = TRUE)[1:10]
sort(rcc_fit_time$loadings$Y[, 1], decreasing = TRUE)[1:10]

# MYOKINES
rcc_fit_myok <- mixOmics::rcc(
  X_cca_myok, Y_cca_myok,
  ncomp  = 2,
  method = "shrinkage"
)

# Optional quick look at top loadings
sort(rcc_fit_myok$loadings$X[, 1], decreasing = TRUE)[1:10]
sort(rcc_fit_myok$loadings$Y[, 1], decreasing = TRUE)[1:10]

# **********************************************************************
## AKTIVITA outputs ----
# **********************************************************************
### Canonical scores scatter (comp1) ----
# **********************************************************************
scores_df_activ <- tibble(
  immet_id   = d12_join$immet_id,
  exam_order = d12_join$exam_order,
  X1 = rcc_fit_activ$variates$X[, 1],
  Y1 = rcc_fit_activ$variates$Y[, 1],
  X2 = rcc_fit_activ$variates$X[, 2],
  Y2 = rcc_fit_activ$variates$Y[, 2]
)

(fig_can_comp_activ_01 <- ggplot(scores_df_activ, aes(X1, Y1, color = exam_order)) +
    geom_point(alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linetype = 2) +
    labs(
      title = "rCCA: Canonical variates (Component 1)",
      x = "X variate 1 (activity block)",
      y = "Y variate 1 (clinical block)",
      color = "Exam"
    ) +
    sjPlot::theme_sjplot2())

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_canon_comp_activ_01.tiff")
  ),
  compression = "lzw",
  res = "300",
  width = 1600,
  height = 1200
)
print(fig_can_comp_activ_01)
dev.off()

# **********************************************************************
### Heatmap: cross-correlation top loadings ----
# **********************************************************************
topX_activ <- tibble(
  var = colnames(X_cca_activ),
  loading = rcc_fit_activ$loadings$X[, 1]
) |>
  arrange(desc(abs(loading))) |>
  slice(1:15)

topX_vars_activ <- topX_activ |> pull(var) |> as.character()

topY_activ <- tibble(
  var = colnames(Y_cca_activ),
  loading = rcc_fit_activ$loadings$Y[, 1]
) |>
  arrange(desc(abs(loading))) |>
  slice(1:20)

topY_vars_activ <- topY_activ |> pull(var) |> as.character()

cor_mat_activ <- cor(
  X_cca_activ[, topX_vars_activ, drop = FALSE],
  Y_cca_activ[, topY_vars_activ, drop = FALSE]
)

cor_df_activ <- as.data.frame(as.table(cor_mat_activ)) |>
  setNames(c("x_var", "y_var", "cor"))

(fig_heatmap_loading_activ_01 <- ggplot(cor_df_activ, aes(x_var, y_var, fill = cor)) +
    geom_tile() +
    coord_equal() +
    labs(
      title = "Cross-correlation heatmap (top X vs top Y)",
      x = "Activity (top loadings)",
      y = "Clinical (top loadings)",
      fill = "cor"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradient2(low="blue", mid="white", high="red", space ="Lab"))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_heatmap_loading_activ_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_heatmap_loading_activ_01)
dev.off()

# **********************************************************************
### Latent axis loadings (barplots) ----
# **********************************************************************
topX_tab_activ <- tibble(
  var = colnames(X_cca_activ),
  loading = rcc_fit_activ$loadings$X[, 1]
) |>
  arrange(desc(abs(loading))) |>
  slice(1:20)

topY_tab_activ <- tibble(
  var = colnames(Y_cca_activ),
  loading = rcc_fit_activ$loadings$Y[, 1]
) |>
  arrange(desc(abs(loading))) |>
  slice(1:30)

load_df_activ <- bind_rows(
  topX_tab_activ |> mutate(type = "Activity"),
  topY_tab_activ |> mutate(type = "Clinical")
)

(fig_act_clin_axis_01 <- ggplot(load_df_activ, aes(x = loading, y = reorder(var, loading), fill = type)) +
    geom_col(alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~type, scales = "free_y") +
    labs(
      x = "Loading (Canonical component 1)",
      y = NULL,
      title = "Opposing poles of the latent activity–clinical axis",
      subtitle = "Variables on opposite sides represent contrasting disease features"
    ) +
    theme_minimal(base_size = 13))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_activity_clinic_axis_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_act_clin_axis_01)
dev.off()

# **********************************************************************
### Correlation with canonical axis (tiles) ----
# **********************************************************************
cor_df_axis_activ <- bind_rows(
  tibble(
    var = colnames(X_cca_activ),
    cor = cor(X_cca_activ, rcc_fit_activ$variates$X[, 1]),
    type = "Activity"
  ),
  tibble(
    var = colnames(Y_cca_activ),
    cor = cor(Y_cca_activ, rcc_fit_activ$variates$Y[, 1]),
    type = "Clinical"
  )
) |>
  filter(abs(cor) > 0.15)

(fig_canon_cor_activ_01 <- ggplot(cor_df_axis_activ, aes(x = type, y = reorder(var, cor), fill = cor)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    labs(
      fill = "Correlation",
      title = "Correlation of individual variables with the latent CCA axis",
      subtitle = "Red vs blue indicate opposing disease patterns"
    ) +
    theme_minimal(base_size = 13))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_canon_cor_activ_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_canon_cor_activ_01)
dev.off()

# **********************************************************************
### Subtype and response visualization ----
# **********************************************************************
plot_df_activ <- d12_join |>
  transmute(
    immet_id, exam_order,
    subtype = as.factor(disease_subtype),
    response = as.factor(response_m0_m6),
    compX1 = rcc_fit_activ$variates$X[, 1],
    compY1 = rcc_fit_activ$variates$Y[, 1]
  )

(fig_rcca_subtype_activ_01 <- ggplot(plot_df_activ, aes(compX1, compY1, color = subtype)) +
    geom_point(size = 2.8, alpha = 0.85) +
    stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
    labs(
      x = "rCCA component 1 (Activity)",
      y = "rCCA component 1 (Clinical)",
      title = "Patients projected onto the shared activity–clinical space",
      subtitle = "Ellipses highlight subtype clustering (80% t-ellipse)"
    ) +
    theme_minimal(base_size = 13))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_rcca_subtype_activ_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_rcca_subtype_activ_01)
dev.off()

(fig_rcca_response_activ_01 <- ggplot(plot_df_activ, aes(compX1, compY1, color = response)) +
    geom_point(size = 2.8, alpha = 0.85) +
    stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
    labs(
      x = "rCCA component 1 (Activity)",
      y = "rCCA component 1 (Clinical)",
      title = "Patients projected onto the shared activity–clinical space",
      subtitle = "Ellipses highlight subtype clustering (80% t-ellipse)"
    ) +
    theme_minimal(base_size = 13))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_rcca_subtype_response_activ_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_rcca_response_activ_01)
dev.off()

# **********************************************************************
## ČAS outputs  ----
# **********************************************************************
### Canonical scores scatter (comp1) ----
# **********************************************************************
scores_df_time <- tibble(
  immet_id   = d13_join$immet_id,
  exam_order = d13_join$exam_order,
  X1 = rcc_fit_time$variates$X[, 1],
  Y1 = rcc_fit_time$variates$Y[, 1],
  X2 = rcc_fit_time$variates$X[, 2],
  Y2 = rcc_fit_time$variates$Y[, 2]
)

(fig_can_comp_time_01 <- ggplot(scores_df_time, aes(X1, Y1, color = exam_order)) +
    geom_point(alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linetype = 2) +
    labs(
      title = "rCCA: Canonical variates (Component 1)",
      x = "X variate 1 (time block)",
      y = "Y variate 1 (clinical block)",
      color = "Exam"
    ) +
    sjPlot::theme_sjplot2())

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_canon_comp_time_01.tiff")
  ),
  compression = "lzw",
  res = "300",
  width = 1600,
  height = 1200
)
print(fig_can_comp_time_01)
dev.off()

# **********************************************************************
### Heatmap: cross-correlation top loadings ----
# **********************************************************************
topX_time <- tibble(
  var = colnames(X_cca_time),
  loading = rcc_fit_time$loadings$X[, 1]
) |>
  arrange(desc(abs(loading))) |>
  slice(1:15)

topX_vars_time <- topX_time |> pull(var) |> as.character()

topY_time <- tibble(
  var = colnames(Y_cca_time),
  loading = rcc_fit_time$loadings$Y[, 1]
) |>
  arrange(desc(abs(loading))) |>
  slice(1:20)

topY_vars_time <- topY_time |> pull(var) |> as.character()

cor_mat_time <- cor(
  X_cca_time[, topX_vars_time, drop = FALSE],
  Y_cca_time[, topY_vars_time, drop = FALSE]
)

cor_df_time <- as.data.frame(as.table(cor_mat_time)) |>
  setNames(c("x_var", "y_var", "cor"))

(fig_heatmap_loading_time_01 <- ggplot(cor_df_time, aes(x_var, y_var, fill = cor)) +
    geom_tile() +
    coord_equal() +
    labs(
      title = "Cross-correlation heatmap (top X vs top Y)",
      x = "Time (top loadings)",
      y = "Clinical (top loadings)",
      fill = "cor"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradient2(low="blue", mid="white", high="red", space ="Lab"))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_heatmap_loading_time_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_heatmap_loading_time_01)
dev.off()

# **********************************************************************
### Latent axis loadings (barplots) ----
# **********************************************************************
topX_tab_time <- tibble(
  var = colnames(X_cca_time),
  loading = rcc_fit_time$loadings$X[, 1]
) |>
  arrange(desc(abs(loading))) |>
  slice(1:20)

topY_tab_time <- tibble(
  var = colnames(Y_cca_time),
  loading = rcc_fit_time$loadings$Y[, 1]
) |>
  arrange(desc(abs(loading))) |>
  slice(1:30)

load_df_time <- bind_rows(
  topX_tab_time |> mutate(type = "Time"),
  topY_tab_time |> mutate(type = "Clinical")
)

(fig_time_clin_axis_01 <- ggplot(load_df_time, aes(x = loading, y = reorder(var, loading), fill = type)) +
    geom_col(alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~type, scales = "free_y") +
    labs(
      x = "Loading (Canonical component 1)",
      y = NULL,
      title = "Opposing poles of the latent time–clinical axis",
      subtitle = "Variables on opposite sides represent contrasting disease features"
    ) +
    theme_minimal(base_size = 13))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_time_clinic_axis_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_time_clin_axis_01)
dev.off()

# **********************************************************************
### Correlation with canonical axis (tiles) ----
# **********************************************************************
cor_df_axis_time <- bind_rows(
  tibble(
    var = colnames(X_cca_time),
    cor = cor(X_cca_time, rcc_fit_time$variates$X[, 1]),
    type = "Time"
  ),
  tibble(
    var = colnames(Y_cca_time),
    cor = cor(Y_cca_time, rcc_fit_time$variates$Y[, 1]),
    type = "Clinical"
  )
) |>
  filter(abs(cor) > 0.15)

(fig_canon_cor_time_01 <- ggplot(cor_df_axis_time, aes(x = type, y = reorder(var, cor), fill = cor)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    labs(
      fill = "Correlation",
      title = "Correlation of individual variables with the latent CCA axis",
      subtitle = "Red vs blue indicate opposing disease patterns"
    ) +
    theme_minimal(base_size = 13))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_canon_cor_time_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_canon_cor_time_01)
dev.off()

# **********************************************************************
### Subtype and response visualization ----
# **********************************************************************
plot_df_time <- d13_join |>
  transmute(
    immet_id, exam_order,
    subtype = as.factor(disease_subtype),
    response = as.factor(response_m0_m6),
    compX1 = rcc_fit_time$variates$X[, 1],
    compY1 = rcc_fit_time$variates$Y[, 1]
  )

(fig_rcca_response_time_01 <- ggplot(plot_df_time, aes(compX1, compY1, color = response)) +
    geom_point(size = 2.8, alpha = 0.85) +
    stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
    labs(
      x = "rCCA component 1 (Time)",
      y = "rCCA component 1 (Clinical)",
      title = "Patients projected onto the shared time–clinical space",
      subtitle = "Ellipses highlight response clustering (80% t-ellipse)"
    ) +
    theme_minimal(base_size = 13))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_rcca_response_time_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_rcca_response_time_01)
dev.off()

(fig_rcca_subtype_time_01 <- ggplot(plot_df_time, aes(compX1, compY1, color = subtype)) +
    geom_point(size = 2.8, alpha = 0.85) +
    stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
    labs(
      x = "rCCA component 1 (Time)",
      y = "rCCA component 1 (Clinical)",
      title = "Patients projected onto the shared time–clinical space",
      subtitle = "Ellipses highlight response clustering (80% t-ellipse)"
    ) +
    theme_minimal(base_size = 13))


# **********************************************************************
## MYOKINY outputs ----
# **********************************************************************
### Canonical scores scatter (comp1) ----
# **********************************************************************

scores_df_myok <- tibble(
  immet_id = common_ids_myok,
  X1 = rcc_fit_myok$variates$X[, 1],
  Y1 = rcc_fit_myok$variates$Y[, 1],
  X2 = rcc_fit_myok$variates$X[, 2],
  Y2 = rcc_fit_myok$variates$Y[, 2]
) |>
  left_join(
    d12_m0 |>
      dplyr::select(immet_id, disease_subtype, response_m0_m6),
    by = "immet_id"
  )

(fig_can_comp_myok_01 <- ggplot(scores_df_myok, aes(X1, Y1)) +
    geom_point(alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, linetype = 2) +
    labs(
      title = "rCCA: Canonical variates (Component 1)",
      x = "X variate 1 (myokines block, M0)",
      y = "Y variate 1 (clinical block, M0)"
    ) +
    sjPlot::theme_sjplot2())

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_canon_comp_myok_01.tiff")
  ),
  compression = "lzw",
  res = "300",
  width = 1600,
  height = 1200
)
print(fig_can_comp_myok_01)
dev.off()

# **********************************************************************
### Heatmap: cross-correlation top loadings ----
# **********************************************************************

topX_myok <- tibble(
  var = colnames(X_cca_myok),
  loading = rcc_fit_myok$loadings$X[, 1]
) |>
  arrange(desc(abs(loading))) |>
  slice(1:15)

topX_vars_myok <- topX_myok |> pull(var) |> as.character()

topY_myok <- tibble(
  var = colnames(Y_cca_myok),
  loading = rcc_fit_myok$loadings$Y[, 1]
) |>
  arrange(desc(abs(loading))) |>
  slice(1:20)

topY_vars_myok <- topY_myok |> pull(var) |> as.character()

cor_mat_myok <- cor(
  X_cca_myok[, topX_vars_myok, drop = FALSE],
  Y_cca_myok[, topY_vars_myok, drop = FALSE]
)

cor_df_myok <- as.data.frame(as.table(cor_mat_myok)) |>
  setNames(c("x_var", "y_var", "cor"))

(fig_heatmap_loading_myok_01 <- ggplot(cor_df_myok, aes(x_var, y_var, fill = cor)) +
    geom_tile() +
    coord_equal() +
    labs(
      title = "Cross-correlation heatmap (top X vs top Y)",
      x = "Myokines (top loadings)",
      y = "Clinical (top loadings)",
      fill = "cor"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradient2(low="blue", mid="white", high="red", space ="Lab"))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_heatmap_loading_myok_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_heatmap_loading_myok_01)
dev.off()

# **********************************************************************
### Latent axis loadings (barplots) ----
# **********************************************************************

topX_tab_myok <- tibble(
  var = colnames(X_cca_myok),
  loading = rcc_fit_myok$loadings$X[, 1]
) |>
  arrange(desc(abs(loading))) |>
  slice(1:20)

topY_tab_myok <- tibble(
  var = colnames(Y_cca_myok),
  loading = rcc_fit_myok$loadings$Y[, 1]
) |>
  arrange(desc(abs(loading))) |>
  slice(1:30)

load_df_myok <- bind_rows(
  topX_tab_myok |> mutate(type = "Myokines"),
  topY_tab_myok |> mutate(type = "Clinical")
)

(fig_myok_clin_axis_01 <- ggplot(load_df_myok, aes(x = loading, y = reorder(var, loading), fill = type)) +
    geom_col(alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~type, scales = "free_y") +
    labs(
      x = "Loading (Canonical component 1)",
      y = NULL,
      title = "Opposing poles of the latent myokines–clinical axis (M0)",
      subtitle = "Variables on opposite sides represent contrasting disease features"
    ) +
    theme_minimal(base_size = 13))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_myokines_clinic_axis_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_myok_clin_axis_01)
dev.off()

# **********************************************************************
### Correlation with canonical axis (tiles) ----
# **********************************************************************

cor_df_axis_myok <- bind_rows(
  tibble(
    var = colnames(X_cca_myok),
    cor = cor(X_cca_myok, rcc_fit_myok$variates$X[, 1]),
    type = "Myokines"
  ),
  tibble(
    var = colnames(Y_cca_myok),
    cor = cor(Y_cca_myok, rcc_fit_myok$variates$Y[, 1]),
    type = "Clinical"
  )
) |>
  filter(abs(cor) > 0.15)

(fig_canon_cor_myok_01 <- ggplot(cor_df_axis_myok, aes(x = type, y = reorder(var, cor), fill = cor)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    labs(
      fill = "Correlation",
      title = "Correlation of individual variables with the latent CCA axis (M0)",
      subtitle = "Red vs blue indicate opposing disease patterns"
    ) +
    theme_minimal(base_size = 13))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_canon_cor_myok_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_canon_cor_myok_01)
dev.off()

# **********************************************************************
### Subtype and response visualization ----
# **********************************************************************

(fig_rcca_subtype_myok_01 <- ggplot(scores_df_myok, aes(X1, Y1, color = as.factor(disease_subtype))) +
   geom_point(size = 2.8, alpha = 0.85) +
   stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
   labs(
     x = "rCCA component 1 (Myokines, M0)",
     y = "rCCA component 1 (Clinical, M0)",
     title = "Patients projected onto the shared myokines–clinical space (M0)",
     subtitle = "Ellipses highlight subtype clustering (80% t-ellipse)",
     color = "Subtype"
   ) +
   theme_minimal(base_size = 13))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_rcca_subtype_myok_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_rcca_subtype_myok_01)
dev.off()

(fig_rcca_response_myok_01 <- ggplot(scores_df_myok, aes(X1, Y1, color = as.factor(response_m0_m6))) +
    geom_point(size = 2.8, alpha = 0.85) +
    stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
    labs(
      x = "rCCA component 1 (Myokines, M0)",
      y = "rCCA component 1 (Clinical, M0)",
      title = "Patients projected onto the shared myokines–clinical space (M0)",
      subtitle = "Ellipses highlight response clustering (80% t-ellipse)",
      color = "Response"
    ) +
    theme_minimal(base_size = 13))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_rcca_response_myok_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_rcca_response_myok_01)
dev.off()

# **********************************************************************
# 4b. rCCA ----
# **********************************************************************
## Data preparation ----
# **********************************************************************

# 0) MYOKINY (X block; měřeno jen M0)

d17_myok_u <- d17_myokines |>
  dplyr::distinct(immet_id, .keep_all = TRUE)

x_cols_myok <- d17_myok_u |>
  dplyr::select(-immet_id) |>
  dplyr::select(where(is.numeric)) |>
  names()


# 1) AKTIVITA v M0 (Y block pro myokiny × aktivita)
exclude_activ_vars <- c("m1_enmo", "enmo_full_recording_mean")

x_cols_activ <- d10_aktivita |>
  dplyr::select(where(is.numeric)) |>
  dplyr::select(-dplyr::any_of(exclude_activ_vars)) |>
  names()

d12_m0 <- d12_join |>
  dplyr::mutate(exam_order = tolower(as.character(exam_order))) |>
  dplyr::filter(exam_order == "m0") |>
  dplyr::distinct(immet_id, .keep_all = TRUE)

common_ids_myok_activ <- intersect(d17_myok_u$immet_id, d12_m0$immet_id)

# X = myokiny (M0)
X_myok_activ <- d17_myok_u |>
  dplyr::filter(immet_id %in% common_ids_myok_activ) |>
  dplyr::arrange(match(immet_id, common_ids_myok_activ)) |>
  dplyr::select(all_of(x_cols_myok)) |>
  as.matrix()

# Y = aktivita (M0)
Y_activ_myok <- d12_m0 |>
  dplyr::filter(immet_id %in% common_ids_myok_activ) |>
  dplyr::arrange(match(immet_id, common_ids_myok_activ)) |>
  dplyr::select(all_of(x_cols_activ)) |>
  as.matrix()

meta_myok_activ <- d12_m0 |>
  dplyr::filter(immet_id %in% common_ids_myok_activ) |>
  dplyr::arrange(match(immet_id, common_ids_myok_activ)) |>
  dplyr::select(immet_id, disease_subtype, response_m0_m6, sex)

stopifnot(nrow(X_myok_activ) == nrow(Y_activ_myok))

X_myok_activ <- drop_const_cols(X_myok_activ)
Y_activ_myok <- drop_const_cols(Y_activ_myok)

# vyhoď řádky s NA (mixOmics::rcc NA neumí)
keep_activ <- complete.cases(X_myok_activ, Y_activ_myok)
X_myok_activ <- X_myok_activ[keep_activ, , drop = FALSE]
Y_activ_myok <- Y_activ_myok[keep_activ, , drop = FALSE]
meta_myok_activ <- meta_myok_activ[keep_activ, , drop = FALSE]

X_myok_activ <- scale(X_myok_activ)
Y_activ_myok <- scale(Y_activ_myok)

# 2) ČAS v M0 (Y block pro myokiny × čas)
x_cols_time <- d10_cas_b |>
  dplyr::select(where(is.numeric)) |>
  names()

d13_m0 <- d13_join |>
  dplyr::mutate(exam_order = tolower(as.character(exam_order))) |>
  dplyr::filter(exam_order == "m0") |>
  dplyr::distinct(immet_id, .keep_all = TRUE)

common_ids_myok_time <- intersect(d17_myok_u$immet_id, d13_m0$immet_id)

# X = myokiny (M0)
X_myok_time <- d17_myok_u |>
  dplyr::filter(immet_id %in% common_ids_myok_time) |>
  dplyr::arrange(match(immet_id, common_ids_myok_time)) |>
  dplyr::select(all_of(x_cols_myok)) |>
  as.matrix()

# Y = čas (M0)
Y_time_myok <- d13_m0 |>
  dplyr::filter(immet_id %in% common_ids_myok_time) |>
  dplyr::arrange(match(immet_id, common_ids_myok_time)) |>
  dplyr::select(all_of(x_cols_time)) |>
  as.matrix()

meta_myok_time <- d13_m0 |>
  dplyr::filter(immet_id %in% common_ids_myok_time) |>
  dplyr::arrange(match(immet_id, common_ids_myok_time)) |>
  dplyr::select(immet_id, disease_subtype, response_m0_m6, sex)

stopifnot(nrow(X_myok_time) == nrow(Y_time_myok))

X_myok_time <- drop_const_cols(X_myok_time)
Y_time_myok <- drop_const_cols(Y_time_myok)

keep_time <- complete.cases(X_myok_time, Y_time_myok)
X_myok_time <- X_myok_time[keep_time, , drop = FALSE]
Y_time_myok <- Y_time_myok[keep_time, , drop = FALSE]
meta_myok_time <- meta_myok_time[keep_time, , drop = FALSE]

X_myok_time <- scale(X_myok_time)
Y_time_myok <- scale(Y_time_myok)

# **********************************************************************
## Analyses (rCCA fits) ----
# **********************************************************************

# 1) MYOKINY × AKTIVITA
rcc_fit_myok_activ <- mixOmics::rcc(
  X_myok_activ, Y_activ_myok,
  ncomp  = 2,
  method = "shrinkage"
)

# 2) MYOKINY × ČAS
rcc_fit_myok_time <- mixOmics::rcc(
  X_myok_time, Y_time_myok,
  ncomp  = 2,
  method = "shrinkage"
)

# Quick top loadings (comp1)
sort(rcc_fit_myok_activ$loadings$X[, 1], decreasing = TRUE)[1:10]
sort(rcc_fit_myok_activ$loadings$Y[, 1], decreasing = TRUE)[1:10]

sort(rcc_fit_myok_time$loadings$X[, 1], decreasing = TRUE)[1:10]
sort(rcc_fit_myok_time$loadings$Y[, 1], decreasing = TRUE)[1:10]

# **********************************************************************
## MYOKINY × AKTIVITA outputs ----
# **********************************************************************

scores_df_myok_activ <- tibble::tibble(
  immet_id = meta_myok_activ$immet_id,
  subtype  = as.factor(meta_myok_activ$disease_subtype),
  response = as.factor(meta_myok_activ$response_m0_m6),
  sex = as.factor(meta_myok_activ$sex),
  X1 = rcc_fit_myok_activ$variates$X[, 1],
  Y1 = rcc_fit_myok_activ$variates$Y[, 1],
  X2 = rcc_fit_myok_activ$variates$X[, 2],
  Y2 = rcc_fit_myok_activ$variates$Y[, 2]
)

(fig_myok_activ_comp1_subtype <- ggplot(scores_df_myok_activ, aes(X1, Y1, color = subtype)) +
    geom_point(size = 2.8, alpha = 0.85) +
    stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
    labs(
      x = "rCCA component 1 (Myokines, M0)",
      y = "rCCA component 1 (Activity, M0)",
      title = "Myokines × Activity (M0): projekce do sdíleného prostoru",
      subtitle = "Barvy = podtyp; elipsy = 80% t-ellipse"
    ) +
    theme_minimal(base_size = 13))

save_tiff(
  fig_myok_activ_comp1_subtype,
  file.path("output/figures", paste0(format(Sys.Date(), "%y%m%d"), "_rcca_myok_activ_subtype_01.tiff")),
  res = 150
)

(fig_myok_activ_comp1_response <- ggplot(scores_df_myok_activ, aes(X1, Y1, color = response)) +
    geom_point(size = 2.8, alpha = 0.85) +
    stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
    labs(
      x = "rCCA component 1 (Myokines, M0)",
      y = "rCCA component 1 (Activity, M0)",
      title = "Myokines × Activity (M0): vztah k léčebné odpovědi",
      subtitle = "Barvy = odpověď; elipsy = 80% t-ellipse"
    ) +
    theme_minimal(base_size = 13))

save_tiff(
  fig_myok_activ_comp1_response,
  file.path("output/figures", paste0(format(Sys.Date(), "%y%m%d"), "_rcca_myok_activ_response_01.tiff")),
  res = 150
)

(fig_myok_activ_comp1_sex <- ggplot(scores_df_myok_activ, aes(X1, Y1, color = sex)) +
    geom_point(size = 2.8, alpha = 0.85) +
    stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
    labs(
      x = "rCCA component 1 (Myokines, M0)",
      y = "rCCA component 1 (Activity, M0)",
      title = "Myokines × Activity (M0): projekce do sdíleného prostoru",
      subtitle = "Barvy = pohlaví; elipsy = 80% t-ellipse"
    ) +
    theme_minimal(base_size = 13))

save_tiff(
  fig_myok_activ_comp1_sex,
  file.path("output/figures", paste0(format(Sys.Date(), "%y%m%d"), "_rcca_myok_activ_sex_01.tiff")),
  res = 150
)

# Heatmap top loadings correlations
topX_myok_activ <- tibble::tibble(var = colnames(X_myok_activ), loading = rcc_fit_myok_activ$loadings$X[, 1]) |>
  dplyr::arrange(desc(abs(loading))) |>
  dplyr::slice(1:15) |>
  dplyr::pull(var) |>
  as.character()

topY_activ_myok <- tibble::tibble(var = colnames(Y_activ_myok), loading = rcc_fit_myok_activ$loadings$Y[, 1]) |>
  dplyr::arrange(desc(abs(loading))) |>
  dplyr::slice(1:20) |>
  dplyr::pull(var) |>
  as.character()

cor_mat_myok_activ <- cor(
  X_myok_activ[, topX_myok_activ, drop = FALSE],
  Y_activ_myok[, topY_activ_myok, drop = FALSE]
)

cor_df_myok_activ <- as.data.frame(as.table(cor_mat_myok_activ)) |>
  setNames(c("x_var", "y_var", "cor"))

(fig_heatmap_myok_activ <- ggplot(cor_df_myok_activ, aes(x_var, y_var, fill = cor)) +
    geom_tile() +
    coord_equal() +
    labs(
      title = "Myokines × Activity: cross-correlation (top loadings)",
      x = "Myokines (top loadings)",
      y = "Activity (top loadings)",
      fill = "cor"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", space = "Lab"))

save_tiff(
  fig_heatmap_myok_activ,
  file.path("output/figures", paste0(format(Sys.Date(), "%y%m%d"), "_heatmap_loading_myok_activ_01.tiff")),
  res = 150
)

# **********************************************************************
### MYOKINY × AKTIVITA: Latent axis loadings (barplots) ----
# **********************************************************************

loadX_myok_activ <- as.matrix(rcc_fit_myok_activ$loadings$X)[, 1]
loadY_activ_myok <- as.matrix(rcc_fit_myok_activ$loadings$Y)[, 1]

topX_tab_myok_activ <- tibble::tibble(
  var = colnames(X_myok_activ),
  loading = loadX_myok_activ
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice_head(n = 20)

topY_tab_myok_activ <- tibble::tibble(
  var = colnames(Y_activ_myok),
  loading = loadY_activ_myok
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice_head(n = 30)

load_df_myok_activ <- dplyr::bind_rows(
  topX_tab_myok_activ |> dplyr::mutate(type = "Myokines"),
  topY_tab_myok_activ |> dplyr::mutate(type = "Activity")
)

(fig_myok_activ_axis_01 <- ggplot2::ggplot(
  load_df_myok_activ,
  ggplot2::aes(x = loading, y = reorder(var, loading), fill = type)
) +
    ggplot2::geom_col(alpha = 0.85) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~type, scales = "free_y") +
    ggplot2::labs(
      x = "Loading (Canonical component 1)",
      y = NULL,
      title = "Opposing poles of the latent myokines–activity axis (M0)",
      subtitle = "Variables on opposite sides represent contrasting biological–functional patterns"
    ) +
    ggplot2::theme_minimal(base_size = 13))

save_tiff(
  fig_myok_activ_axis_01,
  file.path("output/figures", paste0(format(Sys.Date(), "%y%m%d"), "_myokines_activity_axis_01.tiff")),
  res = 150
)


# **********************************************************************
## MYOKINY × ČAS outputs ----
# **********************************************************************

scores_df_myok_time <- tibble::tibble(
  immet_id = meta_myok_time$immet_id,
  subtype  = as.factor(meta_myok_time$disease_subtype),
  response = as.factor(meta_myok_time$response_m0_m6),
  sex = as.factor(meta_myok_time$sex),
  X1 = rcc_fit_myok_time$variates$X[, 1],
  Y1 = rcc_fit_myok_time$variates$Y[, 1],
  X2 = rcc_fit_myok_time$variates$X[, 2],
  Y2 = rcc_fit_myok_time$variates$Y[, 2]
)

(fig_myok_time_comp1_subtype <- ggplot(scores_df_myok_time, aes(X1, Y1, color = subtype)) +
    geom_point(size = 2.8, alpha = 0.85) +
    stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
    labs(
      x = "rCCA component 1 (Myokines, M0)",
      y = "rCCA component 1 (Time, M0)",
      title = "Myokines × Time (M0): projekce do sdíleného prostoru",
      subtitle = "Barvy = podtyp; elipsy = 80% t-ellipse"
    ) +
    theme_minimal(base_size = 13))

save_tiff(
  fig_myok_time_comp1_subtype,
  file.path("output/figures/final", paste0(format(Sys.Date(), "%y%m%d"), "_rcca_myok_time_subtype_01.tiff")),
  res = 150
)

(fig_myok_time_comp1_response <- ggplot(scores_df_myok_time, aes(X1, Y1, color = response)) +
    geom_point(size = 2.8, alpha = 0.85) +
    stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
    labs(
      x = "rCCA component 1 (Myokines, M0)",
      y = "rCCA component 1 (Time, M0)",
      title = "Myokines × Time (M0): vztah k léčebné odpovědi",
      subtitle = "Barvy = odpověď; elipsy = 80% t-ellipse"
    ) +
    theme_minimal(base_size = 13))

save_tiff(
  fig_myok_time_comp1_response,
  file.path("output/figures", paste0(format(Sys.Date(), "%y%m%d"), "_rcca_myok_time_response_01.tiff")),
  res = 150
)

(fig_myok_time_comp1_sex <- ggplot(scores_df_myok_time, aes(X1, Y1, color = sex)) +
    geom_point(size = 2.8, alpha = 0.85) +
    stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
    labs(
      x = "rCCA component 1 (Myokines, M0)",
      y = "rCCA component 1 (Time, M0)",
      title = "Myokines × Time (M0): projekce do sdíleného prostoru",
      subtitle = "Barvy = pohlaví; elipsy = 80% t-ellipse"
    ) +
    theme_minimal(base_size = 13))

save_tiff(
  fig_myok_time_comp1_sex,
  file.path("output/figures/final", paste0(format(Sys.Date(), "%y%m%d"), "_rcca_myok_time_sex_01.tiff")),
  res = 150
)


# Heatmap top loadings correlations
topX_myok_time <- tibble::tibble(var = colnames(X_myok_time), loading = rcc_fit_myok_time$loadings$X[, 1]) |>
  dplyr::arrange(desc(abs(loading))) |>
  dplyr::slice(1:15) |>
  dplyr::pull(var) |>
  as.character()

topY_time_myok <- tibble::tibble(var = colnames(Y_time_myok), loading = rcc_fit_myok_time$loadings$Y[, 1]) |>
  dplyr::arrange(desc(abs(loading))) |>
  dplyr::slice(1:20) |>
  dplyr::pull(var) |>
  as.character()

cor_mat_myok_time <- cor(
  X_myok_time[, topX_myok_time, drop = FALSE],
  Y_time_myok[, topY_time_myok, drop = FALSE]
)

cor_df_myok_time <- as.data.frame(as.table(cor_mat_myok_time)) |>
  setNames(c("x_var", "y_var", "cor"))

(fig_heatmap_myok_time <- ggplot(cor_df_myok_time, aes(x_var, y_var, fill = cor)) +
    geom_tile() +
    coord_equal() +
    labs(
      title = "Myokines × Time: cross-correlation (top loadings)",
      x = "Myokines (top loadings)",
      y = "Time (top loadings)",
      fill = "cor"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", space = "Lab"))

save_tiff(
  fig_heatmap_myok_time,
  file.path("output/figures", paste0(format(Sys.Date(), "%y%m%d"), "_heatmap_loading_myok_time_01.tiff")),
  res = 150
)

# **********************************************************************
### MYOKINY × ČAS: Latent axis loadings (barplots) ----
# **********************************************************************

loadX_myok_time <- as.matrix(rcc_fit_myok_time$loadings$X)[, 1]
loadY_time_myok <- as.matrix(rcc_fit_myok_time$loadings$Y)[, 1]

topX_tab_myok_time <- tibble::tibble(
  var = colnames(X_myok_time),
  loading = loadX_myok_time
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice_head(n = 20)

topY_tab_myok_time <- tibble::tibble(
  var = colnames(Y_time_myok),
  loading = loadY_time_myok
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice_head(n = 30)

load_df_myok_time <- dplyr::bind_rows(
  topX_tab_myok_time |> dplyr::mutate(type = "Myokines"),
  topY_tab_myok_time |> dplyr::mutate(type = "Time")
)

(fig_myok_time_axis_01 <- ggplot2::ggplot(
  load_df_myok_time,
  ggplot2::aes(x = loading, y = reorder(var, loading), fill = type)
) +
    ggplot2::geom_col(alpha = 0.85) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~type, scales = "free_y") +
    ggplot2::labs(
      x = "Loading (Canonical component 1)",
      y = NULL,
      title = "Opposing poles of the latent myokines–time axis (M0)",
      subtitle = "Variables on opposite sides represent contrasting biological–behavioral patterns"
    ) +
    ggplot2::theme_minimal(base_size = 13))

save_tiff(
  fig_myok_time_axis_01,
  file.path("output/figures", paste0(format(Sys.Date(), "%y%m%d"), "_myokines_time_axis_01.tiff")),
  res = 150
)

# **********************************************************************
## Export tabulek pro grantový text ----
# **********************************************************************

# 1) MYOKINY × AKTIVITA
tabs_myok_activ <- build_rcca_tables(
  fit = rcc_fit_myok_activ,
  X   = X_myok_activ,
  Y   = Y_activ_myok,
  meta = meta_myok_activ,
  analysis_name = "myokines_vs_activity_M0",
  block_x = "Myokines (M0)",
  block_y = "Activity (M0)",
  ncomp = 2,
  top_n_load_x = 25,
  top_n_load_y = 25,
  top_n_pairs_x = 15,
  top_n_pairs_y = 20,
  axis_cor_thr = 0.15
)

export_rcca_tables(
  tabs_myok_activ,
  file.path("output/tables", paste0(format(Sys.Date(), "%y%m%d"), "_rcca_myokines_vs_activity_M0.xlsx"))
)

# 2) MYOKINY × ČAS
tabs_myok_time <- build_rcca_tables(
  fit = rcc_fit_myok_time,
  X   = X_myok_time,
  Y   = Y_time_myok,
  meta = meta_myok_time,
  analysis_name = "myokines_vs_time_M0",
  block_x = "Myokines (M0)",
  block_y = "Time (M0)",
  ncomp = 2,
  top_n_load_x = 25,
  top_n_load_y = 25,
  top_n_pairs_x = 15,
  top_n_pairs_y = 20,
  axis_cor_thr = 0.15
)

export_rcca_tables(
  tabs_myok_time,
  file.path("output/tables", paste0(format(Sys.Date(), "%y%m%d"), "_rcca_myokines_vs_time_M0.xlsx"))
)




# **********************************************************************
# 5. DIABLO ----
# **********************************************************************
## Data preparation ----
# **********************************************************************

# AKTIVITA
# outcome (Y)
##subtype (factor)
Y_subtype_activ <- d12_join |>
  dplyr::pull(disease_subtype) |>
  as.factor()
## response (factor)
Y_response_activ <- d12_join |>
  dplyr::pull(response_m0_m6) |>
  as.factor()

# X blocks: activity numeric; clinical = model.matrix one-hot
x_cols_diablo_activ <- d10_aktivita |>
  dplyr::select(where(is.numeric)) |>
  names()

X_act_activ <- d12_join |>
  dplyr::select(all_of(x_cols_diablo_activ)) |>
  as.matrix()

X_cli_df_activ <- d12_join |>
  dplyr::select(-c(immet_id, exam_order), -all_of(x_cols_diablo_activ), 
                -disease_subtype, -response_m0_m6)

X_cli_activ <- model.matrix(~ . - 1, data = X_cli_df_activ)

# drop constant cols + scale each block
X_act_activ <- drop_const_cols(X_act_activ)
X_cli_activ <- drop_const_cols(X_cli_activ)

X_act_activ <- scale(X_act_activ)
X_cli_activ <- scale(X_cli_activ)

# assemble blocks for DIABLO
X_blocks_activ <- list(
  activity = X_act_activ,
  clinical = X_cli_activ
)

# design matrix (strength of association between blocks)
design_activ <- matrix(
  c(0, 1,
    1, 0),
  nrow = 2, byrow = TRUE,
  dimnames = list(names(X_blocks_activ), names(X_blocks_activ))
)

# choose keepX (feature selection per comp, per block)
keepX_activ <- list(
  activity = c(5, 5),
  clinical = c(15, 15)
)

# ČAS
# outcome (Y) 
## subtype (factor)
Y_subtype_time <- d13_join |>
  dplyr::pull(disease_subtype) |>
  as.factor()
## response (factor)
Y_response_time <- d13_join |>
  dplyr::pull(response_m0_m6) |>
  as.factor()


# X blocks: time numeric; clinical = model.matrix one-hot
x_cols_diablo_time <- d10_cas_b |>
  dplyr::select(where(is.numeric)) |>
  names()

X_time_time <- d13_join |>
  dplyr::select(all_of(x_cols_diablo_time)) |>
  as.matrix()

X_cli_df_time <- d13_join |>
  dplyr::select(-c(immet_id, exam_order), -all_of(x_cols_diablo_time), 
                -disease_subtype, -response_m0_m6)

X_cli_time <- model.matrix(~ . - 1, data = X_cli_df_time)

# drop constant cols + scale each block
X_time_time <- drop_const_cols(X_time_time)
X_cli_time  <- drop_const_cols(X_cli_time)

X_time_time <- scale(X_time_time)
X_cli_time  <- scale(X_cli_time)

# assemble blocks for DIABLO
X_blocks_time <- list(
  time     = X_time_time,
  clinical = X_cli_time
)

# design matrix (strength of association between blocks)
design_time <- matrix(
  c(0, 1,
    1, 0),
  nrow = 2, byrow = TRUE,
  dimnames = list(names(X_blocks_time), names(X_blocks_time))
)

# choose keepX (feature selection per comp, per block)
keepX_time <- list(
  time     = c(5, 5),
  clinical = c(15, 15)
)

# MYOKINES
# M0-only klinika
d12_m0 <- d12_join |>
  dplyr::filter(tolower(as.character(exam_order)) == "m0") |>
  dplyr::distinct(immet_id, .keep_all = TRUE)

# společné ID (zarovnání X a kliniky)
common_ids_myok <- intersect(d17_myokines$immet_id, d12_m0$immet_id)

# outcome (Y)
## subtype (factor)
Y_subtype_myok <- d12_m0 |>
  dplyr::filter(immet_id %in% common_ids_myok) |>
  dplyr::arrange(match(immet_id, common_ids_myok)) |>
  dplyr::pull(disease_subtype) |>
  as.factor()

## response (factor)
Y_response_myok <- d12_m0 |>
  dplyr::filter(immet_id %in% common_ids_myok) |>
  dplyr::arrange(match(immet_id, common_ids_myok)) |>
  dplyr::pull(response_m0_m6) |>
  as.factor()

# X blocks:
# myokines = numeric (delta expression), clinical = model.matrix one-hot (M0)
x_cols_diablo_myok <- d17_myokines |>
  dplyr::select(-immet_id) |>
  names()

X_myok_myok <- d17_myokines |>
  dplyr::filter(immet_id %in% common_ids_myok) |>
  dplyr::arrange(match(immet_id, common_ids_myok)) |>
  dplyr::select(all_of(x_cols_diablo_myok)) |>
  as.matrix()

X_cli_df_myok <- d12_m0 |>
  dplyr::filter(immet_id %in% common_ids_myok) |>
  dplyr::arrange(match(immet_id, common_ids_myok)) |>
  dplyr::select(
    -c(immet_id, exam_order),
    -disease_subtype, -response_m0_m6
  )

X_cli_myok <- model.matrix(~ . - 1, data = X_cli_df_myok)

# drop constant cols + scale each block
X_myok_myok <- drop_const_cols(X_myok_myok)
X_cli_myok  <- drop_const_cols(X_cli_myok)

X_myok_myok <- scale(X_myok_myok)
X_cli_myok  <- scale(X_cli_myok)

# assemble blocks for DIABLO
X_blocks_myok <- list(
  myokines = X_myok_myok,
  clinical = X_cli_myok
)

# design matrix
design_myok <- matrix(
  c(0, 1,
    1, 0),
  nrow = 2, byrow = TRUE,
  dimnames = list(names(X_blocks_myok), names(X_blocks_myok))
)

# choose keepX
keepX_myok <- list(
  myokines = c(5, 5),
  clinical = c(15, 15)
)


# **********************************************************************
## Analyses ----
# **********************************************************************
### subtype ----
# **********************************************************************
# activity vs clinical (activity-specific)
set.seed(123)
diablo_fit_activ <- mixOmics::block.splsda(
  X = X_blocks_activ,
  Y = Y_subtype_activ,
  ncomp = 2,
  keepX = keepX_activ,
  near.zero.var = TRUE,
  design = design_activ
)
# basic outputs
diablo_fit_activ$names
diablo_fit_activ$variates$activity[, 1]
diablo_fit_activ$variates$clinical[, 1]
# selected variables (signatures)
mixOmics::selectVar(diablo_fit_activ, comp = 1)$activity$name
mixOmics::selectVar(diablo_fit_activ, comp = 1)$clinical$name
# quick performance (CV)
set.seed(123)
perf_diablo_activ <- mixOmics::perf(
  diablo_fit_activ,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  progressBar = TRUE
)
perf_diablo_activ$error.rate

# čas vs clinical (activity-specific)
set.seed(123)
diablo_fit_time <- mixOmics::block.splsda(
  X = X_blocks_time,
  Y = Y_subtype_time,
  ncomp = 2,
  keepX = keepX_time,
  near.zero.var = TRUE,
  design = design_time
)
# basic outputs
diablo_fit_time$names
diablo_fit_time$variates$time[, 1]
diablo_fit_time$variates$clinical[, 1]
# selected variables (signatures)
mixOmics::selectVar(diablo_fit_time, comp = 1)$time$name
mixOmics::selectVar(diablo_fit_time, comp = 1)$clinical$name

# quick performance (CV)
set.seed(123)
perf_diablo_time <- mixOmics::perf(
  diablo_fit_time,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  progressBar = TRUE
)
perf_diablo_time$error.rate

# MYOKINES 
set.seed(123)
diablo_fit_myok <- mixOmics::block.splsda(
  X = X_blocks_myok,
  Y = Y_subtype_myok,
  ncomp = 2,
  keepX = keepX_myok,
  near.zero.var = TRUE,
  design = design_myok
)

# basic outputs
diablo_fit_myok$names
diablo_fit_myok$variates$myokines[, 1]
diablo_fit_myok$variates$clinical[, 1]

# selected variables (signatures)
mixOmics::selectVar(diablo_fit_myok, comp = 1)$myokines$name
mixOmics::selectVar(diablo_fit_myok, comp = 1)$clinical$name

# quick performance (CV)
set.seed(123)
perf_diablo_myok <- mixOmics::perf(
  diablo_fit_myok,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  progressBar = TRUE
)
perf_diablo_myok$error.rate

# **********************************************************************
### response ----
# **********************************************************************
# activity vs clinical (activity-specific)
set.seed(123)
diablo_fit_response_activ <- mixOmics::block.splsda(
  X = X_blocks_activ,
  Y = Y_response_activ,
  ncomp = 2,
  keepX = keepX_activ,
  near.zero.var = TRUE,
  design = design_activ
)
# basic outputs
diablo_fit_response_activ$names
diablo_fit_response_activ$variates$activity[, 1]
diablo_fit_response_activ$variates$clinical[, 1]
# selected variables (signatures)
mixOmics::selectVar(diablo_fit_response_activ, comp = 1)$activity$name
mixOmics::selectVar(diablo_fit_response_activ, comp = 1)$clinical$name
# quick performance (CV)
set.seed(123)
perf_diablo_response_activ <- mixOmics::perf(
  diablo_fit_response_activ,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  progressBar = TRUE
)
perf_diablo_response_activ$error.rate

# čas vs clinical (activity-specific)
set.seed(123)
diablo_fit_response_time <- mixOmics::block.splsda(
  X = X_blocks_time,
  Y = Y_response_time,
  ncomp = 2,
  keepX = keepX_time,
  near.zero.var = TRUE,
  design = design_time
)
# basic outputs
diablo_fit_response_time$names
diablo_fit_response_time$variates$time[, 1]
diablo_fit_response_time$variates$clinical[, 1]
# selected variables (signatures)
mixOmics::selectVar(diablo_fit_response_time, comp = 1)$time$name
mixOmics::selectVar(diablo_fit_response_time, comp = 1)$clinical$name

# quick performance (CV)
set.seed(123)
perf_diablo_response_time <- mixOmics::perf(
  diablo_fit_response_time,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  progressBar = TRUE
)
perf_diablo_response_time$error.rate

# MYOKINES
set.seed(123)
diablo_fit_response_myok <- mixOmics::block.splsda(
  X = X_blocks_myok,
  Y = Y_response_myok,
  ncomp = 2,
  keepX = keepX_myok,
  near.zero.var = TRUE,
  design = design_myok
)

# basic outputs
diablo_fit_response_myok$names
diablo_fit_response_myok$variates$myokines[, 1]
diablo_fit_response_myok$variates$clinical[, 1]

# selected variables (signatures)
mixOmics::selectVar(diablo_fit_response_myok, comp = 1)$myokines$name
mixOmics::selectVar(diablo_fit_response_myok, comp = 1)$clinical$name

# quick performance (CV)
set.seed(123)
perf_diablo_response_myok <- mixOmics::perf(
  diablo_fit_response_myok,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  progressBar = TRUE
)
perf_diablo_response_myok$error.rate

# **********************************************************************
## AKTIVITA outputs ----
# **********************************************************************
###  clustering/separation ----
# **********************************************************************
### subtype ----
# **********************************************************************
(diablo_indiv_activ_01 <- mixOmics::plotIndiv(
  diablo_fit_activ,
  ellipse = FALSE,
  ind.names = FALSE,
  legend = TRUE,
  title = "DIABLO (ACTIVITA): sample plot"
))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_diablo_indiv_activ_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(diablo_indiv_activ_01$graph)
dev.off()

# **********************************************************************
### response ----
# **********************************************************************
(diablo_indiv_activ_response_01 <- mixOmics::plotIndiv(
  diablo_fit_response_activ,
  ellipse = FALSE,
  ind.names = FALSE,
  legend = TRUE,
  title = "DIABLO (ACTIVITA): sample plot"
))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_diablo_indiv_activ_response_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(diablo_indiv_activ_response_01$graph)
dev.off()

# **********************************************************************
### correlation circle per block ----
# **********************************************************************
tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_diablo_plotVar_activ_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
mixOmics::plotVar(
  diablo_fit_activ,
  block = "activity",
  comp = c(1, 2),
  title = "DIABLO (ACTIVITA): activity variables"
)
dev.off()

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_diablo_plotVar_activ_02.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
mixOmics::plotVar(
  diablo_fit_activ,
  block = "clinical",
  comp = c(1, 2),
  title = "DIABLO (ACTIVITA): clinical variables"
)
dev.off()

# **********************************************************************
### heatmap of selected features ----
# **********************************************************************
# x11()
# diablo_cim_activ_01 <- mixOmics::cimDiablo(
#   diablo_fit_activ,
#   comp = 1,
#   legend = TRUE,
#   title = "DIABLO (ACTIVITA): selected features (comp1)"
# )


# x11()
# diablo_cim_activ_01 <- mixOmics::cimDiablo(
#   diablo_fit_response_activ,
#   comp = 1,
#   legend = TRUE,
#   title = "DIABLO (ACTIVITA): selected features (comp1)"
# )

# **********************************************************************
## ČAS outputs ----
# **********************************************************************
### clustering/separation by subtype ----
# **********************************************************************
(diablo_indiv_time_01 <- mixOmics::plotIndiv(
  diablo_fit_time,
  ellipse = FALSE,
  ind.names = FALSE,
  legend = TRUE,
  title = "DIABLO (ČAS): sample plot"
))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_diablo_indiv_time_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(diablo_indiv_time_01$graph)
dev.off()

# **********************************************************************
### correlation circle per block ----
# **********************************************************************
tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_diablo_plotVar_time_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
mixOmics::plotVar(
  diablo_fit_time,
  block = "time",
  comp = c(1, 2),
  title = "DIABLO (ČAS): time variables"
)
dev.off()

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_diablo_plotVar_time_02.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
mixOmics::plotVar(
  diablo_fit_time,
  block = "clinical",
  comp = c(1, 2),
  title = "DIABLO (ČAS): clinical variables"
)
dev.off()

# **********************************************************************
### heatmap of selected features ----
# **********************************************************************
# x11()
# diablo_cim_time_01 <- mixOmics::cimDiablo(
#   diablo_fit_time,
#   comp = 1,
#   legend = TRUE,
#   title = "DIABLO (ČAS): selected features (comp1)"
# )

# **********************************************************************


# **********************************************************************
## MYOKINY (M0 only) ----
# **********************************************************************
## Data preparation ----
# **********************************************************************



# **********************************************************************
## Analyses ----
# **********************************************************************
### subtype ----
# **********************************************************************

# **********************************************************************
### response ----
# **********************************************************************
set.seed(123)
diablo_fit_response_myok <- mixOmics::block.splsda(
  X = X_blocks_myok,
  Y = Y_response_myok,
  ncomp = 2,
  keepX = keepX_myok,
  near.zero.var = TRUE,
  design = design_myok
)

# basic outputs
diablo_fit_response_myok$names
diablo_fit_response_myok$variates$myokines[, 1]
diablo_fit_response_myok$variates$clinical[, 1]

# selected variables (signatures)
mixOmics::selectVar(diablo_fit_response_myok, comp = 1)$myokines$name
mixOmics::selectVar(diablo_fit_response_myok, comp = 1)$clinical$name

# quick performance (CV)
set.seed(123)
perf_diablo_response_myok <- mixOmics::perf(
  diablo_fit_response_myok,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  progressBar = TRUE
)
perf_diablo_response_myok$error.rate


# **********************************************************************
## MYOKINY outputs ----
# **********************************************************************
###  clustering/separation ----
# **********************************************************************
### subtype ----
# **********************************************************************
(diablo_indiv_myok_01 <- mixOmics::plotIndiv(
  diablo_fit_myok,
  ellipse = FALSE,
  ind.names = FALSE,
  legend = TRUE,
  title = "DIABLO (MYOKINY, M0): sample plot"
))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_diablo_indiv_myok_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(diablo_indiv_myok_01$graph)
dev.off()

# **********************************************************************
### response ----
# **********************************************************************
(diablo_indiv_myok_response_01 <- mixOmics::plotIndiv(
  diablo_fit_response_myok,
  ellipse = FALSE,
  ind.names = FALSE,
  legend = TRUE,
  title = "DIABLO (MYOKINY, M0): sample plot"
))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_diablo_indiv_myok_response_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(diablo_indiv_myok_response_01$graph)
dev.off()

# **********************************************************************
### correlation circle per block ----
# **********************************************************************
tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_diablo_plotVar_myok_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
mixOmics::plotVar(
  diablo_fit_myok,
  block = "myokines",
  comp = c(1, 2),
  title = "DIABLO (MYOKINY, M0): myokines variables"
)
dev.off()

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_diablo_plotVar_myok_02.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
mixOmics::plotVar(
  diablo_fit_myok,
  block = "clinical",
  comp = c(1, 2),
  title = "DIABLO (MYOKINY, M0): clinical variables"
)
dev.off()

# **********************************************************************
### heatmap of selected features ----
# **********************************************************************
# x11()
# diablo_cim_myok_01 <- mixOmics::cimDiablo(
#   diablo_fit_myok,
#   comp = 1,
#   legend = TRUE,
#   title = "DIABLO (MYOKINY, M0): selected features (comp1)"
# )

# x11()
# diablo_cim_myok_response_01 <- mixOmics::cimDiablo(
#   diablo_fit_response_myok,
#   comp = 1,
#   legend = TRUE,
#   title = "DIABLO (MYOKINY, M0): selected features (comp1)"
# )

# **********************************************************************
# 6. MOFA ----
# **********************************************************************
## Data preparation ----
# **********************************************************************

# ACTIVITA
# define "Activity" columns (numeric columns from the ORIGINAL activity table)
mofa_x_cols_activ <- d10_aktivita |>
  dplyr::select(where(is.numeric)) |>
  names()

# build Activity view matrix: features x samples
mofa_X_activ <- d12_mofa |>
  dplyr::select(all_of(mofa_x_cols_activ)) |>
  as.matrix()

rownames(mofa_X_activ) <- d12_mofa$sample
mofa_X_activ <- t(mofa_X_activ)                   # features x samples
colnames(mofa_X_activ) <- d12_mofa$sample          # explicit (safety)

# build Clinical view matrix: one-hot encoded, features x samples
mofa_Y_df_activ <- d12_mofa |>
  dplyr::select(-c(immet_id, exam_order, sample, subtype), -all_of(mofa_x_cols_activ))

mofa_Y_mm_activ <- model.matrix(~ . - 1, data = mofa_Y_df_activ)

rownames(mofa_Y_mm_activ) <- d12_mofa$sample
mofa_Y_activ <- t(mofa_Y_mm_activ)                # features x samples
colnames(mofa_Y_activ) <- d12_mofa$sample

# drop constant features (IMPORTANT; avoids errors in training)
mofa_X_activ <- drop_const_features(mofa_X_activ)
mofa_Y_activ <- drop_const_features(mofa_Y_activ)

# create MOFA object
mofa_views_activ <- list(
  Activity = mofa_X_activ,
  Clinical = mofa_Y_activ
)

mofa_obj_activ <- create_mofa(mofa_views_activ)

# attach metadata (MUST contain column named 'sample')
mofa_meta_activ <- d12_mofa |>
  dplyr::transmute(
    sample     = sample,       # REQUIRED by MOFA2
    immet_id   = immet_id,
    exam_order = exam_order,
    subtype    = subtype
  )

samples_metadata(mofa_obj_activ) <- mofa_meta_activ

# configure model
data_opts_activ  <- get_default_data_options(mofa_obj_activ)
model_opts_activ <- get_default_model_options(mofa_obj_activ)
train_opts_activ <- get_default_training_options(mofa_obj_activ)

model_opts_activ$num_factors <- 10
train_opts_activ$maxiter <- 1000

mofa_obj_activ <- prepare_mofa(
  mofa_obj_activ,
  data_options     = data_opts_activ,
  model_options    = model_opts_activ,
  training_options = train_opts_activ
)

# ČAS
# define "Time" columns (numeric columns from the ORIGINAL time table)
mofa_x_cols_time <- d10_cas_b |>
  dplyr::select(where(is.numeric)) |>
  names()

# build Time view matrix: features x samples
mofa_X_time <- d13_mofa |>
  dplyr::select(all_of(mofa_x_cols_time)) |>
  as.matrix()

rownames(mofa_X_time) <- d13_mofa$sample
mofa_X_time <- t(mofa_X_time)                    # features x samples
colnames(mofa_X_time) <- d13_mofa$sample          # explicit (safety)

# build Clinical view matrix: one-hot encoded, features x samples
mofa_Y_df_time <- d13_mofa |>
  dplyr::select(-c(immet_id, exam_order, sample, subtype), -all_of(mofa_x_cols_time))

mofa_Y_mm_time <- model.matrix(~ . - 1, data = mofa_Y_df_time)

rownames(mofa_Y_mm_time) <- d13_mofa$sample
mofa_Y_time <- t(mofa_Y_mm_time)                 # features x samples
colnames(mofa_Y_time) <- d13_mofa$sample

# drop constant features
mofa_X_time <- drop_const_features(mofa_X_time)
mofa_Y_time <- drop_const_features(mofa_Y_time)

# create MOFA object
mofa_views_time <- list(
  Time     = mofa_X_time,
  Clinical = mofa_Y_time
)

mofa_obj_time <- create_mofa(mofa_views_time)

# attach metadata
mofa_meta_time <- d13_mofa |>
  dplyr::transmute(
    sample     = sample,       # REQUIRED by MOFA2
    immet_id   = immet_id,
    exam_order = exam_order,
    subtype    = subtype
  )

samples_metadata(mofa_obj_time) <- mofa_meta_time

# configure model
data_opts_time  <- get_default_data_options(mofa_obj_time)
model_opts_time <- get_default_model_options(mofa_obj_time)
train_opts_time <- get_default_training_options(mofa_obj_time)

model_opts_time$num_factors <- 10
train_opts_time$maxiter <- 1000

mofa_obj_time <- prepare_mofa(
  mofa_obj_time,
  data_options     = data_opts_time,
  model_options    = model_opts_time,
  training_options = train_opts_time
)

# **********************************************************************
## Analyses ----
# **********************************************************************

# train (use_basilisk solves Python env issues)
mofa_trained_activ <- run_mofa(mofa_obj_activ, use_basilisk = TRUE)
mofa_trained_time  <- run_mofa(mofa_obj_time, use_basilisk = TRUE)


# **********************************************************************
## AKTIVITA outputs ----
# **********************************************************************
### variance plot ----
# **********************************************************************
plot_variance_explained(mofa_trained_activ)

# **********************************************************************
### factors ----
# **********************************************************************
(fig_mofa_fact_subt_activ_01 <- plot_factors(
  mofa_trained_activ,
  factors  = 1:2,
  color_by = "subtype"
))
tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_mofa_activ_clin_fact_subt_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_mofa_fact_subt_activ_01)
dev.off()

(fig_mofa_fact_exam_activ_01 <- plot_factors(
  mofa_trained_activ,
  factors  = 1:2,
  color_by = "exam_order"
))
tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_mofa_activ_clin_fact_exam_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_mofa_fact_exam_activ_01)
dev.off()

(fig_mofa_fact_subt_activ_02 <- plot_factor(
  mofa_trained_activ,
  factors     = 1,
  color_by    = "subtype",
  dodge       = TRUE,
  add_violin  = TRUE
))
tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_mofa_activ_clin_fact_subt_02.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_mofa_fact_subt_activ_02)
dev.off()

(fig_mofa_fact_exam_activ_02 <- plot_factor(
  mofa_trained_activ,
  factors     = 1,
  color_by    = "exam_order",
  dodge       = TRUE,
  add_violin  = TRUE
))
tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_mofa_activ_clin_fact_exam_02.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_mofa_fact_exam_activ_02)
dev.off()

(fig_mofa_weight_activ_01 <- plot_top_weights(
  mofa_trained_activ, view = "Activity", factor = 1, nfeatures = 20
))
(fig_mofa_weight_activ_02 <- plot_top_weights(
  mofa_trained_activ, view = "Clinical", factor = 1, nfeatures = 20
))
tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_mofa_activ_clin_weights_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
gridExtra::grid.arrange(fig_mofa_weight_activ_01, fig_mofa_weight_activ_02, nrow = 1)
dev.off()

# **********************************************************************
### UMAP on factors ----
# **********************************************************************
mofa_trained_activ <- run_umap(mofa_trained_activ)
(p1_act <- plot_dimred(mofa_trained_activ, method = "UMAP", color_by = "subtype") +
    ggtitle("Colored by subtype"))
(p2_act <- plot_dimred(mofa_trained_activ, method = "UMAP", color_by = "exam_order") +
    ggtitle("Colored by exam (M0 vs M6)"))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_mofa_activ_clin_umap_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
gridExtra::grid.arrange(p1_act, p2_act, nrow = 1)
dev.off()


# **********************************************************************
## ČAS outputs ----
# **********************************************************************
# **********************************************************************
### variance plot ----
# **********************************************************************
plot_variance_explained(mofa_trained_time)

# **********************************************************************
### factors ----
# **********************************************************************
(fig_mofa_fact_subt_time_01 <- plot_factors(
  mofa_trained_time,
  factors  = 1:2,
  color_by = "subtype"
))
tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_mofa_time_clin_fact_subt_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_mofa_fact_subt_time_01)
dev.off()

(fig_mofa_fact_exam_time_01 <- plot_factors(
  mofa_trained_time,
  factors  = 1:2,
  color_by = "exam_order"
))
tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_mofa_time_clin_fact_exam_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_mofa_fact_exam_time_01)
dev.off()

(fig_mofa_fact_subt_time_02 <- plot_factor(
  mofa_trained_time,
  factors     = 1,
  color_by    = "subtype",
  dodge       = TRUE,
  add_violin  = TRUE
))
tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_mofa_time_clin_fact_subt_02.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_mofa_fact_subt_time_02)
dev.off()

(fig_mofa_fact_exam_time_02 <- plot_factor(
  mofa_trained_time,
  factors     = 1,
  color_by    = "exam_order",
  dodge       = TRUE,
  add_violin  = TRUE
))
tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_mofa_time_clin_fact_exam_02.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
print(fig_mofa_fact_exam_time_02)
dev.off()

(fig_mofa_weight_time_01 <- plot_top_weights(
  mofa_trained_time, view = "Time", factor = 1, nfeatures = 20
))
(fig_mofa_weight_time_02 <- plot_top_weights(
  mofa_trained_time, view = "Clinical", factor = 1, nfeatures = 20
))
tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_mofa_time_clin_weights_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
gridExtra::grid.arrange(fig_mofa_weight_time_01, fig_mofa_weight_time_02, nrow = 1)
dev.off()

# **********************************************************************
### UMAP on factors ----
# **********************************************************************
mofa_trained_time <- run_umap(mofa_trained_time)
(p1_time <- plot_dimred(mofa_trained_time, method = "UMAP", color_by = "subtype") +
    ggtitle("Colored by subtype"))
(p2_time <- plot_dimred(mofa_trained_time, method = "UMAP", color_by = "exam_order") +
    ggtitle("Colored by exam (M0 vs M6)"))

tiff(
  filename = file.path(
    "output/figures/final",
    paste0(format(Sys.Date(), "%y%m%d"), "_mofa_time_clin_umap_01.tiff")
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
gridExtra::grid.arrange(p1_time, p2_time, nrow = 1)
dev.off()

# **********************************************************************
# 7. REPORT EXPORT: objects for <out_date>_ggir_clinics_01 ----
# **********************************************************************
report_rdata <- file.path("reports", paste0(out_date, "_ggir_clinics_01.RData"))

objs_to_save <- c(
  # --- data (M0/M6 context) ---
  "d10_aktivita", "d10_cas_b", "d11_clinics_imp",
  "d12_join", "d13_join", "d12_join_delta_activ", "d13_join_delta_time",
  "d12_mofa", "d13_mofa",
  # --- correlation outputs ---
  "x_cols_corr_activ", "y_cols_corr_activ",
  "x_cols_corr_time", "y_cols_corr_time",
  "cor_df_all_activ", "cor_df_m0_activ", "cor_df_m6_activ", "cor_df_delta_activ",
  "cor_df_all_time", "cor_df_m0_time", "cor_df_m6_time", "cor_df_delta_time",
  # --- rCCA ---
  "rcc_fit_activ", "rcc_fit_time",
  "scores_df_activ", "scores_df_time",
  "topX_activ", "topY_activ", "topX_time", "topY_time",
  "fig_can_comp_activ_01", "fig_heatmap_loading_activ_01", "fig_rcca_subtype_activ_01",
  "fig_can_comp_time_01", "fig_heatmap_loading_time_01", "fig_rcca_subtype_time_01",
  # --- DIABLO ---
  "diablo_fit_activ", "diablo_fit_time",
  "perf_diablo_activ", "perf_diablo_time",
  "design_activ", "design_time", "keepX_activ", "keepX_time",
  "diablo_indiv_activ_01", "diablo_indiv_time_01",
  # --- MOFA ---
  "mofa_obj_activ", "mofa_obj_time",
  "mofa_trained_activ", "mofa_trained_time",
  "mofa_meta_activ", "mofa_meta_time",
  "fig_mofa_fact_subt_activ_01", "fig_mofa_fact_exam_activ_01",
  "fig_mofa_fact_subt_time_01", "fig_mofa_fact_exam_time_01"
)

objs_to_save <- objs_to_save[objs_to_save %in% ls()]

save(list = objs_to_save, file = report_rdata)
message("Saved report objects to: ", normalizePath(report_rdata))




})

local({
# **********************************************************************
# 1. preparation ----
# **********************************************************************
# **********************************************************************
## libraries ----
# **********************************************************************
library(lme4)
library(broom.mixed)
library(recipes)


# **********************************************************************
## data ----
# **********************************************************************
d14_join <- d10_aktivita |>
  transmute(
    immet_id   = as.character(sample),
    exam_order = tolower(as.character(exam)),
    across(where(is.numeric), identity)
  ) |>
  inner_join(
    d10_cas_b |>
      transmute(
        immet_id   = as.character(sample),
        exam_order = tolower(as.character(exam)),
        across(where(is.numeric), identity)
      ),
    by = c("immet_id", "exam_order"),
    suffix = c("_activ", "_time")
  ) |>
  inner_join(
    d11_clinics_imp |>
      mutate(exam_order = tolower(as.character(exam_order))) |>
      mutate(immet_id = as.character(immet_id)),
    by = c("immet_id", "exam_order")
  ) |>
  filter(exam_order %in% c("m0", "m6")) |>
  distinct()

age_m0 <- d14_join |>
  filter(exam_order == "m0") |>
  select(immet_id, age) |>
  rename(age_m0 = age)

d14_join <- d14_join |>
  left_join(age_m0, by = "immet_id")

# **********************************************************************
# 2. Analyses ----
# **********************************************************************
## Clinics effect on accelerometry ----
# **********************************************************************
### Latent score ----
# **********************************************************************

x_activ <- c(
  "enmo_full_recording_mean","acc_day_spt_wei","m5_enmo","m1_enmo",
  "p9167_enmo","p9375_enmo","p9583_enmo","p9688_enmo","p9792_enmo",
  "p9861_enmo","p9931_enmo","ig_gradient_enmo","ig_intercept_enmo"
)

x_time <- c(
  "m5hr_time","t20_time_min","t20_b1m80_time_min","t20_b5m80_time_min","t20_b10m80_time_min",
  "t30_time_min","t30_b1m80_time_min","t30_b5m80_time_min","t30_b10m80_time_min",
  "t40_time_min","t40_b1m80_time_min","t40_b5m80_time_min","t40_b10m80_time_min",
  "t100_time_min","t100_b1m80_time_min","t100_b5m80_time_min","t100_b10m80_time_min",
  "moderate", "light", "inactivity"
)

prep_scores <- function(df, cols, prefix){
  rec <- recipe(~ ., data = df |> select(all_of(cols))) |>
    step_impute_knn(all_predictors()) |>
    step_zv(all_predictors()) |>
    step_corr(all_predictors(), threshold = 0.9) |>
    step_normalize(all_predictors()) |>
    step_pca(all_predictors(), num_comp = 2) |>
    prep()
  
  scores <- bake(rec, new_data = NULL) |>
    transmute(
      !!paste0(prefix, "_pc1") := PC1,
      !!paste0(prefix, "_pc2") := PC2
    )
  list(recipe = rec, scores = scores)
}

tmp_activ <- prep_scores(d14_join, x_activ, "activ")
tmp_time  <- prep_scores(d14_join, x_time,  "time")

d15_scores <- d14_join |>
  bind_cols(tmp_activ$scores, tmp_time$scores)


clin_num <- clin_cols |> keep(~ is.numeric(d11_clinics_imp[[.x]]))
clin_fac <- clin_cols |> keep(~ is.factor(d11_clinics_imp[[.x]]) || is.character(d11_clinics_imp[[.x]]) || is.logical(d11_clinics_imp[[.x]]))

rec_scaled <- recipe(~ ., data = d15_scores) |>
  step_impute_knn(all_of(clin_num), "age_m0") |>
  step_normalize(all_of(clin_num), "age_m0") |>
  prep()

d15_scaled <- bake(rec_scaled, new_data = d15_scores)

# příklad: klinická proměnná physician_vas, adjustace sex + bmi + age_m0
# fit1 <- lmer(
#   activ_pc1 ~ exam_order + physician_vas + exam_order:physician_vas +
#     sex + bmi + age_m0 + (1 | immet_id),
#   data = d15_scaled
# )
# 
# tidy(fit1, effects = "fixed", conf.int = TRUE)

pc_cols <- c("activ_pc1", "time_pc1")
clin_exclude <- c("immet_id", "exam_order")  
clin_cols <- setdiff(names(d11_clinics_imp), c(clin_exclude, pc_cols))
base_covars <- c("bmi", "age_m0") |> intersect(names(d15_scores))

d15_scaled <- bake(rec_scaled, new_data = d15_scores) |> 
  mutate(exam_order = factor(tolower(as.character(exam_order)), levels = c("m0", "m6"))) |>
  mutate(immet_id = as.character(immet_id))

pc_long <- d15_scaled |>
  select(immet_id, exam_order, all_of(pc_cols)) |>
  pivot_longer(cols = all_of(pc_cols), names_to = "pc_name", values_to = "pc_value")


clin_long_num <- d15_scaled |>
  select(immet_id, exam_order, disease_subtype, all_of(clin_num), all_of(base_covars)) |>
  pivot_longer(cols = setdiff(clin_num, base_covars), 
               names_to = "clin_name", values_to = "clin_value") |>
  mutate(clin_type = "numeric",
         clin_value = as.numeric(clin_value)) |> 
  full_join(pc_long, by = c("immet_id", "exam_order")) |> 
  group_by(clin_name, pc_name) |> 
  nest() |> 
  ungroup()

clin_long_fac <- d15_scaled |>
  select(immet_id, exam_order, disease_subtype, response_m0_m6, all_of(clin_fac), all_of(base_covars)) |>
  mutate(
    response_m0_m6   = factor(response_m0_m6),
    disease_subtype  = factor(disease_subtype)
  ) |>
  pivot_longer(
    cols = -c(immet_id, exam_order, disease_subtype, response_m0_m6, all_of(base_covars)),
    names_to = "clin_name",
    values_to = "clin_value"
  ) |>
  mutate(
    clin_type  = "factor",
    clin_value = factor(clin_value)
  ) |>
  full_join(pc_long, by = c("immet_id", "exam_order")) |>
  group_by(clin_name, pc_name) |>
  nest() |>
  ungroup()


clin_long <- bind_rows(clin_long_num, clin_long_fac)


res_clin_act_time_01 <- clin_long |> 
  mutate(mod = map(data, ~possLMER(data = .x,
                               pc_value ~ exam_order*clin_value + bmi + age_m0 + 
                                 (1 | immet_id))),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE)))

res_clin_act_time_01_tab <- res_clin_act_time_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "clin_value")) |> 
  mutate(term = str_replace(term, "clin_value", clin_name),
         pc_name = str_remove(pc_name, "_pc1")) |> 
  select(-data, -mod, -x, -std.error, -statistic, -df, -clin_name, -effect)

 
res_clin_act_sel_01 <- res_clin_act_time_01_tab |> 
  filter(p.value < 0.06 & pc_name == "activ") |> 
  pull(term) |> 
  str_remove("exam_orderm6:") |> 
  unique() |> 
  str_replace("response_m0_m6Minimal", "response_m0_m6") |> 
  str_replace("mtx1", "mtx") |> 
  str_replace("cpa1", "cpa") |> 
  str_replace("arthritis1", "arthritis") |> 
  str_replace("raynaud1", "raynaud") |> 
  str_replace("oesophageal_motility_disorder1", "oesophageal_motility_disorder")

res_clin_time_sel_01 <- res_clin_act_time_01_tab |> 
  filter(p.value < 0.06 & pc_name == "time") |> 
  pull(term) |> 
  str_remove("exam_orderm6:") |> 
  unique() |> 
  str_replace("disease_subtype2", "disease_subtype") |> 
  str_replace("mtx1", "mtx") |> 
  str_replace("cpa1", "cpa") |> 
  str_replace("arthritis1", "arthritis") |> 
  str_replace("raynaud1", "raynaud") |> 
  str_replace("ild1", "ild") |> 
  str_replace("statins1", "statins")


# **********************************************************************
## Confirmatory analyses ----
# **********************************************************************
### LMER ----
# **********************************************************************

clin_long_num_02 <- d15_scaled |> 
  select(immet_id, exam_order, all_of(x_activ), all_of(x_time), all_of(base_covars),
         intersect(res_clin_act_sel_01, clin_num)) |> 
  pivot_longer(cols = intersect(res_clin_act_sel_01, clin_num), 
               names_to = "clin_name", values_to = "clin_value") |> 
  pivot_longer(cols = all_of(c(x_activ, x_time)), 
               names_to = "var_name", values_to = "var_value") |> 
  mutate(clin_type = "numeric",
         clin_value = as.numeric(clin_value)) |> 
  group_by(clin_name, var_name) |> 
  nest() |> 
  ungroup()


clin_long_fac_02 <- d15_scaled |> 
  select(immet_id, exam_order, all_of(x_activ), all_of(x_time), all_of(base_covars),
         intersect(res_clin_act_sel_01, clin_fac)) |> 
  pivot_longer(cols = intersect(res_clin_act_sel_01, clin_fac), 
               names_to = "clin_name", values_to = "clin_value") |> 
  pivot_longer(cols = all_of(c(x_activ, x_time)), 
               names_to = "var_name", values_to = "var_value") |> 
  mutate(clin_type = "factor",
         clin_value = as.factor(clin_value)) |> 
  group_by(clin_name, var_name) |> 
  nest() |> 
  ungroup()

clin_long_02 <- bind_rows(clin_long_num_02, clin_long_fac_02)


res_clin_act_time_02 <- clin_long_02 |> 
  mutate(mod = map(data, ~possLMER(data = .x,
                                   var_value ~ exam_order*clin_value + bmi + age_m0 + 
                                     (1 | immet_id))),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE)))


res_clin_act_time_02_tab <- res_clin_act_time_02 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "clin_value")) |> 
  mutate(term = str_replace(term, "clin_value", clin_name),
         term = str_replace(term, "exam_orderm6", "exam")) |> 
  select(-data, -mod, -std.error, -statistic, -df,  -effect) |> 
  mutate(ggir_type = case_when(
    var_name %in% x_activ ~ "active",
    var_name %in% x_time ~ "time"
  ))

res_clin_act_time_02_tab_sig <- res_clin_act_time_02_tab |>
  dplyr::group_by(var_name, clin_name) |>
  dplyr::filter(
    any(
      (p.value < 0.06) & (sign(conf.low) == sign(conf.high)),
      na.rm = TRUE)) |>
  dplyr::ungroup()

export(res_clin_act_time_02_tab_sig,
       file.path("output/tables/final", paste0(out_date, "_assoc_clin_activ_time_signif_01.xlsx")))

plot_emm_clin_exam <- function(model, clin_name, var_name) {
  sjPlot::plot_model(
    model = model,
    type  = "emm",
    term  = c("clin_value", "exam_order")
  ) +
    sjPlot::theme_sjplot2() +
    ggplot2::labs(
      title    = "Association between variables",
      x = clin_name,
      y = var_name
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
      axis.title    = ggplot2::element_text(size = 12, face = "bold"),
      axis.text     = ggplot2::element_text(size = 11, face = "bold")
    )
}


res_clin_act_time_02_fig <- res_clin_act_time_02 |>
  dplyr::semi_join(
    res_clin_act_time_02_tab_sig |>
      dplyr::select(clin_name, var_name) |>
      dplyr::distinct(),
    by = c("clin_name", "var_name")) |> 
  mutate(fig = pmap(list(mod, clin_name, var_name),
                    ~plot_emm_clin_exam(model = ..1, clin_name = ..2, var_name = ..3)))


pdf(file.path("output/figures", 
              paste0(format(Sys.Date(), "%y%m%d"),
              "_", 
              "signif_assoc_clin_time_activ_01.pdf")), 
    width = 10, height = 7)
walk(res_clin_act_time_02_fig$fig, print)
dev.off()



# **********************************************************************
### sPLS ----
# **********************************************************************

clin_sel <- intersect(res_clin_act_time_02_tab_sig$clin_name |> unique(), 
                      names(d15_scaled))
y_activ  <- intersect(x_activ,  names(d15_scaled))
y_time   <- intersect(x_time,   names(d15_scaled))

df <- d15_scaled


# ******************************
# 1) Helper: prepare X/Y matrices for a given timepoint
# ******************************
prep_xy <- function(df, exam_level, x_vars, y_vars, id_var = "immet_id") {
  
  dat <- df |>
    filter(exam_order == exam_level) |>
    select(all_of(id_var), all_of(x_vars), all_of(y_vars)) |>
    distinct()
  
  # drop columns with 0 variance (mixOmics will complain)
  nzv <- function(v) sd(v, na.rm = TRUE) > 0
  x_keep <- x_vars |> keep(~ is.numeric(dat[[.x]]) && nzv(dat[[.x]]))
  y_keep <- y_vars |> keep(~ is.numeric(dat[[.x]]) && nzv(dat[[.x]]))
  
  # NOTE: sPLS works on numeric matrices; factors need to be encoded beforehand if you want them in X
  # Here we keep only numeric clinical predictors. If you want binary factors, convert 0/1 to numeric in advance.
  X <- dat |> select(all_of(x_keep)) |> as.matrix()
  Y <- dat |> select(all_of(y_keep)) |> as.matrix()
  
  # simple median imputation (recommended)
  for (j in seq_len(ncol(X))) {
    if (anyNA(X[, j])) X[is.na(X[, j]), j] <- median(X[, j], na.rm = TRUE)
  }
  for (j in seq_len(ncol(Y))) {
    if (anyNA(Y[, j])) Y[is.na(Y[, j]), j] <- median(Y[, j], na.rm = TRUE)
  }
  
  list(X = X, Y = Y, x_keep = x_keep, y_keep = y_keep)
}

# ******************************
# 2) Helper: tune sPLS
# ******************************
tune_spls_xy <- function(X, Y, ncomp = 2,
                         test_keepX = NULL,
                         test_keepY = NULL,
                         folds = 5, nrepeat = 30, seed = 1) {
  
  set.seed(seed)
  
  if (is.null(test_keepX)) {
    # if your X is already a shortlist, you can also set keepX = ncol(X) and skip tuning on X
    test_keepX <- unique(pmin(c(3, 5, 8, 10, 15), ncol(X)))
  } else {
    test_keepX <- unique(pmin(test_keepX, ncol(X)))
  }
  
  if (is.null(test_keepY)) {
    test_keepY <- unique(pmin(c(3, 5, 8, 10, 15), ncol(Y)))
  } else {
    test_keepY <- unique(pmin(test_keepY, ncol(Y)))
  }
  
  tune <- tune.spls(
    X, Y,
    ncomp = ncomp,
    test.keepX = test_keepX,
    test.keepY = test_keepY,
    validation = "Mfold",
    folds = folds,
    nrepeat = nrepeat,
    progressBar = TRUE,
    scale = TRUE   # important unless you are 100% sure everything is already comparable
  )
  
  # choice: one keepX/keepY per component
  # tune$choice.keepX and $choice.keepY are vectors of length ncomp
  list(
    tune = tune,
    keepX = tune$choice.keepX,
    keepY = tune$choice.keepY
  )
}

# ******************************
# 3) Fit + summarise
# ******************************
fit_spls_xy <- function(X, Y, keepX, keepY, ncomp = 2) {
  fit <- spls(X, Y, ncomp = ncomp, keepX = keepX, keepY = keepY, scale = TRUE)
  
  # selected variables
  selX <- lapply(seq_len(ncomp), function(h) selectVar(fit, comp = h)$name)
  selY <- lapply(seq_len(ncomp), function(h) selectVar(fit, comp = h, block = "Y")$name) # sometimes not needed
  names(selX) <- paste0("comp", seq_len(ncomp))
  names(selY) <- paste0("comp", seq_len(ncomp))
  
  list(fit = fit, selX = selX, selY = selY)
}

# ******************************
# 4) Run: AKTIVITA and CAS, separately for M0 and M6
# ******************************

run_one <- function(exam_level, y_vars, label) {
  
  xy <- prep_xy(df, exam_level = exam_level, x_vars = clin_sel, y_vars = y_vars)
  
  # if your clin_sel is already small, you can force keepX = ncol(X) (no selection on X)
  # and only tune Y:
  # keepX_fixed <- rep(ncol(xy$X), 2)
  # keepY_tune <- tune_spls_xy(xy$X, xy$Y, ncomp = 2, test_keepX = keepX_fixed)
  
  tuned <- tune_spls_xy(xy$X, xy$Y, ncomp = 2, folds = 5, nrepeat = 30, seed = 1)
  
  fitted <- fit_spls_xy(xy$X, xy$Y, keepX = tuned$keepX, keepY = tuned$keepY, ncomp = 2)
  
  list(
    exam = exam_level,
    block = label,
    X_vars_used = xy$x_keep,
    Y_vars_used = xy$y_keep,
    tune = tuned$tune,
    keepX = tuned$keepX,
    keepY = tuned$keepY,
    fit = fitted$fit,
    selected_X = fitted$selX,
    selected_Y = fitted$selY
  )
}

# --- Activity ---
spls_activ_m0 <- run_one("m0", y_activ, "AKTIVITA")
spls_activ_m6 <- run_one("m6", y_activ, "AKTIVITA")

# --- Time ---
spls_time_m0  <- run_one("m0", y_time, "CAS")
spls_time_m6  <- run_one("m6", y_time, "CAS")

# ******************************
## Quick outputs / plots ----
# ******************************

# Example: correlation of variates (X vs Y) for M0 activity
plotIndiv(spls_activ_m0$fit, comp = c(1,2), group = NULL, legend = TRUE,
          title = "sPLS (M0) - AKTIVITA - X variates")
plotVar(spls_activ_m0$fit, comp = c(1,2), var.names = TRUE,
        title = "sPLS (M0) - AKTIVITA - selected variables")

plotVar(spls_time_m0$fit, comp = c(1,2), var.names = TRUE,
        title = "sPLS (M0) - TIME - selected variables")


# Network plot (often very useful)
network(spls_activ_m0$fit, comp = 1, cutoff = 0.3, title = "sPLS network (M0) - AKTIVITA - comp1")
network(spls_time_m0$fit,  comp = 1, cutoff = 0.3, title = "sPLS network (M0) - CAS - comp1")

# Selected variables (lists)
spls_activ_m0$selected_X
spls_activ_m0$selected_Y

# **********************************************************************
# 3. Accelerometry - subtype ----
# **********************************************************************




# **********************************************************************
## LMER ----
# **********************************************************************
### subtype ----
# **********************************************************************
res_subtype_act_time_01 <- d15_scores |> 
  select(immet_id, exam_order, all_of(pc_cols), disease_subtype, all_of(base_covars)) |>
  pivot_longer(cols = all_of(pc_cols), names_to = "pc_name", values_to = "pc_value") |> 
  group_by(pc_name) |> 
  nest() |> 
  mutate(mod = map(data, ~possLMER(data = .x,
                                   pc_value ~ exam_order + disease_subtype + bmi + age_m0 + 
                                     (1 | immet_id))),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE)))

res_subtype_act_time_01_tab <- res_subtype_act_time_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "disease_subtype")) |> 
  mutate(pc_name = str_remove(pc_name, "_pc1")) |> 
  select(-data, -mod, -std.error, -statistic, -df, -effect)

res_subtype_act_sel_01 <- res_subtype_act_time_01_tab |> 
  filter(p.value < 0.06 & pc_name == "activ") 

res_subtype_time_sel_01 <- res_subtype_act_time_01_tab |> 
  filter(p.value < 0.06 & pc_name == "time")

# **********************************************************************
### response ----
# **********************************************************************
res_response_act_time_01 <- d15_scores |> 
  select(immet_id, exam_order, all_of(pc_cols), response_m0_m6, all_of(base_covars)) |>
  pivot_longer(cols = all_of(pc_cols), names_to = "pc_name", values_to = "pc_value") |> 
  group_by(pc_name) |> 
  nest() |> 
  mutate(mod = map(data, ~possLMER(data = .x,
                                   pc_value ~ exam_order + response_m0_m6 + bmi + age_m0 + 
                                     (1 | immet_id))),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE)))

res_response_act_time_01_tab <- res_response_act_time_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "response")) |> 
  mutate(pc_name = str_remove(pc_name, "_pc1")) |> 
  select(-data, -mod, -std.error, -statistic, -df, -effect)

res_response_act_sel_01 <- res_response_act_time_01_tab |> 
  filter(p.value < 0.06 & pc_name == "activ") 

res_response_time_sel_01 <- res_response_act_time_01_tab |> 
  filter(p.value < 0.06 & pc_name == "time")

# **********************************************************************
## Confirmatory ----
# **********************************************************************
safe_emm   <- possibly(emmeans, otherwise = NULL)
safe_contr <- possibly(contrast, otherwise = NULL)
safe_tidy  <- possibly(broom::tidy, otherwise = tibble())

res_subtype_act_02 <- d15_scores |> 
  select(immet_id, exam_order, all_of(x_activ), all_of(base_covars), disease_subtype) |> 
  pivot_longer(cols = all_of(x_activ), 
               names_to = "var_name", values_to = "var_value") |> 
  group_by(var_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, ~possLMER(data = .x,
                                   var_value ~ exam_order + disease_subtype + bmi + age_m0 + 
                                     (1 | immet_id))),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE)))

res_subtype_act_02_contrast <- res_subtype_act_02 |> unnest(tidier) |> 
  filter(str_detect(term, "disease_subtype") &
           (p.value < 0.06) & (sign(conf.low) == sign(conf.high)) ) |> 
  mutate(emm_subtype_by_exam = map(mod,
    ~ safe_emm(.x, specs = ~ disease_subtype | exam_order)),
    contr_subtype_by_exam = map(emm_subtype_by_exam,
      ~ safe_contr(.x, method = "pairwise", adjust = "tukey")),
    contr_tidy = map(contr_subtype_by_exam,
      ~ if (is.null(.x)) tibble() else tidy(.x, conf.int = TRUE))) |> 
  select(var_name, contr_tidy) |> 
  unnest(contr_tidy)

res_response_act_02 <- d15_scores |> 
  select(immet_id, exam_order, all_of(x_activ), all_of(base_covars), response_m0_m6) |> 
  pivot_longer(cols = all_of(x_activ), 
               names_to = "var_name", values_to = "var_value") |> 
  group_by(var_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, ~possLMER(data = .x,
                                   var_value ~ exam_order + response_m0_m6 + bmi + age_m0 + 
                                     (1 | immet_id))),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE)))

res_response_act_02_contrast <- res_response_act_02 |> unnest(tidier) |> 
  filter(str_detect(term, "response") &
           (p.value < 0.06) & (sign(conf.low) == sign(conf.high)) ) |> 
  mutate(emm_subtype_by_exam = map(mod,
                                   ~ safe_emm(.x, specs = ~ disease_subtype | exam_order)),
         contr_subtype_by_exam = map(emm_subtype_by_exam,
                                     ~ safe_contr(.x, method = "pairwise", adjust = "tukey")),
         contr_tidy = map(contr_subtype_by_exam,
                          ~ if (is.null(.x)) tibble() else tidy(.x, conf.int = TRUE))) |> 
  select(var_name, contr_tidy) |> 
  unnest(contr_tidy)

# 
res_response_act_02b <- d15_scores |> 
  select(immet_id, exam_order, all_of(x_activ), all_of(base_covars), response_m0_m6) |> 
  mutate(response_m0_m6 = factor(response_m0_m6, ordered = FALSE)) |> 
  pivot_longer(cols = all_of(x_activ), 
               names_to = "var_name", values_to = "var_value") |> 
  group_by(var_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, ~possLMER(data = .x,
                                   var_value ~ exam_order + response_m0_m6 + bmi + age_m0 + 
                                     (1 | immet_id))),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE)))

res_response_act_02b_contrast <- res_response_act_02b |> unnest(tidier) |> 
  filter(str_detect(term, "response") &
           (p.value < 0.06) & (sign(conf.low) == sign(conf.high)) ) |> 
  mutate(emm_subtype_by_exam = map(mod,
                                   ~ safe_emm(.x, specs = ~ disease_subtype | exam_order)),
         contr_subtype_by_exam = map(emm_subtype_by_exam,
                                     ~ safe_contr(.x, method = "pairwise", adjust = "tukey")),
         contr_tidy = map(contr_subtype_by_exam,
                          ~ if (is.null(.x)) tibble() else tidy(.x, conf.int = TRUE))) |> 
  select(var_name, contr_tidy) |> 
  unnest(contr_tidy)

# **********************************************************************
### response ----
# **********************************************************************
res_response_time_01 <- d15_scores |>
  select(immet_id, exam_order, all_of(x_time), all_of(base_covars), response_m0_m6) |>
  mutate(
    exam_order     = factor(exam_order),
    response_m0_m6 = factor(response_m0_m6, ordered = TRUE)  # necháš jako ordered, ale níže testuješ přes emmeans
  ) |>
  pivot_longer(
    cols = all_of(x_time),
    names_to = "var_name",
    values_to = "var_value"
  ) |>
  group_by(var_name) |>
  nest() |>
  ungroup() |>
  mutate(
    mod = map(data, ~ possLMER(
      data = .x,
      var_value ~ exam_order + response_m0_m6 + bmi + age_m0 + (1 | immet_id)
    )),
    tidier = map(mod, ~ tidy(.x, effects = "fixed", conf.int = TRUE))
  )

# 1) Které time parametry mají významný efekt response (v modelu)?
sig_response_time_01 <- res_response_time_01 |>
  transmute(var_name, tidier) |>
  unnest(tidier) |>
  filter(str_detect(term, "^response_m0_m6")) |>
  group_by(var_name) |>
  summarise(
    min_p_response = min(p.value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(min_p_response < 0.05)

sig_response_time_01
# -> Tohle je seznam time proměnných, kde vyšel response efekt (alespoň jeden kontrast)

# 2) Pro významné proměnné: emmeans(response | exam_order) + pairwise kontrasty
res_response_time_01_contrast <- res_response_time_01 |>
  semi_join(sig_response_time_01, by = "var_name") |>
  mutate(
    emm_resp_by_exam = map(mod, ~ safe_emm(.x, specs = ~ response_m0_m6 | exam_order)),
    contr_resp_by_exam = map(emm_resp_by_exam, ~ safe_contr(.x, method = "pairwise", adjust = "tukey")),
    contr_tidy = map(contr_resp_by_exam, ~ if (is.null(.x)) tibble() else tidy(.x, conf.int = TRUE))
  ) |>
  select(var_name, contr_tidy) |>
  unnest(contr_tidy) |>
  # volitelně: jen "jasně významné" kontrasty (p<0.05 a CI mimo 0)
  filter(adj.p.value < 0.06, sign(conf.low) == sign(conf.high)) |>
  arrange(var_name, exam_order, adj.p.value)

res_response_time_01_contrast
# -> konkrétní významné rozdíly mezi response skupinami zvlášť pro M0 a M6


d14_join |> 
  select(immet_id, exam_order, response_m0_m6, disease_subtype, bmi, age_m0) |> 
  mutate(response_m0_m6 = response_m0_m6 |> 
           tolower() |> 
           factor(levels = c("no improvement", "minimal", 
                                "moderate", "major")),
         mod = ordinal::clmm(response_m0_m6 ~ exam_order * disease_subtype + 
                               bmi + age_m0 + (1|immet_id)))
# **********************************************************************
### response - ordered ----
# **********************************************************************

possCLMM <- function(data, formula, ...) {
  tryCatch(
    ordinal::clmm(formula, data = data, Hess = TRUE, ...),
    error = function(e) NULL
  )
}

library(dplyr)
library(tidyr)
library(purrr)
library(ordinal)
library(emmeans)
library(broom)

res_ord_01 <- d15_scores |>
  select(immet_id, exam_order, response_m0_m6, disease_subtype, bmi, age_m0) |>
  mutate(
    response_m0_m6 = response_m0_m6 |>
      tolower() |>
      factor(levels = c("no improvement", "minimal", "moderate", "major")) |>
      ordered(),
    exam_order = factor(exam_order),
    disease_subtype = factor(disease_subtype)
  ) |>
  mutate(.grp = "all") |>
  group_by(.grp) |>
  nest() |>
  ungroup() |>
  mutate(
    mod = map(data, ~ ordinal::clmm(
      data = .x,
      response_m0_m6 ~ exam_order + disease_subtype + bmi + age_m0 + (1 | immet_id),
      link = "logit"
    )),
    mod_m0 = map(data, ~ ordinal::clmm(
      data = .x,
      response_m0_m6 ~ exam_order + bmi + age_m0 + (1 | immet_id),
      link = "logit"
    )),
    emms_lat = map(mod, ~ emmeans(.x, ~ disease_subtype | exam_order, mode = "latent")),
    contrasts_lat = map(emms_lat, ~ pairs(.x, adjust = "tukey")),
    contrasts_lat_tbl = map(contrasts_lat, ~ broom::tidy(.x, conf.int = TRUE)),
    emms_prob = map(mod, ~ emmeans(.x, ~ disease_subtype | exam_order, mode = "prob")),
    emms_prob_tbl = map(emms_prob, as.data.frame)
  ) |>
  select(-.grp)

res_ord_01_tab_summary <- res_ord_01$mod[[1]] |> tidy(conf.int = TRUE)
res_ord_01_tab_contrast <- res_ord_01$contrasts_lat_tbl[[1]]

res_ord_01_anova <- anova(
  res_ord_01$mod_m0[[1]],
  res_ord_01$mod[[1]]
) |> as.data.frame()


export(list(model_sum = res_ord_01_tab_summary,
            model_cotrast = res_ord_01_tab_contrast,
            model_comp_m0 = res_ord_01_anova),
       file.path("output/tables/final", paste0(out_date, "_resp_diseas_type_output_01.xlsx")))

# **********************************************************************
# 4. figures ----
# **********************************************************************
## line graf ----
# **********************************************************************


fig_sel_01 <- c("enmo_full_recording_mean", "inactivity", "light", 
                "moderate", "t400_time_min")

fig_sel_02 <- c("ig_gradient_enmo", "ig_intercept_enmo", "m5_enmo")

fig_sel_03 <- c("t100_b1m80_time_min", "t100_b5m80_time_min", "t100_b10m80_time_min")

fig_sel_04 <- c("p9931_enmo", "p9861_enmo", "p9792_enmo", "p9583_enmo","p9167_enmo")



# ************************************************************
# FIG: line plots M0 vs M6 by response, 4 vars, patchwork 
# ************************************************************

fig_line_m0m6_respcol_build_df <- function(
    data,
    vars_sel,
    id_col = "immet_id",
    response_col = "response_m0_m6",
    exam_col = "exam_order",
    exam_levels = c("m0", "m6"),
    out_var_name = "var_name",
    out_var_value = "var_value"
) {
  data |>
    dplyr::select(all_of(c(vars_sel, id_col, response_col, exam_col))) |>
    dplyr::mutate(
      "{exam_col}" := factor(tolower(as.character(.data[[exam_col]])), levels = exam_levels),
      "{response_col}" := factor(.data[[response_col]])
    ) |>
    tidyr::pivot_longer(
      cols = all_of(vars_sel),
      names_to = out_var_name,
      values_to = out_var_value
    ) |>
    dplyr::mutate("{out_var_name}" := factor(.data[[out_var_name]], levels = vars_sel))
}

fig_line_m0m6_respcol_make_var_column <- function(
    df_var,
    var_label,
    id_col = "immet_id",
    response_col = "response_m0_m6",
    exam_col = "exam_order",
    value_col = "var_value",
    color_values = c(
      "No Improvement" = "blue",
      "Minimal"        = "red",
      "Moderate"       = "darkgreen",
      "Major"          = "#9900FF"
    )
) {
  ggplot(df_var, aes(
    x = .data[[exam_col]],
    y = .data[[value_col]],
    col = .data[[response_col]],
    shape = .data[[response_col]],
    group = .data[[id_col]]
  )) +
    geom_line(linewidth = 0.9, alpha = 0.5) +
    geom_point(
      position = position_jitter(width = 0.08, height = 0),
      size = 2, alpha = 0.5
    ) +
    scale_color_manual(values = color_values) +
    labs(title = var_label) +
    facet_grid(rows = vars(.data[[response_col]]), scales = "fixed") +
    sjPlot::theme_sjplot2() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text  = element_text(size = 10, face = "bold"),
      strip.text.y = element_text(size = 11, face = "bold"),
      legend.position = "none"
    )
}

fig_line_m0m6_respcol_patchwork_4cols <- function(
    data,
    vars_sel,
    id_col = "immet_id",
    response_col = "response_m0_m6",
    exam_col = "exam_order",
    exam_levels = c("m0", "m6"),
    color_values = c(
      "No Improvement" = "blue",
      "Minimal"        = "red",
      "Moderate"       = "darkgreen",
      "Major"          = "#9900FF"
    )
) {
  fig_line_m0m6_respcol_df <- fig_line_m0m6_respcol_build_df(
    data = data,
    vars_sel = vars_sel,
    id_col = id_col,
    response_col = response_col,
    exam_col = exam_col,
    exam_levels = exam_levels
  )
  
  # 1 plot = 1 sloupec (1 var_name), uvnitř 4 řádky (response)
  fig_line_m0m6_respcol_plots <- fig_line_m0m6_respcol_df |>
    group_by(var_name) |>
    group_split() |>
    imap(~ fig_line_m0m6_respcol_make_var_column(
      df_var = .x,
      var_label = as.character(unique(.x$var_name)),
      id_col = id_col,
      response_col = response_col,
      exam_col = exam_col,
      value_col = "var_value",
      color_values = color_values
    ))
  
  # složit do 1 řady (4 sloupce); každý sloupec má vlastní y (protože je to samostatný ggplot)
  fig_line_m0m6_respcol_patch <- patchwork::wrap_plots(fig_line_m0m6_respcol_plots, nrow = 1)
  
  list(
    fig_line_m0m6_respcol_df    = fig_line_m0m6_respcol_df,
    fig_line_m0m6_respcol_plots = fig_line_m0m6_respcol_plots,
    fig_line_m0m6_respcol_patch = fig_line_m0m6_respcol_patch
  )
}

# ---- použít:
fig_line_m0m6_respcol_out <- fig_line_m0m6_respcol_patchwork_4cols(
  data = d14_join,
  vars_sel = fig_sel_01
)

fig_line_01 <- fig_line_m0m6_respcol_out$fig_line_m0m6_respcol_patch
tiff(
  filename = file.path("output/figures/final",
                       paste0(format(Sys.Date(), "%y%m%d"),
                       "_", 
                       "slope_line_chart_01.tiff")),
  compression = "lzw",
  res = "150",
  width = 1400,
  height = 1000
)
print(fig_line_01)
dev.off()

# **********************************************************************
# 5. myokines ----
# **********************************************************************
d17_myokines <- import(paste0(folder_path, "Treatment Response_sérum,expr_pre STAT .xlsx")) |> 
  select(!contains("pg")) |> 
  rename(immet_id = `projekt ID`) |> 
  mutate(immet_id = str_remove_all(immet_id, "_"))


# **********************************************************************
# Data: jen numerické proměnné
# **********************************************************************
# **********************************************************************
# PCA myokines (tidymodels-native)
# **********************************************************************

# 1) Recipe + prep
pca_myokines_rec <- recipe(~ ., data = d17_myokines) |>
  update_role(immet_id, new_role = "id") |>
  step_zv(all_numeric_predictors(), id = "zv") |>
  step_normalize(all_numeric_predictors(), id = "norm") |>
  step_pca(all_numeric_predictors(), threshold = 0.90, id = "pca")

pca_myokines_prep <- prep(pca_myokines_rec)

# 2) Skupiny
pca_myokines_grp <- d14_join |>
  select(immet_id, response_m0_m6) |>
  distinct()

# 3) Scores + skupiny
pca_myokines_scores <- juice(pca_myokines_prep) |>
  left_join(pca_myokines_grp, by = "immet_id")

# 4) Explained variance (%)
pca_myokines_var <- tidy(pca_myokines_prep, id = "pca", type = "variance") |>
  filter(terms == "variance") |>
  mutate(
    percent  = 100 * value / sum(value),
    cum_perc = cumsum(percent)
  )

pc1_perc <- pca_myokines_var |> filter(component == 1) |> pull(percent)
pc2_perc <- pca_myokines_var |> filter(component == 2) |> pull(percent)

# 5) Scree plot
fig_pca_myokines_scree <-
  ggplot(pca_myokines_var, aes(component, percent)) +
  geom_col() +
  geom_line(aes(group = 1)) +
  geom_point() +
  labs(y = "% explained variance", x = "PC") +
  theme_minimal(base_size = 12)

# 6) PC1 vs PC2 plot (filled ellipses)
fig_pca_myokines_pc1_pc2 <-
  ggplot(
    pca_myokines_scores,
    aes(PC1, PC2, color = response_m0_m6, fill = response_m0_m6)
  ) +
  geom_point(size = 2.6, alpha = 0.85) +
  stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.15, linewidth = 1) +
  labs(
    title = "PCA - myokines: PC1 vs PC2",
    x = paste0("PC1 (", round(pc1_perc, 1), "%)"),
    y = paste0("PC2 (", round(pc2_perc, 1), "%)"),
    color = "Response",
    fill  = "Response"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 20),
    axis.title   = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold")
  ) +
  paletteer::scale_color_paletteer_d("MetBrewer::Egypt") +
  paletteer::scale_fill_paletteer_d("MetBrewer::Egypt")

# 7) Export PC1–PC2
tiff(
  filename = file.path(
    "output/figures",
    paste0(format(Sys.Date(), "%y%m%d"), "_pca_myokines_pc1_pc2.tiff")
  ),
  compression = "lzw",
  res = 150,
  width = 1600,
  height = 1200
)
print(fig_pca_myokines_pc1_pc2)
dev.off()

# (volitelně) zobrazit v RStudio
fig_pca_myokines_scree
fig_pca_myokines_pc1_pc2

# **********************************************************************
# 6. lmer - association analyses - PtVAS covariate ----
# **********************************************************************
## clinics ----
# **********************************************************************
### time ----
# **********************************************************************
res_ass_time_clincs_01 <- d14_join |> 
  select(immet_id, exam_order, all_of(sel_time_01), all_of(sel_clinics_01), physician_vas, age_m0) |> 
  pivot_longer(cols = all_of(sel_time_01),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_clinics_01),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order*scale(indep_value) + scale(physician_vas) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_time_clincs_01_tab <- res_ass_time_clincs_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "scale\\(indep_value\\)", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_time_clincs_01_fig <- res_ass_time_clincs_01 |> 
  inner_join(res_ass_time_clincs_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Clinical and Time features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 11, face = "bold")
                      )
  ))


pdf(file.path("output/figures", 
              paste0(format(Sys.Date(), "%y%m%d"),
                     "_", 
                     "assoc_time_clin_01.pdf")), 
    width = 10, height = 7)
walk(res_ass_time_clincs_01_fig$fig, print)
dev.off()

export(res_ass_time_clincs_01_tab, file.path("output/tables", paste0(out_date, "_ass_time_clincs_01.xlsx")))

# **********************************************************************
### activity ----
# **********************************************************************
res_ass_activ_clincs_01 <- d14_join |> 
  select(immet_id, exam_order, all_of(sel_activ_01), all_of(sel_clinics_01), physician_vas, age_m0) |> 
  pivot_longer(cols = all_of(sel_activ_01),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_clinics_01),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order*scale(indep_value) + scale(physician_vas) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_activ_clincs_01_tab <- res_ass_activ_clincs_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "scale\\(indep_value\\)", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_activ_clincs_01_fig <- res_ass_activ_clincs_01 |> 
  inner_join(res_ass_activ_clincs_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Clinical and Activity features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 11, face = "bold")
                      )
  ))


pdf(file.path("output/figures", 
              paste0(format(Sys.Date(), "%y%m%d"),
                     "_", 
                     "assoc_activ_clin_01.pdf")), 
    width = 10, height = 7)
walk(res_ass_activ_clincs_01_fig$fig, print)
dev.off()

export(res_ass_activ_clincs_01_tab, file.path("output/tables", paste0(out_date, "_ass_activ_clincs_01.xlsx")))

# **********************************************************************
##  myokines ----
# **********************************************************************
### time ----
# **********************************************************************
res_ass_time_myokin_01 <- d14_join |> 
  inner_join(d17_myokines, by = "immet_id") |> 
  select(immet_id, exam_order, all_of(sel_time_01), all_of(sel_myokines_01),  physician_vas, age_m0) |> 
  pivot_longer(cols = all_of(sel_time_01),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_myokines_01),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order + indep_value + scale(physician_vas) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_time_myokin_01_tab <- res_ass_time_myokin_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "indep_value", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_time_myokin_01_fig <- res_ass_time_myokin_01 |> 
  inner_join(res_ass_time_myokin_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Myokines and Time features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 11, face = "bold")
                      )
  ))

# # No significant 
# pdf(file.path("output/figures", 
#               paste0(format(Sys.Date(), "%y%m%d"),
#                      "_", 
#                      "assoc_time_myokin_01.pdf")), 
#     width = 10, height = 7)
# walk(res_ass_time_myokin_01_fig$fig, print)
# dev.off()

export(res_ass_time_myokin_01_tab, file.path("output/tables", paste0(out_date, "_ass_time_myokin_01.xlsx")))

# **********************************************************************
### activity ----
# **********************************************************************
res_ass_activ_myokin_01 <- d14_join |> 
  inner_join(d17_myokines, by = "immet_id") |> 
  select(immet_id, exam_order, all_of(sel_activ_01), all_of(sel_myokines_01), physician_vas, age_m0) |> 
  pivot_longer(cols = all_of(sel_activ_01),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_myokines_01),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order + indep_value + scale(physician_vas) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_activ_myokin_01_tab <- res_ass_activ_myokin_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "indep_value", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_activ_myokin_01_fig <- res_ass_activ_myokin_01 |> 
  inner_join(res_ass_activ_myokin_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Myokines and Activity features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 11, face = "bold")
                      )
  ))

# # No significant 
# pdf(file.path("output/figures", 
#               paste0(format(Sys.Date(), "%y%m%d"),
#                      "_", 
#                      "assoc_activ_myokin_01.pdf")), 
#     width = 10, height = 7)
# walk(res_ass_activ_myokin_01_fig$fig, print)
# dev.off()

export(res_ass_activ_myokin_01_tab, file.path("output/tables", paste0(out_date, "_ass_activ_myokin_01.xlsx")))

# **********************************************************************
# 7. lmer - association analyses - bcm covariate ----
# **********************************************************************
# **********************************************************************
## clinics ----
# **********************************************************************
### time ----
# **********************************************************************
res_ass_time_clincs_bcmADJ_01 <- d14_join |> 
  mutate(creatinine = dplyr::if_else(creatinine > 120, NA_real_, creatinine)) |> 
  select(immet_id, exam_order, all_of(sel_time_02), all_of(sel_clinics_02), bcm, age_m0) |> 
  pivot_longer(cols = all_of(sel_time_02),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_clinics_02),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order*scale(indep_value) + scale(bcm) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_time_clincs_bcmADJ_01_tab <- res_ass_time_clincs_bcmADJ_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "scale\\(indep_value\\)", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_time_clincs_bcmADJ_01_fig <- res_ass_time_clincs_bcmADJ_01 |> 
  inner_join(res_ass_time_clincs_bcmADJ_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 show.data = TRUE,
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Clinical and Time features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 18, face = "bold")
                      )
  ))


pdf(file.path("output/figures", 
              paste0(format(Sys.Date(), "%y%m%d"),
                     "_", 
                     "assoc_time_clin_bcmADJ_01.pdf")), 
    width = 10, height = 7)
walk(res_ass_time_clincs_bcmADJ_01_fig$fig, print)
dev.off()

export(res_ass_time_clincs_bcmADJ_01_tab, file.path("output/tables", paste0(out_date, "_ass_time_clincs_bcmADJ_01.xlsx")))

# **********************************************************************
### activity ----
# **********************************************************************
res_ass_activ_clincs_bcmADJ_01 <- d14_join |> 
  mutate(creatinine = dplyr::if_else(creatinine > 120, NA_real_, creatinine)) |> 
  select(immet_id, exam_order, all_of(sel_activ_02), all_of(sel_clinics_02), bcm, age_m0) |> 
  pivot_longer(cols = all_of(sel_activ_02),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_clinics_02),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order*scale(indep_value) + scale(bcm) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_activ_clincs_bcmADJ_01_tab <- res_ass_activ_clincs_bcmADJ_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "scale\\(indep_value\\)", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_activ_clincs_bcmADJ_01_fig <- res_ass_activ_clincs_bcmADJ_01 |> 
  inner_join(res_ass_activ_clincs_bcmADJ_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 show.data = TRUE,
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Clinical and Activity features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 18, face = "bold")
                      )
  ))


pdf(file.path("output/figures", 
              paste0(format(Sys.Date(), "%y%m%d"),
                     "_", 
                     "assoc_activ_clin_bcmADJ_01.pdf")), 
    width = 10, height = 7)
walk(res_ass_activ_clincs_bcmADJ_01_fig$fig, print)
dev.off()

export(res_ass_activ_clincs_bcmADJ_01_tab, file.path("output/tables", paste0(out_date, "_ass_activ_clincs_bcmADJ_01.xlsx")))

# **********************************************************************
## myokines ----
# **********************************************************************
### time ----
# **********************************************************************
res_ass_time_myokin_bcmADJ_01 <- d14_join |> 
  mutate(creatinine = dplyr::if_else(creatinine > 120, NA_real_, creatinine)) |> 
  inner_join(d17_myokines, by = "immet_id") |> 
  select(immet_id, exam_order, all_of(sel_time_02), all_of(sel_myokines_01),  bcm, age_m0) |> 
  pivot_longer(cols = all_of(sel_time_02),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_myokines_01),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order + indep_value + scale(bcm) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_time_myokin_bcmADJ_01_tab <- res_ass_time_myokin_bcmADJ_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "indep_value", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_time_myokin_bcmADJ_01_fig <- res_ass_time_myokin_bcmADJ_01 |> 
  inner_join(res_ass_time_myokin_bcmADJ_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 show.data = TRUE,
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Myokines and Time features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 11, face = "bold")
                      )
  ))


# pdf(file.path("output/figures", 
#               paste0(format(Sys.Date(), "%y%m%d"),
#                      "_", 
#                      "assoc_time_myokin_bcmADJ_01.pdf")), 
#     width = 10, height = 7)
# walk(res_ass_time_myokin_bcmADJ_01_fig$fig, print)
# dev.off()

export(res_ass_time_myokin_bcmADJ_01_tab, file.path("output/tables", paste0(out_date, "_ass_time_myokin_bcmADJ_01.xlsx")))

# **********************************************************************
### activity ----
# **********************************************************************
res_ass_activ_myokin_bcmADJ_01 <- d14_join |> 
  mutate(creatinine = dplyr::if_else(creatinine > 120, NA_real_, creatinine)) |> 
  inner_join(d17_myokines, by = "immet_id") |> 
  select(immet_id, exam_order, all_of(sel_activ_01), all_of(sel_myokines_01), bcm, age_m0) |> 
  pivot_longer(cols = all_of(sel_activ_01),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_myokines_01),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order + indep_value + scale(bcm) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_activ_myokin_bcmADJ_01_tab <- res_ass_activ_myokin_bcmADJ_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "indep_value", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_activ_myokin_bcmADJ_01_fig <- res_ass_activ_myokin_bcmADJ_01 |> 
  inner_join(res_ass_activ_myokin_bcmADJ_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 show.data = TRUE,
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Myokines and Activity features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 11, face = "bold")
                      )
  ))


# pdf(file.path("output/figures", 
#               paste0(format(Sys.Date(), "%y%m%d"),
#                      "_", 
#                      "assoc_activ_myokin_bcmADJ_01.pdf")), 
#     width = 10, height = 7)
# walk(res_ass_activ_myokin_bcmADJ_01_fig$fig, print)
# dev.off()

export(res_ass_activ_myokin_bcmADJ_01_tab, file.path("output/tables", paste0(out_date, "_ass_activ_myokin_bcmADJ_01.xlsx")))

# **********************************************************************
# 8. response to trt - ordinal analysis ----
# **********************************************************************
res_resp_ordinal_01 <- d16_statistika_fi |> 
  filter(exam %in% c("m0", "m6")) |> 
  select(where(~ mean(is.na(.x)) <= 0.20)) |> 
  pivot_longer(cols = -c(immet_id, exam),
               names_to = "var_name",
               values_to = "var_value") |> 
  left_join(d11_clinics |> select(immet_id, disease_subtype, response_m0_m6, bmi), by = "immet_id") |> 
  left_join(age_m0, by = "immet_id") |> 
  group_by(var_name) |> 
  nest() |> 
  mutate(
    mod = map(data, ~ possCLMM(
      data = .x,
      response_m0_m6 ~ exam + var_value + bmi + age_m0 + (1 | immet_id),
      link = "logit"
    )),
    mod_m0 = map(data, ~ possCLMM(
      data = .x,
      response_m0_m6 ~ exam + bmi + age_m0 + (1 | immet_id),
      link = "logit"
    )),
    emms_lat = map(mod, ~ emtrends(.x, ~ var_value | exam, mode = "latent")),
    contrasts_lat = map(emms_lat, ~ pairs(.x, adjust = "tukey")),
    contrasts_lat_tbl = map(contrasts_lat, ~ broom::tidy(.x, conf.int = TRUE)),
    emms_prob = map(mod, ~ emmeans(.x, ~ var_value | exam, mode = "prob")),
    emms_prob_tbl = map(emms_prob, as.data.frame), 
    mod_diff = map2(mod_m0, mod, anova)
  )

})

local({
# -------------------------------------------------------
# d13_join/d13_mofa override: use d10_cas_b
# -------------------------------------------------------
d13_join <- d10_cas_b |>
  transmute(
    immet_id   = as.character(sample),
    exam_order = as.character(exam),
    across(where(is.numeric), identity)
  ) |>
  inner_join(
    d11_clinics_imp |>
      transmute(
        immet_id   = as.character(immet_id),
        exam_order = as.character(exam_order),
        across(-c(immet_id, exam_order), identity)
      ),
    by = c("immet_id", "exam_order")
  ) |>
  distinct() |>
  na.omit() |> 
  group_by(immet_id) |>
  filter(n() == 2) |>
  ungroup()

d13_mofa <- d13_join |>
  mutate(
    sample  = paste0(immet_id, "_", exam_order),
    subtype = as.character(disease_subtype) # convenience col for color_by
  ) |>
  distinct(sample, .keep_all = TRUE)
# -------------------------------------------------------

# **********************************************************************
# 4. rCCA ----   (MYOKINY jako Y blok místo clinics; jen M0 + M6)
# **********************************************************************


# 0) Helpers / config

iim_myo_rcca_myoX_prefix <- paste0(format(Sys.Date(), "%y%m%d"), "_iim_myo_rcca_myoX")
iim_myo_rcca_myoX_outdir <- file.path("output/figures/final")
dir.create(iim_myo_rcca_myoX_outdir, recursive = TRUE, showWarnings = FALSE)

norm_exam_m0m6 <- function(x) {
  x <- tolower(as.character(x))
  dplyr::case_when(
    x %in% c("m0", "0") ~ "m0",
    x %in% c("m6", "6") ~ "m6",
    TRUE ~ x
  )
}

# robustní sjednocení immet_id:
# - pokud je numeric -> vytvoří IMMETxx
# - pokud je character a vypadá jako "IMMET01" -> nechá
# - pokud je character a je jen číslo ("1") -> vytvoří IMMETxx
norm_immet_id <- function(x) {
  if (is.numeric(x)) {
    return(paste0("IMMET", stringr::str_pad(as.integer(x), 2, pad = "0")))
  }
  x_chr <- as.character(x)
  # pokud už začíná IMMET (libovolná velikost písmen), nech
  is_immet <- stringr::str_detect(toupper(x_chr), "^IMMET\\d+")
  out <- x_chr
  # pokud není IMMET a je to čistě číslo -> vytvoř IMMETxx
  is_num_string <- !is_immet & stringr::str_detect(x_chr, "^\\d+$")
  out[is_num_string] <- paste0("IMMET", stringr::str_pad(as.integer(out[is_num_string]), 2, pad = "0"))
  out
}


# 1) MYOKINY (Y block) preparation
#    - bez bbm_id
#    - jen M0/M6
#    - jen numerické myokiny

iim_myo_rcca_myoX_myok_df <- d01_myokiny_2 |>
  dplyr::mutate(
    immet_id   = norm_immet_id(immet_id),
    exam_order = factor(norm_exam_m0m6(exam_order), levels = c("m0", "m6"))
  ) |>
  dplyr::filter(exam_order %in% c("m0", "m6")) |>
  dplyr::select(-bbm_id)

iim_myo_rcca_myoX_myok_cols <- iim_myo_rcca_myoX_myok_df |>
  dplyr::select(-immet_id, -exam_order) |>
  dplyr::select(where(is.numeric)) |>
  names()

# ochrana: kdyby v myokinech nebyl žádný numeric sloupec
stopifnot(length(iim_myo_rcca_myoX_myok_cols) > 0)


# 2) AKTIVITA: X = numeric columns from d10_aktivita
#    Y = myokiny aligned by immet_id + exam_order

iim_myo_rcca_myoX_xcols_activ <- d10_aktivita |>
  dplyr::select(where(is.numeric)) |>
  names()

iim_myo_rcca_myoX_d_activ <- d12_join |>
  dplyr::mutate(
    immet_id   = norm_immet_id(immet_id),
    exam_order = factor(norm_exam_m0m6(exam_order), levels = c("m0", "m6"))
  ) |>
  dplyr::filter(exam_order %in% c("m0", "m6")) |>
  dplyr::inner_join(
    iim_myo_rcca_myoX_myok_df |>
      dplyr::select(immet_id, exam_order, dplyr::all_of(iim_myo_rcca_myoX_myok_cols)),
    by = c("immet_id", "exam_order")
  )

iim_myo_rcca_myoX_X_activ <- iim_myo_rcca_myoX_d_activ |>
  dplyr::select(dplyr::all_of(iim_myo_rcca_myoX_xcols_activ)) |>
  as.matrix()

iim_myo_rcca_myoX_Y_activ <- iim_myo_rcca_myoX_d_activ |>
  dplyr::select(dplyr::all_of(iim_myo_rcca_myoX_myok_cols)) |>
  as.matrix()

stopifnot(nrow(iim_myo_rcca_myoX_X_activ) == nrow(iim_myo_rcca_myoX_Y_activ))

# drop constant cols + scale
iim_myo_rcca_myoX_X_activ <- drop_const_cols(iim_myo_rcca_myoX_X_activ)
iim_myo_rcca_myoX_Y_activ <- drop_const_cols(iim_myo_rcca_myoX_Y_activ)

iim_myo_rcca_myoX_X_activ <- scale(iim_myo_rcca_myoX_X_activ)
iim_myo_rcca_myoX_Y_activ <- scale(iim_myo_rcca_myoX_Y_activ)

# rcc neumí NA -> vyřadit řádky s NA
iim_myo_rcca_myoX_keep_activ <- stats::complete.cases(iim_myo_rcca_myoX_X_activ) &
  stats::complete.cases(iim_myo_rcca_myoX_Y_activ)

iim_myo_rcca_myoX_X_activ <- iim_myo_rcca_myoX_X_activ[iim_myo_rcca_myoX_keep_activ, , drop = FALSE]
iim_myo_rcca_myoX_Y_activ <- iim_myo_rcca_myoX_Y_activ[iim_myo_rcca_myoX_keep_activ, , drop = FALSE]
iim_myo_rcca_myoX_d_activ <- iim_myo_rcca_myoX_d_activ[iim_myo_rcca_myoX_keep_activ, , drop = FALSE]

# po vyřazení řádků ještě jednou konstanty + scale (bezpečné)
iim_myo_rcca_myoX_X_activ <- drop_const_cols(iim_myo_rcca_myoX_X_activ) |> scale()
iim_myo_rcca_myoX_Y_activ <- drop_const_cols(iim_myo_rcca_myoX_Y_activ) |> scale()


# 3) ČAS: X = numeric columns from d10_cas_b
#    Y = myokiny aligned by immet_id + exam_order

iim_myo_rcca_myoX_xcols_time <- d10_cas_b |>
  dplyr::select(where(is.numeric)) |>
  names()

iim_myo_rcca_myoX_d_time <- d13_join |>
  dplyr::mutate(
    immet_id   = norm_immet_id(immet_id),
    exam_order = factor(norm_exam_m0m6(exam_order), levels = c("m0", "m6"))
  ) |>
  dplyr::filter(exam_order %in% c("m0", "m6")) |>
  dplyr::inner_join(
    iim_myo_rcca_myoX_myok_df |>
      dplyr::select(immet_id, exam_order, dplyr::all_of(iim_myo_rcca_myoX_myok_cols)),
    by = c("immet_id", "exam_order")
  )

iim_myo_rcca_myoX_X_time <- iim_myo_rcca_myoX_d_time |>
  dplyr::select(dplyr::all_of(iim_myo_rcca_myoX_xcols_time)) |>
  as.matrix()

iim_myo_rcca_myoX_Y_time <- iim_myo_rcca_myoX_d_time |>
  dplyr::select(dplyr::all_of(iim_myo_rcca_myoX_myok_cols)) |>
  as.matrix()

stopifnot(nrow(iim_myo_rcca_myoX_X_time) == nrow(iim_myo_rcca_myoX_Y_time))

# drop constant cols + scale
iim_myo_rcca_myoX_X_time <- drop_const_cols(iim_myo_rcca_myoX_X_time)
iim_myo_rcca_myoX_Y_time <- drop_const_cols(iim_myo_rcca_myoX_Y_time)

iim_myo_rcca_myoX_X_time <- scale(iim_myo_rcca_myoX_X_time)
iim_myo_rcca_myoX_Y_time <- scale(iim_myo_rcca_myoX_Y_time)

# rcc neumí NA -> vyřadit řádky s NA
iim_myo_rcca_myoX_keep_time <- stats::complete.cases(iim_myo_rcca_myoX_X_time) &
  stats::complete.cases(iim_myo_rcca_myoX_Y_time)

iim_myo_rcca_myoX_X_time <- iim_myo_rcca_myoX_X_time[iim_myo_rcca_myoX_keep_time, , drop = FALSE]
iim_myo_rcca_myoX_Y_time <- iim_myo_rcca_myoX_Y_time[iim_myo_rcca_myoX_keep_time, , drop = FALSE]
iim_myo_rcca_myoX_d_time <- iim_myo_rcca_myoX_d_time[iim_myo_rcca_myoX_keep_time, , drop = FALSE]

# po vyřazení řádků ještě jednou konstanty + scale
iim_myo_rcca_myoX_X_time <- drop_const_cols(iim_myo_rcca_myoX_X_time) |> scale()
iim_myo_rcca_myoX_Y_time <- drop_const_cols(iim_myo_rcca_myoX_Y_time) |> scale()

# **********************************************************************
# 4) Analyses (rCCA fits) ----
# **********************************************************************
iim_myo_rcca_myoX_fit_activ <- mixOmics::rcc(
  X = iim_myo_rcca_myoX_X_activ,
  Y = iim_myo_rcca_myoX_Y_activ,
  ncomp  = 2,
  method = "shrinkage"
)

iim_myo_rcca_myoX_fit_time <- mixOmics::rcc(
  X = iim_myo_rcca_myoX_X_time,
  Y = iim_myo_rcca_myoX_Y_time,
  ncomp  = 2,
  method = "shrinkage"
)

# optional: quick look at top loadings
sort(iim_myo_rcca_myoX_fit_activ$loadings$X[, 1], decreasing = TRUE)[1:10]
sort(iim_myo_rcca_myoX_fit_activ$loadings$Y[, 1], decreasing = TRUE)[1:10]

sort(iim_myo_rcca_myoX_fit_time$loadings$X[, 1], decreasing = TRUE)[1:10]
sort(iim_myo_rcca_myoX_fit_time$loadings$Y[, 1], decreasing = TRUE)[1:10]

# **********************************************************************
# 5) AKTIVITA outputs ---- (Y = MYOKINY)
# **********************************************************************

# --- Canonical scores scatter (comp1)
iim_myo_rcca_myoX_scores_activ <- tibble::tibble(
  immet_id   = iim_myo_rcca_myoX_d_activ$immet_id,
  exam_order = iim_myo_rcca_myoX_d_activ$exam_order,
  X1 = iim_myo_rcca_myoX_fit_activ$variates$X[, 1],
  Y1 = iim_myo_rcca_myoX_fit_activ$variates$Y[, 1],
  X2 = iim_myo_rcca_myoX_fit_activ$variates$X[, 2],
  Y2 = iim_myo_rcca_myoX_fit_activ$variates$Y[, 2]
)

iim_myo_rcca_myoX_fig_can_activ_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_scores_activ,
  ggplot2::aes(X1, Y1, color = exam_order)
) +
  ggplot2::geom_point(alpha = 0.8) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, linetype = 2) +
  ggplot2::labs(
    title = "rCCA: Canonical variates (Component 1)",
    x = "X variate 1 (activity block)",
    y = "Y variate 1 (myokine block)",
    color = "Exam"
  ) +
  sjPlot::theme_sjplot2()

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_canon_comp_activ_01.tiff")),
  compression = "lzw", res = 300, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_can_activ_01)
grDevices::dev.off()

# --- Heatmap: cross-correlation top loadings
iim_myo_rcca_myoX_topX_activ <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_X_activ),
  loading = iim_myo_rcca_myoX_fit_activ$loadings$X[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:15)

iim_myo_rcca_myoX_topY_activ <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_Y_activ),
  loading = iim_myo_rcca_myoX_fit_activ$loadings$Y[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:20)

iim_myo_rcca_myoX_cor_activ <- stats::cor(
  iim_myo_rcca_myoX_X_activ[, iim_myo_rcca_myoX_topX_activ$var, drop = FALSE],
  iim_myo_rcca_myoX_Y_activ[, iim_myo_rcca_myoX_topY_activ$var, drop = FALSE]
)

iim_myo_rcca_myoX_cor_df_activ <- as.data.frame(as.table(iim_myo_rcca_myoX_cor_activ)) |>
  stats::setNames(c("x_var", "y_var", "cor"))

iim_myo_rcca_myoX_fig_hm_activ_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_cor_df_activ,
  ggplot2::aes(x_var, y_var, fill = cor)
) +
  ggplot2::geom_tile() +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Cross-correlation heatmap (top activity vs top myokines)",
    x = "Activity (top loadings)",
    y = "Myokines (top loadings)",
    fill = "cor"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", space = "Lab")

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_heatmap_loading_activ_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_hm_activ_01)
grDevices::dev.off()

# --- Latent axis loadings (barplots)
iim_myo_rcca_myoX_topX_tab_activ <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_X_activ),
  loading = iim_myo_rcca_myoX_fit_activ$loadings$X[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:20)

iim_myo_rcca_myoX_topY_tab_activ <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_Y_activ),
  loading = iim_myo_rcca_myoX_fit_activ$loadings$Y[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:30)

iim_myo_rcca_myoX_load_df_activ <- dplyr::bind_rows(
  iim_myo_rcca_myoX_topX_tab_activ |> dplyr::mutate(type = "Activity"),
  iim_myo_rcca_myoX_topY_tab_activ |> dplyr::mutate(type = "Myokines")
)

iim_myo_rcca_myoX_fig_axis_activ_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_load_df_activ,
  ggplot2::aes(x = loading, y = reorder(var, loading), fill = type)
) +
  ggplot2::geom_col(alpha = 0.8) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::facet_wrap(~type, scales = "free_y") +
  ggplot2::labs(
    x = "Loading (Canonical component 1)",
    y = NULL,
    title = "Opposing poles of the latent activity–myokine axis",
    subtitle = "Variables on opposite sides represent contrasting patterns"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_activity_myokine_axis_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_axis_activ_01)
grDevices::dev.off()

# --- Correlation with canonical axis (tiles)
iim_myo_rcca_myoX_cor_axis_activ <- dplyr::bind_rows(
  tibble::tibble(
    var  = colnames(iim_myo_rcca_myoX_X_activ),
    cor  = stats::cor(iim_myo_rcca_myoX_X_activ, iim_myo_rcca_myoX_fit_activ$variates$X[, 1]),
    type = "Activity"
  ),
  tibble::tibble(
    var  = colnames(iim_myo_rcca_myoX_Y_activ),
    cor  = stats::cor(iim_myo_rcca_myoX_Y_activ, iim_myo_rcca_myoX_fit_activ$variates$Y[, 1]),
    type = "Myokines"
  )
) |>
  dplyr::filter(abs(cor) > 0.15)

iim_myo_rcca_myoX_fig_cor_activ_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_cor_axis_activ,
  ggplot2::aes(x = type, y = reorder(var, cor), fill = cor)
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  ggplot2::labs(
    fill = "Correlation",
    title = "Correlation of individual variables with the latent CCA axis",
    subtitle = "Red vs blue indicate opposing patterns"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_canon_cor_activ_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_cor_activ_01)
grDevices::dev.off()

# --- Subtype and response visualization
iim_myo_rcca_myoX_plot_activ <- iim_myo_rcca_myoX_d_activ |>
  dplyr::transmute(
    immet_id, exam_order,
    subtype  = as.factor(disease_subtype),
    response = as.factor(response_m0_m6),
    compX1   = iim_myo_rcca_myoX_fit_activ$variates$X[, 1],
    compY1   = iim_myo_rcca_myoX_fit_activ$variates$Y[, 1]
  )

iim_myo_rcca_myoX_fig_subtype_activ_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_plot_activ,
  ggplot2::aes(compX1, compY1, color = subtype)
) +
  ggplot2::geom_point(size = 2.8, alpha = 0.85) +
  ggplot2::stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
  ggplot2::labs(
    x = "rCCA component 1 (Activity)",
    y = "rCCA component 1 (Myokines)",
    title = "Patients projected onto the shared activity–myokine space",
    subtitle = "Ellipses highlight subtype clustering (80% t-ellipse)"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_rcca_subtype_activ_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_subtype_activ_01)
grDevices::dev.off()

iim_myo_rcca_myoX_fig_response_activ_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_plot_activ,
  ggplot2::aes(compX1, compY1, color = response)
) +
  ggplot2::geom_point(size = 2.8, alpha = 0.85) +
  ggplot2::stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
  ggplot2::labs(
    x = "rCCA component 1 (Activity)",
    y = "rCCA component 1 (Myokines)",
    title = "Patients projected onto the shared activity–myokine space",
    subtitle = "Ellipses highlight response clustering (80% t-ellipse)"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_rcca_response_activ_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_response_activ_01)
grDevices::dev.off()

# **********************************************************************
# 6) ČAS outputs ---- (Y = MYOKINY)
# **********************************************************************

# --- Canonical scores scatter (comp1)
iim_myo_rcca_myoX_scores_time <- tibble::tibble(
  immet_id   = iim_myo_rcca_myoX_d_time$immet_id,
  exam_order = iim_myo_rcca_myoX_d_time$exam_order,
  X1 = iim_myo_rcca_myoX_fit_time$variates$X[, 1],
  Y1 = iim_myo_rcca_myoX_fit_time$variates$Y[, 1],
  X2 = iim_myo_rcca_myoX_fit_time$variates$X[, 2],
  Y2 = iim_myo_rcca_myoX_fit_time$variates$Y[, 2]
)

iim_myo_rcca_myoX_fig_can_time_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_scores_time,
  ggplot2::aes(X1, Y1, color = exam_order)
) +
  ggplot2::geom_point(alpha = 0.8) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, linetype = 2) +
  ggplot2::labs(
    title = "rCCA: Canonical variates (Component 1)",
    x = "X variate 1 (time block)",
    y = "Y variate 1 (myokine block)",
    color = "Exam"
  ) +
  sjPlot::theme_sjplot2()

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_canon_comp_time_01.tiff")),
  compression = "lzw", res = 300, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_can_time_01)
grDevices::dev.off()

# --- Heatmap: cross-correlation top loadings
iim_myo_rcca_myoX_topX_time <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_X_time),
  loading = iim_myo_rcca_myoX_fit_time$loadings$X[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:15)

iim_myo_rcca_myoX_topY_time <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_Y_time),
  loading = iim_myo_rcca_myoX_fit_time$loadings$Y[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:20)

iim_myo_rcca_myoX_cor_time <- stats::cor(
  iim_myo_rcca_myoX_X_time[, iim_myo_rcca_myoX_topX_time$var, drop = FALSE],
  iim_myo_rcca_myoX_Y_time[, iim_myo_rcca_myoX_topY_time$var, drop = FALSE]
)

iim_myo_rcca_myoX_cor_df_time <- as.data.frame(as.table(iim_myo_rcca_myoX_cor_time)) |>
  stats::setNames(c("x_var", "y_var", "cor"))

iim_myo_rcca_myoX_fig_hm_time_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_cor_df_time,
  ggplot2::aes(x_var, y_var, fill = cor)
) +
  ggplot2::geom_tile() +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Cross-correlation heatmap (top time vs top myokines)",
    x = "Time (top loadings)",
    y = "Myokines (top loadings)",
    fill = "cor"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", space = "Lab")

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_heatmap_loading_time_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_hm_time_01)
grDevices::dev.off()

# --- Latent axis loadings (barplots)
iim_myo_rcca_myoX_topX_tab_time <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_X_time),
  loading = iim_myo_rcca_myoX_fit_time$loadings$X[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:20)

iim_myo_rcca_myoX_topY_tab_time <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_Y_time),
  loading = iim_myo_rcca_myoX_fit_time$loadings$Y[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:30)

iim_myo_rcca_myoX_load_df_time <- dplyr::bind_rows(
  iim_myo_rcca_myoX_topX_tab_time |> dplyr::mutate(type = "Time"),
  iim_myo_rcca_myoX_topY_tab_time |> dplyr::mutate(type = "Myokines")
)

iim_myo_rcca_myoX_fig_axis_time_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_load_df_time,
  ggplot2::aes(x = loading, y = reorder(var, loading), fill = type)
) +
  ggplot2::geom_col(alpha = 0.8) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::facet_wrap(~type, scales = "free_y") +
  ggplot2::labs(
    x = "Loading (Canonical component 1)",
    y = NULL,
    title = "Opposing poles of the latent time–myokine axis",
    subtitle = "Variables on opposite sides represent contrasting patterns"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_time_myokine_axis_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_axis_time_01)
grDevices::dev.off()

# --- Correlation with canonical axis (tiles)
iim_myo_rcca_myoX_cor_axis_time <- dplyr::bind_rows(
  tibble::tibble(
    var  = colnames(iim_myo_rcca_myoX_X_time),
    cor  = stats::cor(iim_myo_rcca_myoX_X_time, iim_myo_rcca_myoX_fit_time$variates$X[, 1]),
    type = "Time"
  ),
  tibble::tibble(
    var  = colnames(iim_myo_rcca_myoX_Y_time),
    cor  = stats::cor(iim_myo_rcca_myoX_Y_time, iim_myo_rcca_myoX_fit_time$variates$Y[, 1]),
    type = "Myokines"
  )
) |>
  dplyr::filter(abs(cor) > 0.15)

iim_myo_rcca_myoX_fig_cor_time_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_cor_axis_time,
  ggplot2::aes(x = type, y = reorder(var, cor), fill = cor)
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  ggplot2::labs(
    fill = "Correlation",
    title = "Correlation of individual variables with the latent CCA axis",
    subtitle = "Red vs blue indicate opposing patterns"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_canon_cor_time_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_cor_time_01)
grDevices::dev.off()

# --- Subtype and response visualization
iim_myo_rcca_myoX_plot_time <- iim_myo_rcca_myoX_d_time |>
  dplyr::transmute(
    immet_id, exam_order,
    subtype  = as.factor(disease_subtype),
    response = as.factor(response_m0_m6),
    compX1   = iim_myo_rcca_myoX_fit_time$variates$X[, 1],
    compY1   = iim_myo_rcca_myoX_fit_time$variates$Y[, 1]
  )

iim_myo_rcca_myoX_fig_response_time_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_plot_time,
  ggplot2::aes(compX1, compY1, color = response)
) +
  ggplot2::geom_point(size = 2.8, alpha = 0.85) +
  ggplot2::stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
  ggplot2::labs(
    x = "rCCA component 1 (Time)",
    y = "rCCA component 1 (Myokines)",
    title = "Patients projected onto the shared time–myokine space",
    subtitle = "Ellipses highlight response clustering (80% t-ellipse)"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_rcca_response_time_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_response_time_01)
grDevices::dev.off()

iim_myo_rcca_myoX_fig_subtype_time_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_plot_time,
  ggplot2::aes(compX1, compY1, color = subtype)
) +
  ggplot2::geom_point(size = 2.8, alpha = 0.85) +
  ggplot2::stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
  ggplot2::labs(
    x = "rCCA component 1 (Time)",
    y = "rCCA component 1 (Myokines)",
    title = "Patients projected onto the shared time–myokine space",
    subtitle = "Ellipses highlight subtype clustering (80% t-ellipse)"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_rcca_subtype_time_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_subtype_time_01)
grDevices::dev.off()


# **********************************************************************
# 10. výstup pro vyhodnocení ----
# **********************************************************************

# ==========================================================
# rCCA -> export tabulek do XLSX (multi-sheet) ----
# ==========================================================

build_rcca_tables <- function(fit, X, Y, meta,
                              analysis_name = "rcca",
                              block_x = "X",
                              block_y = "Y",
                              ncomp = 2,
                              top_n_load_x = 20,
                              top_n_load_y = 30,
                              top_n_pairs_x = 15,
                              top_n_pairs_y = 20,
                              axis_cor_thr = 0.15) {
  
  stopifnot(nrow(X) == nrow(Y))
  stopifnot(nrow(meta) == nrow(X))
  
  # ---- loadings jako matice (robustně)
  loadX <- fit$loadings$X
  loadY <- fit$loadings$Y
  loadX <- if (is.null(dim(loadX))) matrix(loadX, ncol = 1) else as.matrix(loadX)
  loadY <- if (is.null(dim(loadY))) matrix(loadY, ncol = 1) else as.matrix(loadY)
  
  # ---- zajistit jména proměnných (mixOmics někdy nedá rownames)
  if (is.null(rownames(loadX))) rownames(loadX) <- colnames(X)
  if (is.null(rownames(loadY))) rownames(loadY) <- colnames(Y)
  
  # ---- efektivní počet komponent
  ncomp_eff <- min(
    ncomp,
    ncol(loadX), ncol(loadY),
    ncol(fit$variates$X), ncol(fit$variates$Y)
  )
  
  # ---- kanonické korelace
  cc_tbl <- tibble::tibble(
    analysis = analysis_name,
    comp = seq_len(ncomp_eff),
    canonical_cor = purrr::map_dbl(seq_len(ncomp_eff), \(k)
                                   stats::cor(fit$variates$X[, k], fit$variates$Y[, k])
    )
  )
  
  # ---- přehled dat
  summary_tbl <- tibble::tibble(
    analysis = analysis_name,
    n = nrow(X),
    p_x = ncol(X),
    p_y = ncol(Y),
    ncomp = ncomp_eff
  )
  
  # ---- scores + meta
  scores_tbl <- meta |>
    dplyr::mutate(
      X1 = fit$variates$X[, 1],
      Y1 = fit$variates$Y[, 1],
      X2 = if (ncomp_eff >= 2) fit$variates$X[, 2] else NA_real_,
      Y2 = if (ncomp_eff >= 2) fit$variates$Y[, 2] else NA_real_
    )
  
  # ---- loadings long
  load_x_long <- purrr::map_dfr(seq_len(ncomp_eff), \(k) {
    tibble::tibble(
      analysis = analysis_name,
      block = block_x,
      comp = paste0("comp", k),
      var = rownames(loadX),
      loading = loadX[, k],
      abs_loading = abs(loadX[, k])
    )
  })
  
  load_y_long <- purrr::map_dfr(seq_len(ncomp_eff), \(k) {
    tibble::tibble(
      analysis = analysis_name,
      block = block_y,
      comp = paste0("comp", k),
      var = rownames(loadY),
      loading = loadY[, k],
      abs_loading = abs(loadY[, k])
    )
  })
  
  # ---- top loadings (po komponentách)
  top_load_x <- load_x_long |>
    dplyr::group_by(comp) |>
    dplyr::arrange(dplyr::desc(abs_loading), .by_group = TRUE) |>
    dplyr::slice_head(n = top_n_load_x) |>
    dplyr::ungroup()
  
  top_load_y <- load_y_long |>
    dplyr::group_by(comp) |>
    dplyr::arrange(dplyr::desc(abs_loading), .by_group = TRUE) |>
    dplyr::slice_head(n = top_n_load_y) |>
    dplyr::ungroup()
  
  # ---- top proměnné pro cross-correlation (COMP1)
  loadX1 <- loadX[, 1]
  loadY1 <- loadY[, 1]
  
  top_x_vars <- names(sort(abs(loadX1), decreasing = TRUE)) |>
    head(min(top_n_pairs_x, length(loadX1)))
  
  top_y_vars <- names(sort(abs(loadY1), decreasing = TRUE)) |>
    head(min(top_n_pairs_y, length(loadY1)))
  
  if (length(top_x_vars) == 0 || length(top_y_vars) == 0) {
    stop("Top proměnné pro cross-correlation jsou prázdné (loadings mají nulovou délku?).")
  }
  
  cor_raw <- stats::cor(
    X[, top_x_vars, drop = FALSE],
    Y[, top_y_vars, drop = FALSE]
  )
  
  if (is.null(dim(cor_raw))) {
    cor_mat <- matrix(
      cor_raw,
      nrow = length(top_x_vars),
      ncol = length(top_y_vars),
      dimnames = list(top_x_vars, top_y_vars)
    )
  } else {
    cor_mat <- cor_raw
  }
  
  cor_pairs_tbl <- as.data.frame(as.table(cor_mat))
  colnames(cor_pairs_tbl) <- c("x_var", "y_var", "cor")
  
  cor_pairs_tbl <- cor_pairs_tbl |>
    dplyr::mutate(
      analysis = analysis_name,
      x_block = block_x,
      y_block = block_y,
      abs_cor = abs(cor)
    ) |>
    dplyr::arrange(dplyr::desc(abs_cor))
  
  # ---- korelace proměnných s kanonickou osou (comp1)
  axis_cor_tbl <- dplyr::bind_rows(
    tibble::tibble(
      analysis = analysis_name,
      block = block_x,
      var = colnames(X),
      cor_with_axis1 = stats::cor(X, fit$variates$X[, 1])
    ),
    tibble::tibble(
      analysis = analysis_name,
      block = block_y,
      var = colnames(Y),
      cor_with_axis1 = stats::cor(Y, fit$variates$Y[, 1])
    )
  ) |>
    dplyr::mutate(abs_cor = abs(cor_with_axis1)) |>
    dplyr::arrange(dplyr::desc(abs_cor))
  
  axis_cor_tbl_filt <- axis_cor_tbl |>
    dplyr::filter(abs_cor >= axis_cor_thr)
  
  # ---- group summary (pokud existují ve meta)
  group_cols <- intersect(
    c("subtype", "response", "sex", "disease_subtype", "response_m0_m6", "exam_order"),
    colnames(scores_tbl)
  )
  
  group_scores_tbl <- scores_tbl |>
    dplyr::mutate(dplyr::across(dplyr::all_of(group_cols), as.factor)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n = dplyr::n(),
      X1_mean = mean(X1, na.rm = TRUE),
      X1_sd   = sd(X1, na.rm = TRUE),
      Y1_mean = mean(Y1, na.rm = TRUE),
      Y1_sd   = sd(Y1, na.rm = TRUE),
      .groups = "drop"
    )
  
  readme_tbl <- tibble::tibble(
    sheet = c(
      "summary","canonical_correlations","scores",
      "top_loadings_X","top_loadings_Y",
      "loadings_X_all","loadings_Y_all",
      "cor_pairs_top","axis_cor_all","axis_cor_filtered",
      "group_scores"
    ),
    obsah = c(
      "Počet pacientů a proměnných, počet komponent.",
      "Kanonické korelace pro každou komponentu.",
      "Scores (X1/Y1/X2/Y2) + meta.",
      "Top loadings X po komponentách.",
      "Top loadings Y po komponentách.",
      "Všechny loadings X (long).",
      "Všechny loadings Y (long).",
      "Párové korelace top X vs top Y (comp1).",
      "Korelace proměnných s osou comp1.",
      paste0("Filtrované axis korelace |cor| ≥ ", axis_cor_thr, "."),
      "Průměrné scores po skupinách (pokud jsou dostupné)."
    )
  )
  
  list(
    readme = readme_tbl,
    summary = summary_tbl,
    canonical_correlations = cc_tbl,
    scores = scores_tbl,
    top_loadings_X = top_load_x,
    top_loadings_Y = top_load_y,
    loadings_X_all = load_x_long,
    loadings_Y_all = load_y_long,
    cor_pairs_top = cor_pairs_tbl,
    axis_cor_all = axis_cor_tbl,
    axis_cor_filtered = axis_cor_tbl_filt,
    group_scores = group_scores_tbl
  )
}

export_rcca_tables <- function(tables_list, out_xlsx) {
  dir.create(dirname(out_xlsx), showWarnings = FALSE, recursive = TRUE)
  rio::export(tables_list, out_xlsx)
  message("Saved: ", out_xlsx)
}

# ==========================================================
# Export tabulek pro tento projekt (nekolidující názvy) ----
# ==========================================================

# META musí odpovídat řádkům X/Y (po complete.cases filtraci už sedí iim_myo_rcca_myoX_d_*)

iim_myo_rcca_myoX_meta_activ <- iim_myo_rcca_myoX_d_activ |>
  dplyr::select(immet_id, exam_order, dplyr::any_of(c("disease_subtype","response_m0_m6","sex"))) |>
  dplyr::mutate(
    exam_order = factor(tolower(as.character(exam_order)), levels = c("m0","m6"))
  )

iim_myo_rcca_myoX_meta_time <- iim_myo_rcca_myoX_d_time |>
  dplyr::select(immet_id, exam_order, dplyr::any_of(c("disease_subtype","response_m0_m6","sex"))) |>
  dplyr::mutate(
    exam_order = factor(tolower(as.character(exam_order)), levels = c("m0","m6"))
  )

# 1) AKTIVITA × MYOKINY
iim_myo_rcca_myoX_tabs_activ <- build_rcca_tables(
  fit = iim_myo_rcca_myoX_fit_activ,
  X   = iim_myo_rcca_myoX_X_activ,
  Y   = iim_myo_rcca_myoX_Y_activ,
  meta = iim_myo_rcca_myoX_meta_activ,
  analysis_name = "activity_vs_myokines_M0M6",
  block_x = "Activity",
  block_y = "Myokines",
  ncomp = 2,
  top_n_load_x = 25,
  top_n_load_y = 25,
  top_n_pairs_x = 15,
  top_n_pairs_y = 20,
  axis_cor_thr = 0.15
)

export_rcca_tables(
  iim_myo_rcca_myoX_tabs_activ,
  file.path("output/tables", paste0(format(Sys.Date(), "%y%m%d"), "_iim_myo_rcca_activity_vs_myokines_M0M6.xlsx"))
)

# 2) ČAS × MYOKINY
iim_myo_rcca_myoX_tabs_time <- build_rcca_tables(
  fit = iim_myo_rcca_myoX_fit_time,
  X   = iim_myo_rcca_myoX_X_time,
  Y   = iim_myo_rcca_myoX_Y_time,
  meta = iim_myo_rcca_myoX_meta_time,
  analysis_name = "time_vs_myokines_M0M6",
  block_x = "Time",
  block_y = "Myokines",
  ncomp = 2,
  top_n_load_x = 25,
  top_n_load_y = 25,
  top_n_pairs_x = 15,
  top_n_pairs_y = 20,
  axis_cor_thr = 0.15
)

export_rcca_tables(
  iim_myo_rcca_myoX_tabs_time,
  file.path("output/tables", paste0(format(Sys.Date(), "%y%m%d"), "_iim_myo_rcca_time_vs_myokines_M0M6.xlsx"))
)

})



