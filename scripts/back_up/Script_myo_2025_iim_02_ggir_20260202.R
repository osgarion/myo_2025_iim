# 1. Libraries & functions ----
# *******************************************************
# required
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, purrr, conflicted, renv)
# Functions and libraries uploading
list.files("scripts/functions/", pattern="*.*", full.names=TRUE) |>
  map(~source(.))
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
save.image("data/260123_backup.RData")
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
d10_cas |> 
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
d10_cas |> 
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
x_cols_corr_time <- d10_cas |>
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

export(res_adjust_01_tab, "output/tables/260129_adjust_model_performance.xlsx")

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
       "output/tables/260129_stacked_bout_data_01.xlsx")

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

# ČAS: X = numeric columns from d10_cas only
x_cols_cca_time <- d10_cas |>
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
x_cols_time <- d10_cas |>
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
x_cols_diablo_time <- d10_cas |>
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
mofa_x_cols_time <- d10_cas |>
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
# 7. REPORT EXPORT: objects for 260124_ggir_clinics_01 ----
# **********************************************************************
report_rdata <- file.path("reports", "260124_ggir_clinics_01.RData")

objs_to_save <- c(
  # --- data (M0/M6 context) ---
  "d10_aktivita", "d10_cas", "d11_clinics_imp",
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



