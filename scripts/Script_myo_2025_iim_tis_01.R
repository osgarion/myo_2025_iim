
# **********************************************************************
# 1. data ----
# **********************************************************************
## Variables ----
# **********************************************************************
sel_tis_d18_data_01 <- c(
  "extramuscular_global_assessment", "pt_vas", "ph_vas", "mmt8", "haq", "ast", "alt", "ck", "ld",
  "borg", "da_muscle", "crp", "myoglobin",
  "phase_angle", "bcm", "sf36_pcs",  "sf36_mcs", "sf36pf", "sf36rp", "sf36bp",
  "sf36gh", "sf36vt", "sf36sf", "sf36re", "sf36mf"
)

sel_tis_d18_data_02 <- c("extramuscular_global_assessment", "pt_vas", "ph_vas", "mmt8", "haq",
              "ast", "alt", "ck", "ld")

# **********************************************************************
## Import ----
# **********************************************************************
d18_ggir_tis_legend <-import("R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/GENEActiv/GGIR/merged_table_final_klinická.xlsx",
                             which = "legenda") 

d18_ggir_tis_data <- import(
  "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/GENEActiv/GGIR/merged_table_final_klinická.xlsx",
  which = 2
) |>
  rename(any_of(setNames(d18_ggir_tis_legend$names_orig,
                         d18_ggir_tis_legend$names_new))) |>
  mutate(
    exam = factor(
      exam,
      levels = unique(exam)[
        order(as.numeric(sub("^M", "", unique(exam))))
      ]
    )
  ) |>
  droplevels()

# **********************************************************************
## Missing data ----
# **********************************************************************
### Imputation ----
# **********************************************************************
rec_knn_imp <- recipe(~ ., data = d18_ggir_tis_data) |>
  
  # identifikátory
  update_role(file_id, new_role = "id") |>
  update_role(exam, new_role = "group") |>
  
  # škálování numerických proměnných
  # step_normalize(all_numeric_predictors()) |>
  
  # KNN imputace
  step_impute_knn(
    all_of(sel_tis_d18_data_02),
    neighbors = 7
  )

prep_knn <- rec_knn_imp |>
  prep()

d18_tis_knn_imp <- prep_knn |>
  bake(new_data = NULL)

d18_tis_knn_imp

# **********************************************************************
### Figure ----
# **********************************************************************
na_mask <- d18_ggir_tis_data |>
  mutate(
    across(
      all_of(sel_tis_d18_data_02),
      ~ is.na(.x),
      .names = "{.col}_imp"
    )
  )

d18_tis_knn_imp <- d18_tis_knn_imp |>
  bind_cols(
    select(na_mask, ends_with("_imp"))
  )

plot_long <- d18_tis_knn_imp |>
  mutate(
    across(
      all_of(sel_tis_d18_data_02),
      ~ as.numeric(as.character(.x))
    )
  ) |>
  select(file_id, exam, any_of(sel_tis_d18_data_02), ends_with("_imp")) |>
  pivot_longer(
    cols = any_of(sel_tis_d18_data_02),
    names_to = "variable",
    values_to = "value"
  ) |>
  rowwise() |>
  mutate(
    imp_flag = ifelse(
      get(paste0(variable, "_imp")),
      "imputed",
      "observed"
    )
  ) |>
  ungroup()

ggplot(plot_long,
       aes(exam, value, colour = imp_flag)) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  facet_wrap(~variable, scales = "free_y") +
  scale_colour_manual(values = c(
    observed = "black",
    imputed = "red"
  )) +
  ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 30, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 25, face = "italic"),
      axis.title = ggplot2::element_text(size = 22, face = "bold"),
      axis.text = ggplot2::element_text(size = 21),
      strip.text = ggplot2::element_text(size = 21),
  
      panel.background = ggplot2::element_rect(fill = "white"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
  
      panel.grid.major = ggplot2::element_line(color = "grey80", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_line(color = "grey92", linewidth = 0.2),
  
      legend.background = ggplot2::element_rect(fill = "white"),
      legend.box.background = ggplot2::element_rect(fill = "white")
    )

# **********************************************************************
# 2. TIS computation ----
# **********************************************************************
# **********************************************************************
## Settings ----
# **********************************************************************

exam_keep <- c("M0", "M3", "M6", "M18", "M30")
enz_cols  <- c("ast", "alt", "ck", "ld")

# rozsahy škál pro percentuální změnu
measure_range <- c(
  extramuscular_global_assessment = 100,
  pt_vas = 100,
  ph_vas = 100,
  mmt8 = 80,
  haq = 3
)

# DOPLŇTE své ULN hodnoty
# příkladové placeholdery, upravte podle své laboratoře / pohlaví / metodiky
enzyme_uln <- c(
  ast = 0.6,
  alt = 0.6,
  ck  = 2.83,
  ld  = 4.13
)

# adult DM/PM enzyme range multipliers podle vašeho návodu
enzyme_range_mult_adult <- c(
  ast = 3,
  alt = 3,
  ck  = 15,
  ld  = 3
)

# **********************************************************************
## Conversion - percentage -> TIS scale ----
# **********************************************************************
score_tis <- function(x, cuts, scores) {
  out <- cut(
    x,
    breaks = c(-Inf, cuts, Inf),
    labels = scores,
    right = TRUE
  )
  as.numeric(as.character(out))
}

# **********************************************************************
## Input data ----
# **********************************************************************

d_tis_base <- d18_tis_knn_imp |>
  select(immet_id, exam, all_of(sel_tis_d18_data_02)) |>
  filter(exam %in% exam_keep) |>
  arrange(immet_id, exam)

# **********************************************************************
## Baseline information ----
# **********************************************************************
#    - odstraní pacienty bez M0
#    - určí nejabundantnější enzym v M0
#    - uloží baseline hodnoty


m0_info <- d_tis_base |>
  filter(exam == "M0") |>
  rowwise() |>
  mutate(
    enzyme_name = enz_cols[which.max(c_across(all_of(enz_cols)))],
    enzyme_m0 = case_when(
      enzyme_name == "ast" ~ ast,
      enzyme_name == "alt" ~ alt,
      enzyme_name == "ck"  ~ ck,
      enzyme_name == "ld"  ~ ld
    )
  ) |>
  ungroup() |>
  transmute(
    immet_id,
    m0_extramuscular_global_assessment = extramuscular_global_assessment,
    m0_pt_vas = pt_vas,
    m0_ph_vas = ph_vas,
    m0_mmt8 = mmt8,
    m0_haq = haq,
    enzyme_name,
    enzyme_m0
  )

# **********************************************************************
## Main TIS pipeline ----
# **********************************************************************
d18_tis_comp_01 <- d_tis_base |>
  semi_join(m0_info, by = "immet_id") |>
  left_join(m0_info, by = "immet_id") |>
  mutate(
    # koncentrace enzymu vybraného podle M0
    enzyme = case_when(
      enzyme_name == "ast" ~ ast,
      enzyme_name == "alt" ~ alt,
      enzyme_name == "ck"  ~ ck,
      enzyme_name == "ld"  ~ ld
    ),
    
    # změny ve směru "improvement = kladné číslo"
    d_extramuscular_global_assessment =
      m0_extramuscular_global_assessment - extramuscular_global_assessment,
    
    d_pt_vas = m0_pt_vas - pt_vas,
    d_ph_vas = m0_ph_vas - ph_vas,
    
    # u MMT je improvement opačně
    d_mmt8 = mmt8 - m0_mmt8,
    
    d_haq = m0_haq - haq,
    d_enzyme = enzyme_m0 - enzyme
  ) |>
  mutate(
    # percent improvement
    pct_extramuscular_global_assessment =
      100 * d_extramuscular_global_assessment /
      measure_range["extramuscular_global_assessment"],
    
    pct_pt_vas =
      100 * d_pt_vas /
      measure_range["pt_vas"],
    
    pct_ph_vas =
      100 * d_ph_vas /
      measure_range["ph_vas"],
    
    pct_mmt8 =
      100 * d_mmt8 /
      measure_range["mmt8"],
    
    pct_haq =
      100 * d_haq /
      measure_range["haq"],
    
    pct_enzyme =
      100 * d_enzyme /
      (enzyme_uln[enzyme_name] * enzyme_range_mult_adult[enzyme_name])
  ) |>
  mutate(
    # component TIS scores
    
    # Physician Global Activity
    tis_ph_vas = score_tis(
      pct_ph_vas,
      cuts = c(5, 15, 25, 40),
      scores = c(0, 7.5, 15, 17.5, 20)
    ),
    
    # Patient Global Activity
    tis_pt_vas = score_tis(
      pct_pt_vas,
      cuts = c(5, 15, 25, 40),
      scores = c(0, 2.5, 5, 7.5, 10)
    ),
    
    # Manual Muscle Testing
    tis_mmt8 = score_tis(
      pct_mmt8,
      cuts = c(2, 10, 20, 30),
      scores = c(0, 10, 20, 27.5, 32.5)
    ),
    
    # HAQ
    tis_haq = score_tis(
      pct_haq,
      cuts = c(5, 15, 25, 40),
      scores = c(0, 5, 7.5, 7.5, 10)
    ),
    
    # Enzyme
    tis_enzyme = score_tis(
      pct_enzyme,
      cuts = c(5, 15, 25, 40),
      scores = c(0, 2.5, 5, 7.5, 7.5)
    ),
    
    # Extramuscular Activity
    tis_extramuscular = score_tis(
      pct_extramuscular_global_assessment,
      cuts = c(5, 15, 25, 40),
      scores = c(0, 7.5, 12.5, 15, 20)
    )
  ) |>
  mutate(
    tis_total =
      tis_ph_vas +
      tis_pt_vas +
      tis_mmt8 +
      tis_haq +
      tis_enzyme +
      tis_extramuscular,
    
    tis_response = case_when(
      is.na(tis_total) ~ NA_character_,
      tis_total >= 60 ~ "major",
      tis_total >= 40 ~ "moderate",
      tis_total >= 20 ~ "minimal",
      TRUE ~ "none"
    )
  )

# **********************************************************************
## Export ----
# **********************************************************************
export(d18_tis_comp_01, "data/processed/260310_d18_tis_scale_01.xlsx")

# **********************************************************************
## Data compilation ----
# **********************************************************************
d18_ggir_tis_data_02 <- d18_ggir_tis_data |> 
  full_join(d18_tis_comp_01 |> select(immet_id, exam, tis_total, tis_response))

# **********************************************************************
## EDA: TIS total vs activity parameters ----
# **********************************************************************
### Correlation with TIS ----
# **********************************************************************
activity_candidates_01 <- list(
  ig_gradient_enmo = c("ig_gradient_enmo", "ig_gradient"),
  ig_intercept_enmo = c("ig_intercept_enmo", "ig_intercept"),
  m5_enmo = c("m5_enmo"),
  mvpa_t100_time_min = c("mvpa_t100_time_min", "mvpa_t100"),
  p9931_enmo = c("p9931_enmo", "enmo_p99.31"),
  acc_day_spt_wei = c("acc_day_spt_wei", "enmo_day_spt_min")
)

activity_label_map_01 <- c(
  ig_gradient_enmo = "IG gradient ENMO",
  ig_intercept_enmo = "IG intercept ENMO",
  m5_enmo = "M5 ENMO",
  mvpa_t100_time_min = "MVPA t100 time (min)",
  p9931_enmo = "P99.31 ENMO",
  acc_day_spt_wei = "Acc/day SPT weighted"
)

resolve_activity_col_01 <- function(candidates, data_names) {
  match_cols <- intersect(candidates, data_names)
  if (length(match_cols) == 0) {
    return(NA_character_)
  }
  match_cols[[1]]
}

calc_kendall_tis_activity_01 <- function(df) {
  dat_ok <- df |>
    filter(complete.cases(tis_total, activity_value))

  if (nrow(dat_ok) < 3 ||
      n_distinct(dat_ok$tis_total) < 2 ||
      n_distinct(dat_ok$activity_value) < 2) {
    return(tibble(
      n_complete = nrow(dat_ok),
      tau = NA_real_,
      p_value = NA_real_
    ))
  }

  tst <- suppressWarnings(
    cor.test(
      dat_ok$tis_total,
      dat_ok$activity_value,
      method = "kendall",
      exact = FALSE
    )
  )

  tibble(
    n_complete = nrow(dat_ok),
    tau = unname(tst$estimate[[1]]),
    p_value = tst$p.value
  )
}

activity_col_map_01 <- imap_chr(
  activity_candidates_01,
  ~ resolve_activity_col_01(.x, names(d18_ggir_tis_data_02))
)

activity_col_map_01 <- activity_col_map_01[!is.na(activity_col_map_01)]

d18_tis_activity_plot_long_01 <- d18_ggir_tis_data_02 |>
  select(immet_id, exam, tis_total, any_of(unname(activity_col_map_01))) |>
  filter(!is.na(tis_total), !is.na(exam)) |>
  pivot_longer(
    cols = any_of(unname(activity_col_map_01)),
    names_to = "activity_source",
    values_to = "activity_value"
  ) |>
  filter(!is.na(activity_value)) |>
  mutate(
    exam = factor(
      exam,
      levels = exam_keep
    ),
    activity_var = names(activity_col_map_01)[match(activity_source, activity_col_map_01)],
    activity_label = recode(activity_var, !!!activity_label_map_01)
  )

# **********************************************************************
### Correlation over time ----
# **********************************************************************
d18_tis_activity_time_01 <- d18_tis_activity_plot_long_01 |>
  mutate(
    exam_month = readr::parse_number(as.character(exam))
  ) |>
  filter(!is.na(exam_month))

plot_tis_activity_time_01 <- function(activity_i) {
  dat_plot <- d18_tis_activity_time_01 |>
    filter(activity_var == activity_i)

  if (nrow(dat_plot) == 0) {
    return(NULL)
  }

  dat_med <- dat_plot |>
    group_by(exam_month, exam) |>
    summarise(
      tis_median = median(tis_total, na.rm = TRUE),
      activity_median = median(activity_value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(exam_month)

  if (nrow(dat_med) == 0) {
    return(NULL)
  }

  range_tis <- range(dat_med$tis_median, na.rm = TRUE)
  range_activity <- range(dat_med$activity_median, na.rm = TRUE)

  if (!all(is.finite(c(range_tis, range_activity)))) {
    return(NULL)
  }

  if ((range_tis[2] - range_tis[1]) <= 0) {
    range_tis <- range_tis + c(-0.5, 0.5)
  }

  if ((range_activity[2] - range_activity[1]) <= 0) {
    range_activity <- range_activity + c(-0.5, 0.5)
  }

  a <- (range_tis[2] - range_tis[1]) / (range_activity[2] - range_activity[1])
  b <- range_tis[1] - a * range_activity[1]

  dat_med <- dat_med |>
    mutate(activity_scaled = a * activity_median + b)

  ggplot(dat_med, aes(x = exam_month)) +
    geom_point(aes(y = tis_median, color = "TIS total"), size = 3.8) +
    geom_line(aes(y = tis_median, color = "TIS total", group = 1), linewidth = 2.2) +
    geom_point(
      aes(y = activity_scaled, color = activity_label_map_01[[activity_i]]),
      shape = 1,
      size = 3.8,
      stroke = 1.4
    ) +
    geom_line(
      aes(
        y = activity_scaled,
        color = activity_label_map_01[[activity_i]],
        linetype = activity_label_map_01[[activity_i]],
        group = 1
      ),
      linewidth = 2.2
    ) +
    scale_x_continuous(
      breaks = dat_med$exam_month,
      labels = as.character(dat_med$exam)
    ) +
    scale_y_continuous(
      name = "Median TIS total",
      sec.axis = sec_axis(
        ~ (. - b) / a,
        name = paste0("Median ", activity_label_map_01[[activity_i]])
      )
    ) +
    scale_color_manual(
      name = NULL,
      values = c("TIS total" = "black", setNames("darkred", activity_label_map_01[[activity_i]]))
    ) +
    scale_linetype_manual(
      name = NULL,
      values = setNames("dashed", activity_label_map_01[[activity_i]])
    ) +
    labs(
      title = paste0("TIS total and ", activity_label_map_01[[activity_i]], " over time"),
      subtitle = "Exam-wise medians",
      x = "Exam"
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 40, hjust = 0.5),
      plot.subtitle = element_text(face = "bold", size = 34, hjust = 0.5),
      axis.title = element_text(size = 34, face = "bold"),
      axis.text = element_text(size = 32, face = "bold"),
      legend.text = element_text(size = 26, face = "bold")
    )
}

res_tis_activity_time_plot_list_01 <- names(activity_col_map_01) |>
  set_names() |>
  map(plot_tis_activity_time_01) |>
  purrr::compact()

walk2(
  res_tis_activity_time_plot_list_01,
  names(res_tis_activity_time_plot_list_01),
  \(plot_obj, activity_nm) {
    file_nm <- file.path(
      "output",
      "figures",
      "tis",
      paste0(
        format(Sys.Date(), "%y%m%d"),
        "_tis_total_",
        activity_nm,
        "_over_time_01.png"
      )
    )

    if (file.exists(file_nm)) {
      removed_old_01 <- isTRUE(file.remove(file_nm))
      if (!removed_old_01) {
        warning(
          paste0(
            "Existing file could not be removed before overwrite: ",
            file_nm
          )
        )
      }
    }

    ggsave(
      filename = file_nm,
      plot = plot_obj,
      width = 14,
      height = 10,
      units = "in",
      dpi = 300,
      bg = "white",
      limitsize = FALSE
    )
  }
)

plot_tis_activity_multi_01 <- function(activity_i) {
  dat_plot <- d18_tis_activity_plot_long_01 |>
    filter(activity_var == activity_i)

  if (nrow(dat_plot) == 0) {
    return(NULL)
  }

  kendall_stats_01 <- dat_plot |>
    group_by(exam) |>
    group_modify(~ calc_kendall_tis_activity_01(.x)) |>
    ungroup() |>
    mutate(
      kendall_label = case_when(
        is.na(tau) | is.na(p_value) ~ paste0(
          "n = ", n_complete,
          "\nKendall tau = NA",
          "\np = NA"
        ),
        TRUE ~ paste0(
          "n = ", n_complete,
          "\nKendall tau = ", format(round(tau, 3), nsmall = 3),
          "\np = ", format.pval(p_value, digits = 3, eps = 0.001)
        )
      )
    )

  ggplot(dat_plot, aes(x = tis_total, y = activity_value)) +
    geom_point(alpha = 0.7, size = 2.4, colour = "#2b8cbe") +
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 1.1,
      colour = "#d95f0e"
    ) +
    geom_text(
      data = kendall_stats_01,
      aes(x = Inf, y = Inf, label = kendall_label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = 9.2,
      fontface = "bold"
    ) +
    facet_wrap(~exam, ncol = 2, scales = "free_y") +
    labs(
      title = paste0("TIS total vs ", activity_label_map_01[[activity_i]]),
      subtitle = "Facet: exam",
      x = "tis_total",
      y = activity_label_map_01[[activity_i]]
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 40, hjust = 0.5),
      plot.subtitle = element_text(face = "bold", size = 34, hjust = 0.5),
      axis.title = element_text(size = 34, face = "bold"),
      axis.text = element_text(size = 32, face = "bold"),
      strip.text = element_text(face = "bold", size = 34),
      legend.text = element_text(size = 26, face = "bold")
    )
}

res_tis_activity_plot_list_01 <- names(activity_col_map_01) |>
  set_names() |>
  map(plot_tis_activity_multi_01) |>
  purrr::compact()

dir.create(file.path("output", "figures", "tis"), recursive = TRUE, showWarnings = FALSE)

walk2(
  res_tis_activity_plot_list_01,
  names(res_tis_activity_plot_list_01),
  \(plot_obj, activity_nm) {
    file_nm <- file.path(
      "output",
      "figures",
      "tis",
      paste0(
        format(Sys.Date(), "%y%m%d"),
        "_tis_total_",
        activity_nm,
        "_by_exam_01.png"
      )
    )

    if (file.exists(file_nm)) {
      removed_old_01 <- isTRUE(file.remove(file_nm))
      if (!removed_old_01) {
        warning(
          paste0(
            "Existing file could not be removed before overwrite: ",
            file_nm
          )
        )
      }
    }

    ggsave(
      filename = file_nm,
      plot = plot_obj,
      width = 14,
      height = 10,
      units = "in",
      dpi = 300,
      bg = "white",
      limitsize = FALSE
    )
  }
)

# **********************************************************************
### Correlation over time 2 ----
# **********************************************************************
plot_tis_activity_time_dual_02 <- function(activity_i) {
  dat_i <- d18_tis_activity_time_01 |>
    filter(activity_var == activity_i)

  if (nrow(dat_i) == 0) {
    return(NULL)
  }

  dat_med <- dat_i |>
    group_by(exam_month, exam) |>
    summarise(
      tis_median = median(tis_total, na.rm = TRUE),
      activity_median = median(activity_value, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) |>
    arrange(exam_month)

  if (nrow(dat_med) == 0) {
    return(NULL)
  }

  range_tis <- range(dat_med$tis_median, na.rm = TRUE)
  range_activity <- range(dat_med$activity_median, na.rm = TRUE)

  if (!all(is.finite(c(range_tis, range_activity)))) {
    return(NULL)
  }

  if ((range_tis[2] - range_tis[1]) <= 0) {
    range_tis <- range_tis + c(-0.5, 0.5)
  }

  if ((range_activity[2] - range_activity[1]) <= 0) {
    range_activity <- range_activity + c(-0.5, 0.5)
  }

  a <- (range_tis[2] - range_tis[1]) / (range_activity[2] - range_activity[1])
  b <- range_tis[1] - a * range_activity[1]

  dat_med <- dat_med |>
    mutate(activity_scaled = a * activity_median + b)

  p_main <- ggplot(dat_med, aes(x = exam_month)) +
    geom_point(aes(y = tis_median, color = "TIS total"), size = 3.8) +
    geom_line(aes(y = tis_median, color = "TIS total"), linewidth = 2.2) +
    geom_point(
      aes(y = activity_scaled, color = activity_label_map_01[[activity_i]]),
      shape = 1,
      size = 3.8,
      stroke = 1.4
    ) +
    geom_line(
      aes(
        y = activity_scaled,
        color = activity_label_map_01[[activity_i]],
        linetype = activity_label_map_01[[activity_i]]
      ),
      linewidth = 2.2
    ) +
    scale_x_continuous(
      breaks = dat_med$exam_month,
      labels = paste0("M", dat_med$exam_month)
    ) +
    scale_y_continuous(
      name = "Median TIS total",
      sec.axis = sec_axis(
        ~ (. - b) / a,
        name = paste0("Median ", activity_label_map_01[[activity_i]])
      )
    ) +
    scale_color_manual(
      name = "Measure",
      values = c("TIS total" = "black", setNames("darkred", activity_label_map_01[[activity_i]]))
    ) +
    scale_linetype_manual(
      name = "Measure",
      values = setNames("dashed", activity_label_map_01[[activity_i]])
    ) +
    labs(
      title = paste0("TIS total + ", activity_label_map_01[[activity_i]], " over time"),
      subtitle = "Exam-wise medians without error bars",
      x = "Exam"
    ) +
    theme_bw(base_size = 24) +
    theme(
      plot.title = element_text(face = "bold", size = 30),
      plot.subtitle = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 22),
      axis.text = element_text(face = "bold", size = 18),
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(face = "bold", size = 16)
    )

  tab_counts <- dat_i |>
    distinct(immet_id, exam_month, exam) |>
    count(exam_month, exam, name = "n_patients") |>
    arrange(exam_month)

  exam_levels_tab <- tab_counts |>
    pull(exam) |>
    as.character()

  tab_plot_df <- bind_rows(
    tab_counts |>
      transmute(
        exam = factor(as.character(exam), levels = exam_levels_tab),
        row = "Exam",
        value = as.character(exam)
      ),
    tab_counts |>
      transmute(
        exam = factor(as.character(exam), levels = exam_levels_tab),
        row = "N patients",
        value = as.character(n_patients)
      )
  ) |>
    mutate(row = factor(row, levels = c("N patients", "Exam")))

  p_table <- ggplot(tab_plot_df, aes(x = exam, y = row)) +
    geom_tile(fill = "white", color = "grey50", linewidth = 0.6) +
    geom_text(aes(label = value), size = 8.4, fontface = "bold") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 20) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(face = "bold", size = 20),
      axis.ticks = element_blank(),
      plot.margin = margin(t = 0, r = 20, b = 10, l = 20)
    )

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_main / p_table + patchwork::plot_layout(heights = c(8.5, 4))
  } else {
    p_main
  }
}

res_tis_activity_time_dual_plot_list_02 <- names(activity_col_map_01) |>
  set_names() |>
  map(plot_tis_activity_time_dual_02) |>
  purrr::compact()

walk2(
  res_tis_activity_time_dual_plot_list_02,
  names(res_tis_activity_time_dual_plot_list_02),
  \(plot_obj, activity_nm) {
    file_nm <- file.path(
      "output",
      "figures",
      "tis",
      paste0(
        format(Sys.Date(), "%y%m%d"),
        "_tis_total_",
        activity_nm,
        "_over_time_01.png"
      )
    )

    if (file.exists(file_nm)) {
      removed_old_01 <- isTRUE(file.remove(file_nm))
      if (!removed_old_01) {
        warning(
          paste0(
            "Existing file could not be removed before overwrite: ",
            file_nm
          )
        )
      }
    }

    ggsave(
      filename = file_nm,
      plot = plot_obj,
      width = 14,
      height = 12,
      units = "in",
      dpi = 300,
      bg = "white",
      limitsize = FALSE
    )
  }
)

# **********************************************************************
### Comparison TIS vs activity ----
# **********************************************************************
#### BCM ----
# **********************************************************************
d18_tis_activity_bcm_01 <- d18_tis_activity_plot_long_01 |>
  left_join(
    d18_ggir_tis_data_02 |>
      select(immet_id, exam, bcm) |>
      distinct(),
    by = c("immet_id", "exam")
  ) |>
  mutate(bcm = as.numeric(bcm)) |>
  filter(!is.na(bcm))

plot_tis_activity_bcm_dual_01 <- function(activity_i) {
  dat_i <- d18_tis_activity_bcm_01 |>
    filter(activity_var == activity_i) |>
    arrange(exam, bcm)

  if (nrow(dat_i) == 0) {
    return(NULL)
  }

  range_tis <- range(dat_i$tis_total, na.rm = TRUE)
  range_activity <- range(dat_i$activity_value, na.rm = TRUE)

  if (!all(is.finite(c(range_tis, range_activity)))) {
    return(NULL)
  }

  if ((range_tis[2] - range_tis[1]) <= 0) {
    range_tis <- range_tis + c(-0.5, 0.5)
  }

  if ((range_activity[2] - range_activity[1]) <= 0) {
    range_activity <- range_activity + c(-0.5, 0.5)
  }

  a <- (range_tis[2] - range_tis[1]) / (range_activity[2] - range_activity[1])
  b <- range_tis[1] - a * range_activity[1]

  dat_i <- dat_i |>
    mutate(activity_scaled = a * activity_value + b)

  kendall_stats_01 <- dat_i |>
    group_by(exam) |>
    group_modify(
      ~ {
        dat_ok_tis <- .x |>
          filter(complete.cases(bcm, tis_total))

        dat_ok_activity <- .x |>
          filter(complete.cases(bcm, activity_value))

        tau_tis <- NA_real_
        p_tis <- NA_real_
        if (nrow(dat_ok_tis) >= 3 &&
            n_distinct(dat_ok_tis$bcm) >= 2 &&
            n_distinct(dat_ok_tis$tis_total) >= 2) {
          tst_tis <- suppressWarnings(
            cor.test(
              dat_ok_tis$bcm,
              dat_ok_tis$tis_total,
              method = "kendall",
              exact = FALSE
            )
          )
          tau_tis <- unname(tst_tis$estimate[[1]])
          p_tis <- tst_tis$p.value
        }

        tau_activity <- NA_real_
        p_activity <- NA_real_
        if (nrow(dat_ok_activity) >= 3 &&
            n_distinct(dat_ok_activity$bcm) >= 2 &&
            n_distinct(dat_ok_activity$activity_value) >= 2) {
          tst_activity <- suppressWarnings(
            cor.test(
              dat_ok_activity$bcm,
              dat_ok_activity$activity_value,
              method = "kendall",
              exact = FALSE
            )
          )
          tau_activity <- unname(tst_activity$estimate[[1]])
          p_activity <- tst_activity$p.value
        }

        tibble(
          tau_tis = tau_tis,
          p_tis = p_tis,
          tau_activity = tau_activity,
          p_activity = p_activity
        )
      }
    ) |>
    ungroup() |>
    mutate(
      kendall_label = paste0(
        "TIS: tau = ",
        if_else(is.na(tau_tis), "NA", format(round(tau_tis, 3), nsmall = 3)),
        ", p = ",
        if_else(is.na(p_tis), "NA", format.pval(p_tis, digits = 3, eps = 0.001)),
        "\n",
        activity_label_map_01[[activity_i]],
        ": tau = ",
        if_else(is.na(tau_activity), "NA", format(round(tau_activity, 3), nsmall = 3)),
        ", p = ",
        if_else(is.na(p_activity), "NA", format.pval(p_activity, digits = 3, eps = 0.001))
      )
    )

  ggplot(dat_i, aes(x = bcm)) +
    geom_point(aes(y = tis_total, color = "TIS total"), size = 2.8, alpha = 0.85) +
    geom_smooth(
      aes(y = tis_total, color = "TIS total"),
      method = "lm",
      se = FALSE,
      linewidth = 1.6
    ) +
    geom_point(
      aes(y = activity_scaled, color = activity_label_map_01[[activity_i]]),
      shape = 1,
      size = 2.8,
      stroke = 1.1,
      alpha = 0.85
    ) +
    geom_smooth(
      aes(
        y = activity_scaled,
        color = activity_label_map_01[[activity_i]]
      ),
      method = "lm",
      se = FALSE,
      linewidth = 1.6
    ) +
    geom_text(
      data = kendall_stats_01,
      aes(x = Inf, y = Inf, label = kendall_label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = 9.2,
      fontface = "bold"
    ) +
    facet_wrap(~exam, ncol = 2, scales = "free_x") +
    scale_y_continuous(
      name = "TIS total",
      sec.axis = sec_axis(
        ~ (. - b) / a,
        name = activity_label_map_01[[activity_i]]
      )
    ) +
    scale_color_manual(
      name = "Measure",
      values = c("TIS total" = "black", setNames("darkred", activity_label_map_01[[activity_i]]))
    ) +
    labs(
      title = paste0("TIS total + ", activity_label_map_01[[activity_i]], " by BCM"),
      subtitle = "Faceted by exam",
      x = "BCM"
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 40, hjust = 0.5),
      plot.subtitle = element_text(face = "bold", size = 34, hjust = 0.5),
      axis.title = element_text(size = 34, face = "bold"),
      axis.text = element_text(size = 32, face = "bold"),
      strip.text = element_text(face = "bold", size = 34),
      legend.title = element_text(face = "bold", size = 26),
      legend.text = element_text(size = 26, face = "bold")
    )
}

res_tis_activity_bcm_plot_list_01 <- names(activity_col_map_01) |>
  set_names() |>
  map(plot_tis_activity_bcm_dual_01) |>
  purrr::compact()

walk2(
  res_tis_activity_bcm_plot_list_01,
  names(res_tis_activity_bcm_plot_list_01),
  \(plot_obj, activity_nm) {
    file_nm <- file.path(
      "output",
      "figures",
      "tis",
      paste0(
        format(Sys.Date(), "%y%m%d"),
        "_tis_total_",
        activity_nm,
        "_by_bcm_01.png"
      )
    )

    if (file.exists(file_nm)) {
      removed_old_01 <- isTRUE(file.remove(file_nm))
      if (!removed_old_01) {
        warning(
          paste0(
            "Existing file could not be removed before overwrite: ",
            file_nm
          )
        )
      }
    }

    ggsave(
      filename = file_nm,
      plot = plot_obj,
      width = 16,
      height = 12,
      units = "in",
      dpi = 300,
      bg = "white",
      limitsize = FALSE
    )
  }
)

# **********************************************************************
#### Ph-VAS ----
# **********************************************************************
d18_tis_activity_phvas_01 <- d18_tis_activity_plot_long_01 |>
  left_join(
    d18_ggir_tis_data_02 |>
      select(immet_id, exam, ph_vas) |>
      distinct(),
    by = c("immet_id", "exam")
  ) |>
  mutate(ph_vas = as.numeric(ph_vas)) |>
  filter(!is.na(ph_vas))

plot_tis_activity_phvas_dual_01 <- function(activity_i) {
  dat_i <- d18_tis_activity_phvas_01 |>
    filter(activity_var == activity_i) |>
    arrange(exam, ph_vas)

  if (nrow(dat_i) == 0) {
    return(NULL)
  }

  range_tis <- range(dat_i$tis_total, na.rm = TRUE)
  range_activity <- range(dat_i$activity_value, na.rm = TRUE)

  if (!all(is.finite(c(range_tis, range_activity)))) {
    return(NULL)
  }

  if ((range_tis[2] - range_tis[1]) <= 0) {
    range_tis <- range_tis + c(-0.5, 0.5)
  }

  if ((range_activity[2] - range_activity[1]) <= 0) {
    range_activity <- range_activity + c(-0.5, 0.5)
  }

  a <- (range_tis[2] - range_tis[1]) / (range_activity[2] - range_activity[1])
  b <- range_tis[1] - a * range_activity[1]

  dat_i <- dat_i |>
    mutate(activity_scaled = a * activity_value + b)

  kendall_stats_01 <- dat_i |>
    group_by(exam) |>
    group_modify(
      ~ {
        dat_ok_tis <- .x |>
          filter(complete.cases(ph_vas, tis_total))

        dat_ok_activity <- .x |>
          filter(complete.cases(ph_vas, activity_value))

        tau_tis <- NA_real_
        p_tis <- NA_real_
        if (nrow(dat_ok_tis) >= 3 &&
            n_distinct(dat_ok_tis$ph_vas) >= 2 &&
            n_distinct(dat_ok_tis$tis_total) >= 2) {
          tst_tis <- suppressWarnings(
            cor.test(
              dat_ok_tis$ph_vas,
              dat_ok_tis$tis_total,
              method = "kendall",
              exact = FALSE
            )
          )
          tau_tis <- unname(tst_tis$estimate[[1]])
          p_tis <- tst_tis$p.value
        }

        tau_activity <- NA_real_
        p_activity <- NA_real_
        if (nrow(dat_ok_activity) >= 3 &&
            n_distinct(dat_ok_activity$ph_vas) >= 2 &&
            n_distinct(dat_ok_activity$activity_value) >= 2) {
          tst_activity <- suppressWarnings(
            cor.test(
              dat_ok_activity$ph_vas,
              dat_ok_activity$activity_value,
              method = "kendall",
              exact = FALSE
            )
          )
          tau_activity <- unname(tst_activity$estimate[[1]])
          p_activity <- tst_activity$p.value
        }

        tibble(
          tau_tis = tau_tis,
          p_tis = p_tis,
          tau_activity = tau_activity,
          p_activity = p_activity
        )
      }
    ) |>
    ungroup() |>
    mutate(
      kendall_label = paste0(
        "TIS: tau = ",
        if_else(is.na(tau_tis), "NA", format(round(tau_tis, 3), nsmall = 3)),
        ", p = ",
        if_else(is.na(p_tis), "NA", format.pval(p_tis, digits = 3, eps = 0.001)),
        "\n",
        activity_label_map_01[[activity_i]],
        ": tau = ",
        if_else(is.na(tau_activity), "NA", format(round(tau_activity, 3), nsmall = 3)),
        ", p = ",
        if_else(is.na(p_activity), "NA", format.pval(p_activity, digits = 3, eps = 0.001))
      )
    )

  ggplot(dat_i, aes(x = ph_vas)) +
    geom_point(aes(y = tis_total, color = "TIS total"), size = 2.8, alpha = 0.85) +
    geom_smooth(
      aes(y = tis_total, color = "TIS total"),
      method = "lm",
      se = FALSE,
      linewidth = 1.6
    ) +
    geom_point(
      aes(y = activity_scaled, color = activity_label_map_01[[activity_i]]),
      shape = 1,
      size = 2.8,
      stroke = 1.1,
      alpha = 0.85
    ) +
    geom_smooth(
      aes(
        y = activity_scaled,
        color = activity_label_map_01[[activity_i]]
      ),
      method = "lm",
      se = FALSE,
      linewidth = 1.6
    ) +
    geom_text(
      data = kendall_stats_01,
      aes(x = Inf, y = Inf, label = kendall_label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = 9.2,
      fontface = "bold"
    ) +
    facet_wrap(~exam, ncol = 2, scales = "free_x") +
    scale_y_continuous(
      name = "TIS total",
      sec.axis = sec_axis(
        ~ (. - b) / a,
        name = activity_label_map_01[[activity_i]]
      )
    ) +
    scale_color_manual(
      name = "Measure",
      values = c("TIS total" = "black", setNames("darkred", activity_label_map_01[[activity_i]]))
    ) +
    labs(
      title = paste0("TIS total + ", activity_label_map_01[[activity_i]], " by Ph-VAS"),
      subtitle = "Faceted by exam",
      x = "Ph-VAS"
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 40, hjust = 0.5),
      plot.subtitle = element_text(face = "bold", size = 34, hjust = 0.5),
      axis.title = element_text(size = 34, face = "bold"),
      axis.text = element_text(size = 32, face = "bold"),
      strip.text = element_text(face = "bold", size = 34),
      legend.title = element_text(face = "bold", size = 26),
      legend.text = element_text(size = 26, face = "bold")
    )
}

res_tis_activity_phvas_plot_list_01 <- names(activity_col_map_01) |>
  set_names() |>
  map(plot_tis_activity_phvas_dual_01) |>
  purrr::compact()

walk2(
  res_tis_activity_phvas_plot_list_01,
  names(res_tis_activity_phvas_plot_list_01),
  \(plot_obj, activity_nm) {
    file_nm <- file.path(
      "output",
      "figures",
      "tis",
      paste0(
        format(Sys.Date(), "%y%m%d"),
        "_tis_total_",
        activity_nm,
        "_by_phvas_01.png"
      )
    )

    if (file.exists(file_nm)) {
      removed_old_01 <- isTRUE(file.remove(file_nm))
      if (!removed_old_01) {
        warning(
          paste0(
            "Existing file could not be removed before overwrite: ",
            file_nm
          )
        )
      }
    }

    ggsave(
      filename = file_nm,
      plot = plot_obj,
      width = 16,
      height = 12,
      units = "in",
      dpi = 300,
      bg = "white",
      limitsize = FALSE
    )
  }
)

# **********************************************************************
### Model Selection: BCM ----
# **********************************************************************
bcm_predictor_candidates_01 <- c(
  tis_total = "tis_total",
  unname(activity_col_map_01)
)
names(bcm_predictor_candidates_01) <- c("tis_total", names(activity_col_map_01))

predictor_model_label_map_01 <- c(
  tis_total = "TIS total",
  activity_label_map_01
)

age_m0_tis_01 <- d18_ggir_tis_data_02 |>
  filter(exam == "M0") |>
  select(immet_id, age) |>
  distinct() |>
  rename(age_m0 = age)

d18_tis_bcm_model_base_01 <- d18_ggir_tis_data_02 |>
  select(immet_id, exam, bcm, sex, age, tis_total, any_of(unname(activity_col_map_01))) |>
  distinct() |>
  left_join(age_m0_tis_01, by = "immet_id") |>
  mutate(
    bcm = as.numeric(bcm),
    age_m0 = as.numeric(age_m0),
    exam = factor(exam, levels = exam_keep),
    sex = factor(sex),
    immet_id = factor(immet_id)
  )

fit_bcm_models_one_predictor_01 <- function(dat) {
  out <- list(
    model_lmer_linear = NULL,
    model_glmer_gamma = NULL,
    model_gam_smooth = NULL
  )

  dat <- dat |>
    filter(
      !is.na(bcm),
      !is.na(predictor_value),
      !is.na(age_m0),
      !is.na(sex),
      !is.na(exam),
      !is.na(immet_id)
    ) |>
    droplevels()

  if (nrow(dat) < 20 || dplyr::n_distinct(dat$immet_id) < 5 || dplyr::n_distinct(dat$predictor_value) < 3) {
    return(out)
  }

  out$model_lmer_linear <- tryCatch(
    lme4::lmer(
      bcm ~ scale(predictor_value) + exam + scale(age_m0) + sex + (1 | immet_id),
      data = dat,
      REML = FALSE
    ),
    error = function(e) NULL
  )

  if (all(dat$bcm > 0, na.rm = TRUE)) {
    out$model_glmer_gamma <- tryCatch(
      lme4::glmer(
        bcm ~ scale(predictor_value) + exam + scale(age_m0) + sex + (1 | immet_id),
        data = dat,
        family = Gamma(link = "log"),
        control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
      ),
      error = function(e) NULL
    )
  }

  if (requireNamespace("mgcv", quietly = TRUE)) {
    k_use <- max(3, min(5, dplyr::n_distinct(dat$predictor_value) - 1))
    out$model_gam_smooth <- tryCatch(
      mgcv::gam(
        bcm ~ s(predictor_value, k = k_use) + exam + s(age_m0, k = 4) + sex + s(immet_id, bs = "re"),
        data = dat,
        method = "REML"
      ),
      error = function(e) NULL
    )
  }

  out
}

summarise_model_compare_01 <- function(models_named, predictor_nm) {
  if (length(models_named) == 0) {
    return(tibble(
      predictor = predictor_nm,
      model = NA_character_,
      AIC = NA_real_,
      BIC = NA_real_,
      Performance_Score = NA_real_,
      best_model = NA_character_,
      note = "No valid models"
    ))
  }

  perf_tbl <- NULL
  if (length(models_named) >= 2 && requireNamespace("performance", quietly = TRUE)) {
    perf_tbl <- tryCatch(
      do.call(performance::compare_performance, models_named) |>
        as_tibble(),
      error = function(e) NULL
    )
  }

  if (is.null(perf_tbl)) {
    perf_tbl <- imap_dfr(
      models_named,
      \(mod, nm) tibble(
        model = nm,
        AIC = tryCatch(AIC(mod), error = function(e) NA_real_),
        BIC = tryCatch(BIC(mod), error = function(e) NA_real_)
      )
    )
  }

  name_col <- intersect(c("Name", "Model", "model", "name"), names(perf_tbl))
  if (length(name_col) > 0) {
    perf_tbl <- perf_tbl |>
      rename(model = all_of(name_col[[1]]))
  }

  if (!("AIC" %in% names(perf_tbl))) {
    perf_tbl <- perf_tbl |>
      mutate(AIC = map_dbl(model, ~ tryCatch(AIC(models_named[[.x]]), error = function(e) NA_real_)))
  }

  if (!("BIC" %in% names(perf_tbl))) {
    perf_tbl <- perf_tbl |>
      mutate(BIC = map_dbl(model, ~ tryCatch(BIC(models_named[[.x]]), error = function(e) NA_real_)))
  }

  perf_tbl <- perf_tbl |>
    mutate(predictor = predictor_nm, .before = 1)

  best_model_name <- NA_character_
  if ("Performance_Score" %in% names(perf_tbl) && any(is.finite(perf_tbl$Performance_Score))) {
    best_model_name <- perf_tbl |>
      filter(is.finite(Performance_Score)) |>
      arrange(desc(Performance_Score), AIC) |>
      slice_head(n = 1) |>
      pull(model)
  } else if (any(is.finite(perf_tbl$AIC))) {
    best_model_name <- perf_tbl |>
      filter(is.finite(AIC)) |>
      arrange(AIC, BIC) |>
      slice_head(n = 1) |>
      pull(model)
  }

  perf_tbl |>
    mutate(
      best_model = best_model_name,
      note = if_else(model == best_model_name, "Selected", NA_character_)
    )
}

res_bcm_model_selection_01 <- imap(
  bcm_predictor_candidates_01,
  \(predictor_col, predictor_nm) {
    dat_i <- d18_tis_bcm_model_base_01 |>
      transmute(
        immet_id,
        exam,
        bcm,
        sex,
        age_m0,
        predictor_value = as.numeric(.data[[predictor_col]])
      )

    fits <- fit_bcm_models_one_predictor_01(dat_i)

    models_named <- list(
      "LMER linear" = fits$model_lmer_linear,
      "GLMER Gamma" = fits$model_glmer_gamma,
      "GAM smooth + RE" = fits$model_gam_smooth
    )
    models_named <- models_named[!map_lgl(models_named, is.null)]

    compare_tbl <- summarise_model_compare_01(
      models_named = models_named,
      predictor_nm = predictor_model_label_map_01[[predictor_nm]]
    )

    status_tbl <- tibble(
      predictor = predictor_model_label_map_01[[predictor_nm]],
      n = nrow(dat_i |> filter(!is.na(bcm), !is.na(predictor_value), !is.na(age_m0), !is.na(sex), !is.na(exam))),
      n_patients = dplyr::n_distinct(dat_i$immet_id[!is.na(dat_i$predictor_value) & !is.na(dat_i$bcm)]),
      has_lmer = !is.null(fits$model_lmer_linear),
      has_glmer = !is.null(fits$model_glmer_gamma),
      has_gam = !is.null(fits$model_gam_smooth)
    )

    list(compare = compare_tbl, status = status_tbl, fits = fits, data = dat_i)
  }
)

res_bcm_model_compare_01 <- map_dfr(res_bcm_model_selection_01, "compare")
res_bcm_model_status_01 <- map_dfr(res_bcm_model_selection_01, "status")

res_bcm_model_best_01 <- res_bcm_model_compare_01 |>
  filter(!is.na(best_model)) |>
  group_by(predictor) |>
  slice_head(n = 1) |>
  ungroup() |>
  select(predictor, best_model, AIC, BIC, any_of("Performance_Score"))

dir.create(file.path("output", "tables", "tis"), recursive = TRUE, showWarnings = FALSE)

export(
  list(
    model_compare = res_bcm_model_compare_01,
    model_best = res_bcm_model_best_01,
    model_status = res_bcm_model_status_01
  ),
  file.path(
    "output",
    "tables",
    "tis",
    paste0(format(Sys.Date(), "%y%m%d"), "_bcm_model_selection_01.xlsx")
  )
)

# **********************************************************************
## Model evaluation ----
# **********************************************************************
### Regression table ----
# **********************************************************************
regression_xlsx_file_01 <- file.path(
  "output",
  "tables",
  "tis",
  paste0(format(Sys.Date(), "%y%m%d"), "_bcm_gam_regression_tables_01.xlsx")
)

fit_bcm_gam_table_01 <- function(predictor_col, predictor_label) {
  dat_i <- d18_tis_bcm_model_base_01 |>
    transmute(
      immet_id,
      exam,
      bcm,
      sex,
      age_m0,
      predictor_value = as.numeric(.data[[predictor_col]])
    ) |>
    filter(
      !is.na(bcm),
      !is.na(predictor_value),
      !is.na(age_m0),
      !is.na(sex),
      !is.na(exam),
      !is.na(immet_id)
    ) |>
    droplevels()

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(data.frame(
      section = "meta",
      term = "model",
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      note = "Package mgcv unavailable",
      stringsAsFactors = FALSE
    ))
  }

  if (nrow(dat_i) < 20 || dplyr::n_distinct(dat_i$immet_id) < 5 || dplyr::n_distinct(dat_i$predictor_value) < 3) {
    return(data.frame(
      section = "meta",
      term = "model",
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      note = "Insufficient data for GAM model",
      stringsAsFactors = FALSE
    ))
  }

  k_use <- max(3, min(5, dplyr::n_distinct(dat_i$predictor_value) - 1))
  gam_model_obj <- tryCatch(
    mgcv::gam(
      bcm ~ s(predictor_value, k = k_use) + exam + s(age_m0, k = 4) + sex + s(immet_id, bs = "re"),
      data = dat_i,
      method = "REML"
    ),
    error = function(e) NULL
  )

  if (is.null(gam_model_obj)) {
    return(data.frame(
      section = "meta",
      term = "model",
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      note = "GAM model fit failed",
      stringsAsFactors = FALSE
    ))
  }

  sm <- summary(gam_model_obj)

  param_tbl <- tryCatch(
    {
      pt <- as.data.frame(sm$p.table, stringsAsFactors = FALSE)
      pt$term <- rownames(pt)
      rownames(pt) <- NULL
      data.frame(
        section = "parametric",
        term = pt$term,
        estimate = as.numeric(pt[[1]]),
        std_error = if (ncol(pt) >= 2) as.numeric(pt[[2]]) else NA_real_,
        statistic = if (ncol(pt) >= 3) as.numeric(pt[[3]]) else NA_real_,
        p_value = if (ncol(pt) >= 4) as.numeric(pt[[4]]) else NA_real_,
        note = NA_character_,
        stringsAsFactors = FALSE
      )
    },
    error = function(e) data.frame()
  )

  smooth_tbl <- tryCatch(
    {
      st <- as.data.frame(sm$s.table, stringsAsFactors = FALSE)
      st$term <- rownames(st)
      rownames(st) <- NULL
      data.frame(
        section = "smooth",
        term = st$term,
        estimate = if (ncol(st) >= 1) as.numeric(st[[1]]) else NA_real_,
        std_error = if (ncol(st) >= 2) as.numeric(st[[2]]) else NA_real_,
        statistic = if (ncol(st) >= 3) as.numeric(st[[3]]) else NA_real_,
        p_value = if (ncol(st) >= 4) as.numeric(st[[4]]) else NA_real_,
        note = NA_character_,
        stringsAsFactors = FALSE
      )
    },
    error = function(e) data.frame()
  )

  meta_tbl <- data.frame(
    section = "meta",
    term = c("predictor", "formula", "n_obs", "AIC", "BIC", "adj_r_squared", "deviance_explained"),
    estimate = c(
      NA_real_,
      NA_real_,
      as.numeric(stats::nobs(gam_model_obj)),
      tryCatch(as.numeric(AIC(gam_model_obj)), error = function(e) NA_real_),
      tryCatch(as.numeric(BIC(gam_model_obj)), error = function(e) NA_real_),
      if (!is.null(sm$r.sq) && length(sm$r.sq) == 1) as.numeric(sm$r.sq) else NA_real_,
      if (!is.null(sm$dev.expl) && length(sm$dev.expl) == 1) as.numeric(sm$dev.expl) else NA_real_
    ),
    std_error = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    note = c(
      predictor_label,
      paste(deparse(stats::formula(gam_model_obj)), collapse = ""),
      NA_character_,
      NA_character_,
      NA_character_,
      NA_character_,
      NA_character_
    ),
    stringsAsFactors = FALSE
  )

  out_tbl <- dplyr::bind_rows(meta_tbl, param_tbl, smooth_tbl)

  if (nrow(out_tbl) == 0) {
    out_tbl <- data.frame(
      section = "meta",
      term = "model",
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      note = "No rows returned from GAM summary",
      stringsAsFactors = FALSE
    )
  }

  out_tbl
}

fit_bcm_gam_model_01 <- function(predictor_col) {
  dat_i <- d18_tis_bcm_model_base_01 |>
    transmute(
      immet_id,
      exam,
      bcm,
      sex,
      age_m0,
      predictor_value = as.numeric(.data[[predictor_col]])
    ) |>
    filter(
      !is.na(bcm),
      !is.na(predictor_value),
      !is.na(age_m0),
      !is.na(sex),
      !is.na(exam),
      !is.na(immet_id)
    ) |>
    droplevels()

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(NULL)
  }

  if (nrow(dat_i) < 20 || dplyr::n_distinct(dat_i$immet_id) < 5 || dplyr::n_distinct(dat_i$predictor_value) < 3) {
    return(NULL)
  }

  k_use <- max(3, min(5, dplyr::n_distinct(dat_i$predictor_value) - 1))

  tryCatch(
    mgcv::gam(
      bcm ~ s(predictor_value, k = k_use) + exam + s(age_m0, k = 4) + sex + s(immet_id, bs = "re"),
      data = dat_i,
      method = "REML"
    ),
    error = function(e) NULL
  )
}

res_bcm_gam_regression_tables_01 <- imap(
  bcm_predictor_candidates_01,
  \(predictor_col, nm) {
    fit_bcm_gam_table_01(
      predictor_col = predictor_col,
      predictor_label = predictor_model_label_map_01[[nm]]
    )
  }
)

names(res_bcm_gam_regression_tables_01) <- names(res_bcm_gam_regression_tables_01) |>
  stringr::str_replace_all("[^A-Za-z0-9]", "_") |>
  stringr::str_sub(1, 31)

if (file.exists(regression_xlsx_file_01)) {
  removed_old_01 <- isTRUE(file.remove(regression_xlsx_file_01))
  if (!removed_old_01) {
    warning(
      paste0(
        "Existing file could not be removed before overwrite: ",
        regression_xlsx_file_01
      )
    )
  }
}

rio::export(
  x = stats::setNames(res_bcm_gam_regression_tables_01, names(res_bcm_gam_regression_tables_01)),
  file = regression_xlsx_file_01
)

# **********************************************************************
### Model plots ----
# **********************************************************************
fit_bcm_gam_model_plot_01 <- function(predictor_col) {
  dat_i <- d18_tis_bcm_model_base_01 |>
    transmute(
      immet_id,
      exam,
      bcm,
      sex,
      age_m0,
      predictor_value = as.numeric(.data[[predictor_col]])
    ) |>
    filter(
      !is.na(bcm),
      !is.na(predictor_value),
      !is.na(age_m0),
      !is.na(sex),
      !is.na(exam),
      !is.na(immet_id)
    ) |>
    droplevels()

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(NULL)
  }

  if (nrow(dat_i) < 20 || dplyr::n_distinct(dat_i$immet_id) < 5 || dplyr::n_distinct(dat_i$predictor_value) < 3) {
    return(NULL)
  }

  k_use <- max(3, min(5, dplyr::n_distinct(dat_i$predictor_value) - 1))

  tryCatch(
    mgcv::gam(
      bcm ~ s(predictor_value, k = k_use) + exam + s(age_m0, k = 4) + sex + s(immet_id, bs = "re"),
      data = dat_i,
      method = "REML"
    ),
    error = function(e) NULL
  )
}

res_bcm_gam_plot_manifest_01 <- imap_dfr(
  bcm_predictor_candidates_01,
  \(predictor_col, nm) {
    predictor_label <- predictor_model_label_map_01[[nm]]
    mod <- fit_bcm_gam_model_plot_01(predictor_col)
    file_nm <- file.path(
      "output",
      "figures",
      "tis",
      paste0(
        format(Sys.Date(), "%y%m%d"),
        "_bcm_gam_plot_",
        nm,
        "_01.png"
      )
    )

    saved <- FALSE

    if (!is.null(mod) && requireNamespace("ggeffects", quietly = TRUE)) {
      pred_obj <- tryCatch(
        ggeffects::ggpredict(mod, terms = c("predictor_value", "exam"), ci_level = 0.95),
        error = function(e) NULL
      )

      plot_obj <- NULL
      if (!is.null(pred_obj)) {
        plot_obj <- tryCatch(
          plot(pred_obj, show_data = FALSE, show_residuals = FALSE, colors = "bw") +
            labs(
              y = "BCM",
              x = predictor_label,
              title = ""
            ) +
            theme_bw(base_size = 16) +
            theme(
              axis.title = element_text(size = 20, face = "bold", colour = "black"),
              axis.text = element_text(size = 15, colour = "black"),
              legend.title = element_text(size = 16, face = "bold", colour = "black"),
              legend.text = element_text(size = 14, colour = "black")
            ),
          error = function(e) NULL
        )
      }

      if (!is.null(plot_obj)) {
        if (file.exists(file_nm)) {
          removed_old_01 <- isTRUE(file.remove(file_nm))
          if (!removed_old_01) {
            warning(
              paste0(
                "Existing file could not be removed before overwrite: ",
                file_nm
              )
            )
          }
        }

        saved <- tryCatch(
          {
            ggsave(
              filename = file_nm,
              plot = plot_obj,
              width = 12,
              height = 8,
              units = "in",
              dpi = 300,
              bg = "white",
              limitsize = FALSE
            )
            TRUE
          },
          error = function(e) FALSE
        )
      }
    }

    tibble(
      predictor = predictor_label,
      file_path = file_nm,
      saved = saved
    )
  }
)

# **********************************************************************
### Validation of TIS and activity parameters ----
# **********************************************************************
validation_predictors_01 <- c(
  tis_total = "tis_total",
  unname(activity_col_map_01)
)
names(validation_predictors_01) <- c("tis_total", names(activity_col_map_01))

validation_label_map_01 <- c(
  tis_total = "TIS total",
  activity_label_map_01
)

d18_tis_validation_01 <- d18_ggir_tis_data_02 |>
  select(
    immet_id,
    exam,
    sex,
    mda_category,
    phase_angle,
    bcm,
    tis_total,
    any_of(unname(activity_col_map_01))
  ) |>
  mutate(
    mda_category_num = suppressWarnings(as.numeric(as.character(mda_category))),
    phase_angle = as.numeric(phase_angle),
    bcm = as.numeric(bcm)
  ) |>
  filter(!is.na(exam)) |>
  mutate(
    mda_category_factor = if_else(
      is.na(mda_category_num),
      NA_character_,
      as.character(mda_category_num)
    ),
    mda_category_factor = factor(
      mda_category_factor,
      levels = sort(unique(mda_category_factor[!is.na(mda_category_factor)]))
    )
  )

mda_min_level_01 <- suppressWarnings(min(d18_tis_validation_01$mda_category_num, na.rm = TRUE))
if (!is.finite(mda_min_level_01)) {
  mda_min_level_01 <- NA_real_
}

d18_tis_validation_01 <- d18_tis_validation_01 |>
  mutate(
    mda_active_bin = case_when(
      is.na(mda_category_num) ~ NA_integer_,
      is.na(mda_min_level_01) ~ NA_integer_,
      mda_category_num > mda_min_level_01 ~ 1L,
      TRUE ~ 0L
    )
  )

fit_validation_gam_binomial_01 <- function(dat) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(NULL)
  }

  if (nrow(dat) < 20 || dplyr::n_distinct(dat$immet_id) < 5 || dplyr::n_distinct(dat$predictor_value) < 3 || dplyr::n_distinct(dat$mda_active_bin) < 2) {
    return(NULL)
  }

  k_use <- max(3, min(5, dplyr::n_distinct(dat$predictor_value) - 1))

  tryCatch(
    mgcv::gam(
      mda_active_bin ~ s(predictor_value, k = k_use) + sex + exam + s(immet_id, bs = "re"),
      data = dat,
      family = binomial(),
      method = "REML"
    ),
    error = function(e) NULL
  )
}

fit_validation_gam_gaussian_01 <- function(dat, outcome_col) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(NULL)
  }

  if (nrow(dat) < 20 || dplyr::n_distinct(dat$immet_id) < 5 || dplyr::n_distinct(dat$predictor_value) < 3) {
    return(NULL)
  }

  k_use <- max(3, min(5, dplyr::n_distinct(dat$predictor_value) - 1))
  form <- as.formula(
    paste0(outcome_col, " ~ s(predictor_value, k = ", k_use, ") + sex + exam + s(immet_id, bs = 're')")
  )

  tryCatch(
    mgcv::gam(
      form,
      data = dat,
      method = "REML"
    ),
    error = function(e) NULL
  )
}

calc_roc_from_model_01 <- function(model_obj, dat) {
  if (is.null(model_obj)) {
    return(tibble(
      auc = NA_real_,
      threshold = NA_real_,
      sensitivity = NA_real_,
      specificity = NA_real_,
      accuracy = NA_real_,
      n = nrow(dat),
      note = "GAM binomial fit failed"
    ))
  }

  if (!requireNamespace("pROC", quietly = TRUE)) {
    return(tibble(
      auc = NA_real_,
      threshold = NA_real_,
      sensitivity = NA_real_,
      specificity = NA_real_,
      accuracy = NA_real_,
      n = nrow(dat),
      note = "Package pROC unavailable"
    ))
  }

  pred_prob <- tryCatch(
    predict(model_obj, newdata = dat, type = "response"),
    error = function(e) rep(NA_real_, nrow(dat))
  )

  dat_eval <- dat |>
    mutate(pred_prob = pred_prob) |>
    filter(!is.na(pred_prob), !is.na(mda_active_bin))

  if (nrow(dat_eval) < 10 || n_distinct(dat_eval$mda_active_bin) < 2) {
    return(tibble(
      auc = NA_real_,
      threshold = NA_real_,
      sensitivity = NA_real_,
      specificity = NA_real_,
      accuracy = NA_real_,
      n = nrow(dat_eval),
      note = "Insufficient data for ROC"
    ))
  }

  roc_obj <- tryCatch(
    pROC::roc(dat_eval$mda_active_bin, dat_eval$pred_prob, quiet = TRUE),
    error = function(e) NULL
  )

  if (is.null(roc_obj)) {
    return(tibble(
      auc = NA_real_,
      threshold = NA_real_,
      sensitivity = NA_real_,
      specificity = NA_real_,
      accuracy = NA_real_,
      n = nrow(dat_eval),
      note = "ROC fit failed"
    ))
  }

  coords_best <- tryCatch(
    pROC::coords(
      roc_obj,
      x = "best",
      best.method = "youden",
      ret = c("threshold", "sensitivity", "specificity")
    ),
    error = function(e) c(threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_)
  )

  thr <- as.numeric(coords_best[["threshold"]])
  pred_class <- ifelse(dat_eval$pred_prob >= thr, 1L, 0L)
  acc <- mean(pred_class == dat_eval$mda_active_bin)

  tibble(
    auc = as.numeric(pROC::auc(roc_obj)),
    threshold = thr,
    sensitivity = as.numeric(coords_best[["sensitivity"]]),
    specificity = as.numeric(coords_best[["specificity"]]),
    accuracy = acc,
    n = nrow(dat_eval),
    note = NA_character_
  )
}

calc_continuous_metrics_01 <- function(model_obj, dat, outcome_col) {
  if (is.null(model_obj)) {
    return(tibble(
      r2_marginal = NA_real_,
      r2_conditional = NA_real_,
      rmse = NA_real_,
      mae = NA_real_,
      n = nrow(dat),
      note = "GAM gaussian fit failed"
    ))
  }

  pred_val <- tryCatch(
    predict(model_obj, newdata = dat, type = "response"),
    error = function(e) rep(NA_real_, nrow(dat))
  )

  obs_val <- dat[[outcome_col]]
  keep <- is.finite(obs_val) & is.finite(pred_val)
  obs_val <- obs_val[keep]
  pred_val <- pred_val[keep]

  r2_obj <- tryCatch(performance::r2(model_obj), error = function(e) NULL)
  r2_marginal <- if (!is.null(r2_obj) && "R2" %in% names(r2_obj)) as.numeric(r2_obj$R2[[1]]) else if (!is.null(r2_obj) && "R2_marginal" %in% names(r2_obj)) as.numeric(r2_obj$R2_marginal[[1]]) else NA_real_
  r2_conditional <- if (!is.null(r2_obj) && "R2_conditional" %in% names(r2_obj)) as.numeric(r2_obj$R2_conditional[[1]]) else NA_real_

  tibble(
    r2_marginal = r2_marginal,
    r2_conditional = r2_conditional,
    rmse = if (length(obs_val) > 0) sqrt(mean((obs_val - pred_val)^2)) else NA_real_,
    mae = if (length(obs_val) > 0) mean(abs(obs_val - pred_val)) else NA_real_,
    n = length(obs_val),
    note = NA_character_
  )
}

res_validation_anchor_overview_01 <- d18_tis_validation_01 |>
  summarise(
    n_rows = n(),
    n_patients = n_distinct(immet_id),
    n_mda_non_missing = sum(!is.na(mda_category_num)),
    n_phase_angle_non_missing = sum(!is.na(phase_angle)),
    n_bcm_non_missing = sum(!is.na(bcm)),
    mda_levels = paste(sort(unique(mda_category_num[!is.na(mda_category_num)])), collapse = ", ")
  )

res_validation_mda_roc_01 <- imap_dfr(
  validation_predictors_01,
  \(predictor_col, nm) {
    dat_i <- d18_tis_validation_01 |>
      transmute(
        immet_id = factor(immet_id),
        exam = factor(exam, levels = exam_keep),
        sex = factor(sex),
        mda_active_bin = as.integer(mda_active_bin),
        predictor_value = as.numeric(.data[[predictor_col]])
      ) |>
      filter(!is.na(immet_id), !is.na(exam), !is.na(sex), !is.na(mda_active_bin), !is.na(predictor_value))

    fit_validation_gam_binomial_01(dat_i) |>
      (\(mod) calc_roc_from_model_01(mod, dat_i))() |>
      mutate(
        predictor = validation_label_map_01[[nm]],
        .before = 1
      )
  }
)

res_validation_phase_angle_model_01 <- imap_dfr(
  validation_predictors_01,
  \(predictor_col, nm) {
    dat_i <- d18_tis_validation_01 |>
      transmute(
        immet_id = factor(immet_id),
        exam = factor(exam, levels = exam_keep),
        sex = factor(sex),
        phase_angle = as.numeric(phase_angle),
        predictor_value = as.numeric(.data[[predictor_col]])
      ) |>
      filter(!is.na(immet_id), !is.na(exam), !is.na(sex), !is.na(phase_angle), !is.na(predictor_value))

    fit_validation_gam_gaussian_01(dat_i, "phase_angle") |>
      (\(mod) calc_continuous_metrics_01(mod, dat_i, "phase_angle"))() |>
      mutate(
        predictor = validation_label_map_01[[nm]],
        anchor = "phase_angle",
        .before = 1
      )
  }
)

res_validation_bcm_model_01 <- imap_dfr(
  validation_predictors_01,
  \(predictor_col, nm) {
    dat_i <- d18_tis_validation_01 |>
      transmute(
        immet_id = factor(immet_id),
        exam = factor(exam, levels = exam_keep),
        sex = factor(sex),
        bcm = as.numeric(bcm),
        predictor_value = as.numeric(.data[[predictor_col]])
      ) |>
      filter(!is.na(immet_id), !is.na(exam), !is.na(sex), !is.na(bcm), !is.na(predictor_value))

    fit_validation_gam_gaussian_01(dat_i, "bcm") |>
      (\(mod) calc_continuous_metrics_01(mod, dat_i, "bcm"))() |>
      mutate(
        predictor = validation_label_map_01[[nm]],
        anchor = "bcm",
        .before = 1
      )
  }
)

validation_xlsx_file_01 <- file.path(
  "output",
  "tables",
  "tis",
  paste0(format(Sys.Date(), "%y%m%d"), "_tis_activity_validation_01.xlsx")
)

if (file.exists(validation_xlsx_file_01)) {
  removed_old_01 <- isTRUE(file.remove(validation_xlsx_file_01))
  if (!removed_old_01) {
    warning(
      paste0(
        "Existing file could not be removed before overwrite: ",
        validation_xlsx_file_01
      )
    )
  }
}

rio::export(
  x = list(
    anchor_overview = as.data.frame(res_validation_anchor_overview_01),
    mda_roc = as.data.frame(res_validation_mda_roc_01),
    phase_angle_model = as.data.frame(res_validation_phase_angle_model_01),
    bcm_model = as.data.frame(res_validation_bcm_model_01)
  ),
  file = validation_xlsx_file_01
)

# **********************************************************************
## Model evaluation 2 ----
# **********************************************************************
### Regression table ----
# **********************************************************************

dir.create(file.path("output", "tables", "tis"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("output", "figures", "tis"), recursive = TRUE, showWarnings = FALSE)

lmer_regression_xlsx_file_01 <- file.path(
  "output",
  "tables",
  "tis",
  paste0(format(Sys.Date(), "%y%m%d"), "_bcm_lmer_regression_tables_01.xlsx")
)

fit_bcm_lmer_table_01 <- function(predictor_col, predictor_label) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(data.frame(
      section = "meta",
      term = "model",
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      note = "Package lme4 unavailable"
    ))
  }

  dat_i <- d18_tis_bcm_model_base_01 |>
    transmute(
      immet_id = factor(immet_id),
      exam = factor(exam, levels = exam_keep),
      bcm = as.numeric(bcm),
      sex = factor(sex),
      age_m0 = as.numeric(age_m0),
      predictor_value = as.numeric(.data[[predictor_col]])
    ) |>
    filter(
      !is.na(immet_id),
      !is.na(exam),
      !is.na(bcm),
      !is.na(sex),
      !is.na(age_m0),
      !is.na(predictor_value)
    )

  if (nrow(dat_i) < 20 || dplyr::n_distinct(dat_i$immet_id) < 5 || dplyr::n_distinct(dat_i$predictor_value) < 3) {
    return(data.frame(
      section = "meta",
      term = "model",
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      note = "Insufficient data for LMER model"
    ))
  }

  mod <- tryCatch(
    {
      if (requireNamespace("lmerTest", quietly = TRUE)) {
        lmerTest::lmer(
          bcm ~ predictor_value + exam + scale(age_m0) + sex + (1 | immet_id),
          data = dat_i,
          REML = TRUE
        )
      } else {
        lme4::lmer(
          bcm ~ predictor_value + exam + scale(age_m0) + sex + (1 | immet_id),
          data = dat_i,
          REML = TRUE
        )
      }
    },
    error = function(e) NULL
  )

  if (is.null(mod)) {
    return(data.frame(
      section = "meta",
      term = "model",
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      note = "LMER model fit failed"
    ))
  }

  coef_tbl <- tryCatch(
    {
      if (requireNamespace("broom.mixed", quietly = TRUE)) {
        coef_tidy <- broom.mixed::tidy(
          mod,
          effects = "fixed",
          conf.int = TRUE,
          conf.method = "Wald"
        ) |>
          as_tibble()

        coef_tidy |>
          transmute(
            section = "fixed_effects",
            term = term,
            estimate = estimate,
            std_error = if ("std.error" %in% names(coef_tidy)) .data[["std.error"]] else NA_real_,
            statistic = if ("statistic" %in% names(coef_tidy)) .data[["statistic"]] else if ("t.value" %in% names(coef_tidy)) .data[["t.value"]] else NA_real_,
            p_value = if ("p.value" %in% names(coef_tidy)) .data[["p.value"]] else NA_real_,
            conf_low = if ("conf.low" %in% names(coef_tidy)) .data[["conf.low"]] else NA_real_,
            conf_high = if ("conf.high" %in% names(coef_tidy)) .data[["conf.high"]] else NA_real_
          )
      } else {
        sm <- summary(mod)
        coef_mat <- as.data.frame(sm$coefficients)
        coef_mat$term <- rownames(coef_mat)
        rownames(coef_mat) <- NULL

        ci_mat <- tryCatch(
          suppressMessages(confint(mod, parm = "beta_", method = "Wald")),
          error = function(e) NULL
        )

        coef_tbl_fallback <- coef_mat |>
          as_tibble() |>
          mutate(
            estimate = .data[["Estimate"]],
            std_error = .data[["Std. Error"]],
            statistic = dplyr::coalesce(
              if ("t value" %in% names(coef_mat)) .data[["t value"]] else NA_real_,
              if ("z value" %in% names(coef_mat)) .data[["z value"]] else NA_real_
            ),
            p_value = if ("Pr(>|t|)" %in% names(.)) .data[["Pr(>|t|)"]] else if ("Pr(>|z|)" %in% names(.)) .data[["Pr(>|z|)"]] else NA_real_
          ) |>
          select(term, estimate, std_error, statistic, p_value)

        if (!is.null(ci_mat)) {
          ci_tbl <- as.data.frame(ci_mat)
          ci_tbl$term <- rownames(ci_tbl)
          rownames(ci_tbl) <- NULL
          names(ci_tbl)[1:2] <- c("conf_low", "conf_high")
          coef_tbl_fallback <- coef_tbl_fallback |>
            left_join(ci_tbl, by = "term")
        } else {
          coef_tbl_fallback <- coef_tbl_fallback |>
            mutate(conf_low = NA_real_, conf_high = NA_real_)
        }

        coef_tbl_fallback |>
          mutate(section = "fixed_effects", .before = 1)
      }
    },
    error = function(e) tibble()
  )

  re_var <- tryCatch(as.data.frame(lme4::VarCorr(mod)), error = function(e) NULL)
  re_tbl <- if (!is.null(re_var) && nrow(re_var) > 0) {
    re_var |>
      transmute(
        section = "random_effects",
        term = paste(grp, var1, sep = ":"),
        estimate = vcov,
        std_error = sdcor,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_
      )
  } else {
    tibble()
  }

  r2_obj <- tryCatch(performance::r2(mod), error = function(e) NULL)
  meta_tbl <- tibble(
    section = "meta",
    term = c("predictor", "formula", "n_obs", "AIC", "BIC", "R2_marginal", "R2_conditional"),
    estimate = c(
      NA_real_,
      NA_real_,
      nobs(mod),
      tryCatch(AIC(mod), error = function(e) NA_real_),
      tryCatch(BIC(mod), error = function(e) NA_real_),
      if (!is.null(r2_obj) && "R2_marginal" %in% names(r2_obj)) as.numeric(r2_obj$R2_marginal[[1]]) else NA_real_,
      if (!is.null(r2_obj) && "R2_conditional" %in% names(r2_obj)) as.numeric(r2_obj$R2_conditional[[1]]) else NA_real_
    ),
    std_error = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    note = c(
      predictor_label,
      paste(deparse(stats::formula(mod)), collapse = ""),
      NA_character_,
      NA_character_,
      NA_character_,
      NA_character_,
      NA_character_
    )
  )

  bind_rows(
    meta_tbl,
    coef_tbl |> mutate(note = NA_character_),
    re_tbl |> mutate(note = NA_character_)
  ) |>
    as.data.frame()
}

res_bcm_lmer_regression_tables_01 <- imap(
  bcm_predictor_candidates_01,
  \(predictor_col, nm) {
    fit_bcm_lmer_table_01(
      predictor_col = predictor_col,
      predictor_label = predictor_model_label_map_01[[nm]]
    )
  }
)

names(res_bcm_lmer_regression_tables_01) <- names(res_bcm_lmer_regression_tables_01) |>
  stringr::str_replace_all("[^A-Za-z0-9]", "_") |>
  stringr::str_sub(1, 31)

if (file.exists(lmer_regression_xlsx_file_01)) {
  removed_old_01 <- isTRUE(file.remove(lmer_regression_xlsx_file_01))
  if (!removed_old_01) {
    warning(
      paste0(
        "Existing file could not be removed before overwrite: ",
        lmer_regression_xlsx_file_01
      )
    )
  }
}

rio::export(
  x = stats::setNames(res_bcm_lmer_regression_tables_01, names(res_bcm_lmer_regression_tables_01)),
  file = lmer_regression_xlsx_file_01
)

# **********************************************************************
### Model plots ----
# **********************************************************************

fit_bcm_lmer_model_plot_01 <- function(predictor_col) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(NULL)
  }

  dat_i <- d18_tis_bcm_model_base_01 |>
    transmute(
      immet_id = factor(immet_id),
      exam = factor(exam, levels = exam_keep),
      bcm = as.numeric(bcm),
      sex = factor(sex),
      age_m0 = as.numeric(age_m0),
      predictor_value = as.numeric(.data[[predictor_col]])
    ) |>
    filter(
      !is.na(immet_id),
      !is.na(exam),
      !is.na(bcm),
      !is.na(sex),
      !is.na(age_m0),
      !is.na(predictor_value)
    )

  if (nrow(dat_i) < 20 || dplyr::n_distinct(dat_i$immet_id) < 5 || dplyr::n_distinct(dat_i$predictor_value) < 3) {
    return(NULL)
  }

  tryCatch(
    lme4::lmer(
      bcm ~ predictor_value + exam + scale(age_m0) + sex + (1 | immet_id),
      data = dat_i,
      REML = TRUE
    ),
    error = function(e) NULL
  )
}

res_bcm_lmer_plot_manifest_01 <- imap_dfr(
  bcm_predictor_candidates_01,
  \(predictor_col, nm) {
    predictor_label <- predictor_model_label_map_01[[nm]]
    mod <- fit_bcm_lmer_model_plot_01(predictor_col)
    file_nm <- file.path(
      "output",
      "figures",
      "tis",
      paste0(
        format(Sys.Date(), "%y%m%d"),
        "_bcm_lmer_plot_",
        nm,
        "_01.png"
      )
    )

    saved <- FALSE

    if (!is.null(mod) && requireNamespace("ggeffects", quietly = TRUE)) {
      pred_obj <- tryCatch(
        ggeffects::ggpredict(mod, terms = c("predictor_value", "exam"), ci_level = 0.95),
        error = function(e) NULL
      )

      plot_obj <- NULL
      if (!is.null(pred_obj)) {
        plot_obj <- tryCatch(
          plot(pred_obj, show_data = FALSE, show_residuals = FALSE, colors = "bw") +
            labs(
              y = "BCM",
              x = predictor_label,
              title = ""
            ) +
            theme_bw(base_size = 16) +
            theme(
              axis.title = element_text(size = 20, face = "bold", colour = "black"),
              axis.text = element_text(size = 15, colour = "black"),
              legend.title = element_text(size = 16, face = "bold", colour = "black"),
              legend.text = element_text(size = 14, colour = "black")
            ),
          error = function(e) NULL
        )
      }

      if (!is.null(plot_obj)) {
        if (file.exists(file_nm)) {
          removed_old_01 <- isTRUE(file.remove(file_nm))
          if (!removed_old_01) {
            warning(
              paste0(
                "Existing file could not be removed before overwrite: ",
                file_nm
              )
            )
          }
        }

        saved <- tryCatch(
          {
            ggsave(
              filename = file_nm,
              plot = plot_obj,
              width = 12,
              height = 8,
              units = "in",
              dpi = 300,
              bg = "white",
              limitsize = FALSE
            )
            TRUE
          },
          error = function(e) FALSE
        )
      }
    }

    tibble(
      predictor = predictor_label,
      file_path = file_nm,
      saved = saved
    )
  }
)

# **********************************************************************
## Model comparison: TIS + activity parameters ----
# **********************************************************************

activity_only_candidates_01 <- unname(activity_col_map_01)
names(activity_only_candidates_01) <- names(activity_col_map_01)

make_tis_activity_model_labels_01 <- function(activity_label) {
  c(
    model_a = "Model 1: TIS",
    model_b = paste0("Model 2: TIS + ", activity_label),
    model_c = paste0("Model 3: ", activity_label)
  )
}

relabel_compare_performance_tis_01 <- function(perf_obj, model_labels) {
  if (is.null(perf_obj)) {
    return(perf_obj)
  }

  nm_col <- intersect(c("Name", "name", "Model", "model"), names(perf_obj))
  if (length(nm_col) == 0) {
    return(perf_obj)
  }
  nm_col <- nm_col[[1]]

  old_vals <- as.character(perf_obj[[nm_col]])
  idx <- suppressWarnings(as.integer(stringr::str_extract(old_vals, "\\d+")))
  valid_idx <- !is.na(idx) & idx >= 1 & idx <= length(model_labels)

  new_vals <- old_vals
  new_vals[valid_idx] <- model_labels[idx[valid_idx]]
  perf_obj[[nm_col]] <- new_vals

  chr_cols <- names(perf_obj)[vapply(perf_obj, is.character, logical(1))]
  for (cc in chr_cols) {
    vals_i <- as.character(perf_obj[[cc]])
    vals_i[valid_idx] <- new_vals[valid_idx]
    perf_obj[[cc]] <- vals_i
  }

  perf_obj
}

fit_bcm_compare_models_lmer_01 <- function(dat) {
  out <- list(model_a = NULL, model_b = NULL, model_c = NULL)

  dat <- dat |>
    filter(
      !is.na(bcm),
      !is.na(tis_total),
      !is.na(activity_value),
      !is.na(age_m0),
      !is.na(sex),
      !is.na(exam),
      !is.na(immet_id)
    ) |>
    droplevels()

  if (nrow(dat) < 20 || dplyr::n_distinct(dat$immet_id) < 5 || dplyr::n_distinct(dat$tis_total) < 3 || dplyr::n_distinct(dat$activity_value) < 3) {
    return(out)
  }

  out$model_a <- tryCatch(
    lme4::lmer(
      bcm ~ scale(tis_total) + exam + scale(age_m0) + sex + (1 | immet_id),
      data = dat,
      REML = FALSE
    ),
    error = function(e) NULL
  )

  out$model_b <- tryCatch(
    lme4::lmer(
      bcm ~ scale(tis_total) + scale(activity_value) + exam + scale(age_m0) + sex + (1 | immet_id),
      data = dat,
      REML = FALSE
    ),
    error = function(e) NULL
  )

  out$model_c <- tryCatch(
    lme4::lmer(
      bcm ~ scale(activity_value) + exam + scale(age_m0) + sex + (1 | immet_id),
      data = dat,
      REML = FALSE
    ),
    error = function(e) NULL
  )

  out
}

fit_bcm_compare_models_gam_01 <- function(dat) {
  out <- list(model_a = NULL, model_b = NULL, model_c = NULL)

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(out)
  }

  dat <- dat |>
    filter(
      !is.na(bcm),
      !is.na(tis_total),
      !is.na(activity_value),
      !is.na(age_m0),
      !is.na(sex),
      !is.na(exam),
      !is.na(immet_id)
    ) |>
    droplevels()

  if (nrow(dat) < 20 || dplyr::n_distinct(dat$immet_id) < 5 || dplyr::n_distinct(dat$tis_total) < 3 || dplyr::n_distinct(dat$activity_value) < 3) {
    return(out)
  }

  k_tis <- max(3, min(5, dplyr::n_distinct(dat$tis_total) - 1))
  k_act <- max(3, min(5, dplyr::n_distinct(dat$activity_value) - 1))

  out$model_a <- tryCatch(
    mgcv::gam(
      bcm ~ s(tis_total, k = k_tis) + exam + s(age_m0, k = 4) + sex + s(immet_id, bs = "re"),
      data = dat,
      method = "REML"
    ),
    error = function(e) NULL
  )

  out$model_b <- tryCatch(
    mgcv::gam(
      bcm ~ s(tis_total, k = k_tis) + s(activity_value, k = k_act) + exam + s(age_m0, k = 4) + sex + s(immet_id, bs = "re"),
      data = dat,
      method = "REML"
    ),
    error = function(e) NULL
  )

  out$model_c <- tryCatch(
    mgcv::gam(
      bcm ~ s(activity_value, k = k_act) + exam + s(age_m0, k = 4) + sex + s(immet_id, bs = "re"),
      data = dat,
      method = "REML"
    ),
    error = function(e) NULL
  )

  out
}

summarise_tis_activity_compare_01 <- function(models_named, activity_label) {
  if (length(models_named) == 0) {
    return(tibble(
      activity_parameter = activity_label,
      model = NA_character_,
      AIC = NA_real_,
      BIC = NA_real_,
      Performance_Score = NA_real_,
      best_model = NA_character_,
      note = "No valid models"
    ))
  }

  perf_tbl <- NULL
  if (length(models_named) >= 2 && requireNamespace("performance", quietly = TRUE)) {
    perf_tbl <- tryCatch(
      do.call(performance::compare_performance, models_named) |>
        relabel_compare_performance_tis_01(model_labels = names(models_named)) |>
        as_tibble(),
      error = function(e) NULL
    )
  }

  if (is.null(perf_tbl)) {
    perf_tbl <- imap_dfr(
      models_named,
      \(mod, nm) tibble(
        model = nm,
        AIC = tryCatch(AIC(mod), error = function(e) NA_real_),
        BIC = tryCatch(BIC(mod), error = function(e) NA_real_)
      )
    )
  }

  name_col <- intersect(c("Name", "Model", "model", "name"), names(perf_tbl))
  if (length(name_col) > 0) {
    perf_tbl <- perf_tbl |>
      rename(model = all_of(name_col[[1]]))
  }

  if (!("AIC" %in% names(perf_tbl))) {
    perf_tbl <- perf_tbl |>
      mutate(AIC = map_dbl(model, ~ tryCatch(AIC(models_named[[.x]]), error = function(e) NA_real_)))
  }

  if (!("BIC" %in% names(perf_tbl))) {
    perf_tbl <- perf_tbl |>
      mutate(BIC = map_dbl(model, ~ tryCatch(BIC(models_named[[.x]]), error = function(e) NA_real_)))
  }

  best_model_name <- NA_character_
  if ("Performance_Score" %in% names(perf_tbl) && any(is.finite(perf_tbl$Performance_Score))) {
    best_model_name <- perf_tbl |>
      filter(is.finite(Performance_Score)) |>
      arrange(desc(Performance_Score), AIC) |>
      slice_head(n = 1) |>
      pull(model)
  } else if (any(is.finite(perf_tbl$AIC))) {
    best_model_name <- perf_tbl |>
      filter(is.finite(AIC)) |>
      arrange(AIC, BIC) |>
      slice_head(n = 1) |>
      pull(model)
  }

  perf_tbl |>
    mutate(
      activity_parameter = activity_label,
      .before = 1
    ) |>
    mutate(
      best_model = best_model_name,
      note = if_else(model == best_model_name, "Selected", NA_character_)
    )
}

res_bcm_compare_lmer_01 <- imap(
  activity_only_candidates_01,
  \(predictor_col, nm) {
    activity_label <- activity_label_map_01[[nm]]

    dat_i <- d18_tis_bcm_model_base_01 |>
      transmute(
        immet_id,
        exam,
        bcm,
        sex,
        age_m0,
        tis_total = as.numeric(tis_total),
        activity_value = as.numeric(.data[[predictor_col]])
      )

    fit_obj <- fit_bcm_compare_models_lmer_01(dat_i)
    model_labels <- make_tis_activity_model_labels_01(activity_label)
    models_named <- list(
      fit_obj$model_a,
      fit_obj$model_b,
      fit_obj$model_c
    )
    names(models_named) <- unname(model_labels[c("model_a", "model_b", "model_c")])
    models_named <- models_named[!map_lgl(models_named, is.null)]

    list(
      compare = summarise_tis_activity_compare_01(models_named, activity_label),
      status = tibble(
        activity_parameter = activity_label,
        n = nrow(dat_i |> filter(!is.na(bcm), !is.na(tis_total), !is.na(activity_value), !is.na(age_m0), !is.na(sex), !is.na(exam))),
        n_patients = dplyr::n_distinct(dat_i$immet_id[!is.na(dat_i$bcm) & !is.na(dat_i$tis_total) & !is.na(dat_i$activity_value)]),
        has_model_a = !is.null(fit_obj$model_a),
        has_model_b = !is.null(fit_obj$model_b),
        has_model_c = !is.null(fit_obj$model_c)
      )
    )
  }
)

res_bcm_compare_gam_01 <- imap(
  activity_only_candidates_01,
  \(predictor_col, nm) {
    activity_label <- activity_label_map_01[[nm]]

    dat_i <- d18_tis_bcm_model_base_01 |>
      transmute(
        immet_id,
        exam,
        bcm,
        sex,
        age_m0,
        tis_total = as.numeric(tis_total),
        activity_value = as.numeric(.data[[predictor_col]])
      )

    fit_obj <- fit_bcm_compare_models_gam_01(dat_i)
    model_labels <- make_tis_activity_model_labels_01(activity_label)
    models_named <- list(
      fit_obj$model_a,
      fit_obj$model_b,
      fit_obj$model_c
    )
    names(models_named) <- unname(model_labels[c("model_a", "model_b", "model_c")])
    models_named <- models_named[!map_lgl(models_named, is.null)]

    list(
      compare = summarise_tis_activity_compare_01(models_named, activity_label),
      status = tibble(
        activity_parameter = activity_label,
        n = nrow(dat_i |> filter(!is.na(bcm), !is.na(tis_total), !is.na(activity_value), !is.na(age_m0), !is.na(sex), !is.na(exam))),
        n_patients = dplyr::n_distinct(dat_i$immet_id[!is.na(dat_i$bcm) & !is.na(dat_i$tis_total) & !is.na(dat_i$activity_value)]),
        has_model_a = !is.null(fit_obj$model_a),
        has_model_b = !is.null(fit_obj$model_b),
        has_model_c = !is.null(fit_obj$model_c)
      )
    )
  }
)

res_bcm_compare_lmer_tbl_01 <- map_dfr(res_bcm_compare_lmer_01, "compare")
res_bcm_compare_lmer_status_01 <- map_dfr(res_bcm_compare_lmer_01, "status")
res_bcm_compare_gam_tbl_01 <- map_dfr(res_bcm_compare_gam_01, "compare")
res_bcm_compare_gam_status_01 <- map_dfr(res_bcm_compare_gam_01, "status")

res_bcm_compare_lmer_best_01 <- res_bcm_compare_lmer_tbl_01 |>
  filter(!is.na(best_model)) |>
  group_by(activity_parameter) |>
  slice_head(n = 1) |>
  ungroup() |>
  select(activity_parameter, best_model, AIC, BIC, any_of("Performance_Score"))

res_bcm_compare_gam_best_01 <- res_bcm_compare_gam_tbl_01 |>
  filter(!is.na(best_model)) |>
  group_by(activity_parameter) |>
  slice_head(n = 1) |>
  ungroup() |>
  select(activity_parameter, best_model, AIC, BIC, any_of("Performance_Score"))

comparison_xlsx_file_01 <- file.path(
  "output",
  "tables",
  "tis",
  paste0(format(Sys.Date(), "%y%m%d"), "_bcm_tis_activity_model_comparison_01.xlsx")
)

if (file.exists(comparison_xlsx_file_01)) {
  removed_old_01 <- isTRUE(file.remove(comparison_xlsx_file_01))
  if (!removed_old_01) {
    warning(
      paste0(
        "Existing file could not be removed before overwrite: ",
        comparison_xlsx_file_01
      )
    )
  }
}

rio::export(
  x = list(
    lmer_compare = as.data.frame(res_bcm_compare_lmer_tbl_01),
    lmer_best = as.data.frame(res_bcm_compare_lmer_best_01),
    lmer_status = as.data.frame(res_bcm_compare_lmer_status_01),
    gam_compare = as.data.frame(res_bcm_compare_gam_tbl_01),
    gam_best = as.data.frame(res_bcm_compare_gam_best_01),
    gam_status = as.data.frame(res_bcm_compare_gam_status_01)
  ),
  file = comparison_xlsx_file_01
)

build_bcm_compare_perf_plot_01 <- function(model_entry, activity_label, model_family = c("LMER", "GAM")) {
  model_family <- match.arg(model_family)

  if (
    is.null(model_entry) ||
      is.null(model_entry$fits) ||
      !requireNamespace("performance", quietly = TRUE) ||
      !requireNamespace("see", quietly = TRUE)
  ) {
    return(NULL)
  }

  models_named <- model_entry$fits
  models_named <- models_named[!map_lgl(models_named, is.null)]

  if (length(models_named) < 2) {
    return(NULL)
  }

  cmp_obj <- tryCatch(
    do.call(performance::compare_performance, models_named) |>
      relabel_compare_performance_tis_01(model_labels = names(models_named)),
    error = function(e) NULL
  )

  if (is.null(cmp_obj)) {
    return(NULL)
  }

  p_obj <- tryCatch(plot(cmp_obj), error = function(e) NULL)

  if (is.null(p_obj) || !inherits(p_obj, "ggplot")) {
    return(NULL)
  }

  p_obj +
    ggplot2::labs(title = paste0(model_family, ": ", activity_label)) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 34, face = "bold"),
      axis.title = ggplot2::element_text(size = 31, face = "bold"),
      axis.text = ggplot2::element_text(size = 29, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 32, hjust = 1),
      legend.title = ggplot2::element_text(size = 30, face = "bold"),
      legend.text = ggplot2::element_text(size = 28, face = "bold")
    )
}

res_bcm_compare_lmer_plot_input_01 <- imap(
  activity_only_candidates_01,
  \(predictor_col, nm) {
    activity_label <- activity_label_map_01[[nm]]
    dat_i <- d18_tis_bcm_model_base_01 |>
      transmute(
        immet_id,
        exam,
        bcm,
        sex,
        age_m0,
        tis_total = as.numeric(tis_total),
        activity_value = as.numeric(.data[[predictor_col]])
      )

    fit_obj <- fit_bcm_compare_models_lmer_01(dat_i)
    model_labels <- make_tis_activity_model_labels_01(activity_label)
    fits_named <- list(
      fit_obj$model_a,
      fit_obj$model_b,
      fit_obj$model_c
    )
    names(fits_named) <- unname(model_labels[c("model_a", "model_b", "model_c")])

    list(activity_label = activity_label, fits = fits_named)
  }
)

res_bcm_compare_gam_plot_input_01 <- imap(
  activity_only_candidates_01,
  \(predictor_col, nm) {
    activity_label <- activity_label_map_01[[nm]]
    dat_i <- d18_tis_bcm_model_base_01 |>
      transmute(
        immet_id,
        exam,
        bcm,
        sex,
        age_m0,
        tis_total = as.numeric(tis_total),
        activity_value = as.numeric(.data[[predictor_col]])
      )

    fit_obj <- fit_bcm_compare_models_gam_01(dat_i)
    model_labels <- make_tis_activity_model_labels_01(activity_label)
    fits_named <- list(
      fit_obj$model_a,
      fit_obj$model_b,
      fit_obj$model_c
    )
    names(fits_named) <- unname(model_labels[c("model_a", "model_b", "model_c")])

    list(activity_label = activity_label, fits = fits_named)
  }
)

res_bcm_compare_lmer_plot_list_01 <- imap(
  res_bcm_compare_lmer_plot_input_01,
  \(entry, nm) build_bcm_compare_perf_plot_01(
    model_entry = entry,
    activity_label = entry$activity_label,
    model_family = "LMER"
  )
) |>
  purrr::compact()

res_bcm_compare_gam_plot_list_01 <- imap(
  res_bcm_compare_gam_plot_input_01,
  \(entry, nm) build_bcm_compare_perf_plot_01(
    model_entry = entry,
    activity_label = entry$activity_label,
    model_family = "GAM"
  )
) |>
  purrr::compact()

if (length(res_bcm_compare_lmer_plot_list_01) > 0) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    res_bcm_compare_lmer_multi_plot_01 <- patchwork::wrap_plots(
      res_bcm_compare_lmer_plot_list_01,
      ncol = 2
    ) +
      patchwork::plot_annotation(
        title = "BCM models: compare_performance (Model 1/2/3)",
        subtitle = "LMER: TIS vs activity parameters",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 18, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = 12, face = "bold")
        )
      ) &
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 31, face = "bold"),
        axis.text = ggplot2::element_text(size = 29, face = "bold"),
        strip.text = ggplot2::element_text(size = 31, face = "bold"),
        legend.title = ggplot2::element_text(size = 30, face = "bold"),
        legend.text = ggplot2::element_text(size = 28, face = "bold")
      )
  } else {
    res_bcm_compare_lmer_multi_plot_01 <- NULL
  }
} else {
  res_bcm_compare_lmer_multi_plot_01 <- NULL
}

if (length(res_bcm_compare_gam_plot_list_01) > 0) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    res_bcm_compare_gam_multi_plot_01 <- patchwork::wrap_plots(
      res_bcm_compare_gam_plot_list_01,
      ncol = 2
    ) +
      patchwork::plot_annotation(
        title = "BCM models: compare_performance (Model 1/2/3)",
        subtitle = "GAM: TIS vs activity parameters",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 18, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = 12, face = "bold")
        )
      ) &
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 31, face = "bold"),
        axis.text = ggplot2::element_text(size = 29, face = "bold"),
        strip.text = ggplot2::element_text(size = 31, face = "bold"),
        legend.title = ggplot2::element_text(size = 30, face = "bold"),
        legend.text = ggplot2::element_text(size = 28, face = "bold")
      )
  } else {
    res_bcm_compare_gam_multi_plot_01 <- NULL
  }
} else {
  res_bcm_compare_gam_multi_plot_01 <- NULL
}

save_compare_perf_plot_01 <- function(plot_obj, file_nm) {
  if (is.null(plot_obj)) {
    return(FALSE)
  }

  if (file.exists(file_nm)) {
    removed_old_01 <- isTRUE(file.remove(file_nm))
    if (!removed_old_01) {
      warning(
        paste0(
          "Existing file could not be removed before overwrite: ",
          file_nm
        )
      )
    }
  }

  tryCatch(
    {
      ggsave(
        filename = file_nm,
        plot = plot_obj,
        width = 20,
        height = 12,
        units = "in",
        dpi = 300,
        bg = "white",
        limitsize = FALSE
      )
      TRUE
    },
    error = function(e) FALSE
  )
}

res_bcm_compare_lmer_plot_file_01 <- file.path(
  "output",
  "figures",
  "tis",
  paste0(format(Sys.Date(), "%y%m%d"), "_bcm_tis_activity_compare_performance_lmer_01.png")
)

res_bcm_compare_gam_plot_file_01 <- file.path(
  "output",
  "figures",
  "tis",
  paste0(format(Sys.Date(), "%y%m%d"), "_bcm_tis_activity_compare_performance_gam_01.png")
)

res_bcm_compare_lmer_plot_saved_01 <- save_compare_perf_plot_01(
  plot_obj = res_bcm_compare_lmer_multi_plot_01,
  file_nm = res_bcm_compare_lmer_plot_file_01
)

res_bcm_compare_gam_plot_saved_01 <- save_compare_perf_plot_01(
  plot_obj = res_bcm_compare_gam_multi_plot_01,
  file_nm = res_bcm_compare_gam_plot_file_01
)


