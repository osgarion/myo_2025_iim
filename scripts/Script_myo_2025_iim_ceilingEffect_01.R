# **********************************************************************
# 1. preparation ----
# **********************************************************************
## libraries ----
# **********************************************************************
# Functions and objects for ceiling/floor analysis only
source("scripts/functions/FUN_ceil_flo_01.R")
source("scripts/functions/OBJ_ceil_flo_01.R")
if (isTRUE(getOption("myo_run_renv_status", FALSE))) {
  renv::status()
}

# Default analysis configuration.
# NOTE: wrappers/background-job scripts can override these objects before sourcing this file.
if (!exists("sel_dep_var_ceil_01", inherits = FALSE)) {
  sel_dep_var_ceil_01 <- c("mmt8", "mmt10", "fi2", "fi3", "haq")
}
if (!exists("sel_indep_var_ceil_01", inherits = FALSE)) {
  sel_indep_var_ceil_01 <- c("ig_gradient_enmo")
}
if (!exists("output_root_ceil_01", inherits = FALSE)) {
  output_root_ceil_01 <- file.path("output", "ceil_flo")
}
if (!exists("run_repro_backup_block_ceil_01", inherits = FALSE)) {
  run_repro_backup_block_ceil_01 <- FALSE
}
output_fig_dir_ceil_01 <- file.path(output_root_ceil_01, "figures")
output_tab_dir_ceil_01 <- file.path(output_root_ceil_01, "tables")

walk(
  c(output_root_ceil_01, output_fig_dir_ceil_01, output_tab_dir_ceil_01),
  ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE)
)

# **********************************************************************
## Data ----
# **********************************************************************
d18_ggir_data_long <- d18_ggir_data |>
  select(any_of(sel_dep_var_ceil_01), any_of(sel_indep_var_ceil_01), exam) |>
  pivot_longer(
    cols = -exam,
    names_to = "var_name",
    values_to = "var_value"
  ) |>
  mutate(
    cei_flo = case_when(
      var_name == "mmt8" & var_value == 800 ~ 1,
      var_name == "mmt10" & var_value == 100 ~ 1,
      var_name == "fi2" & var_value == 100 ~ 1,
      var_name == "fi3" & var_value == 100 ~ 1,
      var_name == "haq" & var_value == 0 ~ 1,
      TRUE ~ 0
    )
  )

# Canonical exam ordering used across all multi-figure facets.
exam_levels_ref_01 <- c("M0", "M3", "M6", "M18", "M30", "M48")

# Normalize exam labels to canonical "M<number>" factor ordering.
normalize_exam_factor_01 <- function(x, levels_ref = exam_levels_ref_01) {
  x_chr <- stringr::str_to_upper(stringr::str_trim(as.character(x)))
  x_num <- suppressWarnings(readr::parse_number(x_chr))
  x_lab <- if_else(is.na(x_num), x_chr, paste0("M", as.integer(x_num)))
  factor(x_lab, levels = levels_ref)
}

# **********************************************************************
# 2. Environment reproducibility ----
# **********************************************************************
# RENV
# renv::init()        # once
# renv::status()      # sync?
# renv::snapshot()    # after changing deps
# renv::restore()     # on a new machine

if (isTRUE(run_repro_backup_block_ceil_01)) {
  # Back-up
  back_up("scripts/functions/FUN_01.R") # the the destination subdirectory specify using 'path_dest'
  back_up("scripts/functions/OBJ_01.R") # the the destination subdirectory specify using 'path_dest'
  back_up("scripts/Script_myo_2025_iim_02_ggir.R") # the the destination subdirectory specify using 'path_dest'

  # Save to .RData
  save(xx1,
    file = "reports/markD_01.RData"
  )
  save.image("reports/markD_01.RData")
  save.image("data/260123_backup.RData")
  # Save the tables into data/tables.RData by listing them individually
  cgwtools::resave(xx2, xx3, file = "reports/markD_01.RData") # resave a list of tables that I'll use in the .Rmd file.
  # Save the tables into data/tables.RData using "patterns"
  cgwtools::resave(list = ls(pattern = "tbl"), file = "reports/markD_01.RData")
}

# *******************************************************
# 3. EDA ----
# *******************************************************
tiff(
  filename = file.path(
    output_fig_dir_ceil_01,
    paste0(
      format(Sys.Date(), "%y%m%d"),
      "_",
      "ggpairs_01.tiff"
    )
  ),
  compression = "lzw",
  res = "150",
  width = 1600,
  height = 1200
)
d18_ggir_data |>
  select(
    contains("enmo"),
    contains("mvpa"),
    starts_with("ig"), starts_with("mmt"),
    exam, haq, fi2, fi3,
    -matches("_b\\d+")
  ) |>
  ggpairs(aes(color = exam, alpha = 0.5),
    lower = list(continuous = "smooth")
  )
dev.off()

library(dplyr)
library(purrr)
library(corrplot)

d18_examlevel_ggir_clinic_corr_input <- d18_ggir_data |>
  select(
    contains("enmo"),
    contains("mvpa"),
    starts_with("ig"),
    starts_with("mmt"),
    haq, fi2, fi3,
    exam,
    -matches("_b\\d+")
  )

split(
  d18_examlevel_ggir_clinic_corr_input,
  d18_examlevel_ggir_clinic_corr_input$exam
) |>
  walk(function(df_exam) {
    df_num <- df_exam |>
      select(-exam) |>
      select(where(is.numeric)) |>
      select(where(~ {
        s <- suppressWarnings(stats::sd(.x, na.rm = TRUE))
        is.finite(s) && (s > 0)
      }))

    if (ncol(df_num) > 1) {
      M <- suppressWarnings(
        cor(
          df_num,
          method = "spearman",
          use = "pairwise.complete.obs"
        )
      )

      # Pairwise-complete correlations can produce NA/Inf blocks.
      # Keep only variables with fully finite correlations to avoid hclust failure.
      if (is.matrix(M) && ncol(M) > 1) {
        diag(M) <- 1
        keep_cols <- apply(M, 2, function(v) all(is.finite(v)))
        M <- M[keep_cols, keep_cols, drop = FALSE]
      }

      if (is.matrix(M) && ncol(M) > 1) {
        corrplot(
          M,
          type = "upper",
          order = "hclust",
          tl.col = "black",
          tl.srt = 45,
          main = paste("Exam:", unique(df_exam$exam))
        )
      } else {
        message(
          "Skipping corrplot for exam ",
          unique(df_exam$exam),
          " (insufficient finite correlation matrix)."
        )
      }
    }
  })


# **********************************************************************
## Histogram & ECDF ----
# **********************************************************************
# histograms
d18_ggir_data_long |>
  ggplot(aes(x = var_value)) +
  geom_density(aes(fill = exam),
    bins = 30, alpha = 0.5,
    position = position_dodge()
  ) +
  scale_fill_viridis_d() +
  facet_wrap(~var_name, scales = "free") +
  labs(
    title = "Histogram by exam",
    subtitle = "All cohort"
  )
sjPlot::theme_sjplot2() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
    axis.title = ggplot2::element_text(size = 12, face = "bold"),
    axis.text = ggplot2::element_text(size = 11, face = "bold")
  )

d18_ggir_data_long |>
  filter(cei_flo == 0) |>
  ggplot(aes(x = var_value)) +
  geom_density(aes(fill = exam),
    bins = 30, alpha = 0.5,
    position = position_dodge()
  ) +
  scale_fill_viridis_d() +
  facet_wrap(~var_name, scales = "free") +
  labs(
    title = "Histogram by exam",
    subtitle = "Without ceiling/floor cohort"
  )
sjPlot::theme_sjplot2() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
    axis.title = ggplot2::element_text(size = 12, face = "bold"),
    axis.text = ggplot2::element_text(size = 11, face = "bold")
  )

# eCDF
d18_ggir_data_long |>
  ggplot(aes(x = var_value)) +
  stat_ecdf(aes(color = exam, linetype = exam),
    geom = "point", size = 1.5
  ) +
  scale_fill_viridis_d() +
  facet_wrap(~var_name, scales = "free") +
  labs(
    title = "eCDF by exam",
    subtitle = "All cohort"
  )
sjPlot::theme_sjplot2() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
    axis.title = ggplot2::element_text(size = 12, face = "bold"),
    axis.text = ggplot2::element_text(size = 11, face = "bold")
  )


d18_ggir_data_long |>
  filter(cei_flo == 0) |>
  ggplot(aes(x = var_value)) +
  stat_ecdf(aes(color = exam, linetype = exam),
    geom = "point", size = 1.5
  ) +
  scale_fill_viridis_d() +
  facet_wrap(~var_name, scales = "free") +
  labs(
    title = "eCDF by exam",
    subtitle = "Without ceiling/floor cohort"
  )
sjPlot::theme_sjplot2() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
    axis.title = ggplot2::element_text(size = 12, face = "bold"),
    axis.text = ggplot2::element_text(size = 11, face = "bold")
  )

# **********************************************************************
## Kendall correlation EDA: ig_gradient vs clinical outcomes ----
# **********************************************************************
x_col_eda_ceil_01 <- sel_indep_var_ceil_01[[1]]
x_lab_eda_ceil_01 <- stringr::str_remove(x_col_eda_ceil_01, "_enmo$")

# Safe Kendall test returning n, tau and p-value for one x/y vector pair.
kendall_test_safe_eda_01 <- function(x, y) {
  dat <- tibble(x = x, y = y) |>
    filter(!is.na(x), !is.na(y))

  if (nrow(dat) < 3 || dplyr::n_distinct(dat$x) < 2 || dplyr::n_distinct(dat$y) < 2) {
    return(tibble(n = nrow(dat), tau = NA_real_, p_value = NA_real_))
  }

  tst <- suppressWarnings(
    cor.test(dat$x, dat$y, method = "kendall", exact = FALSE)
  )

  tibble(
    n = nrow(dat),
    tau = unname(tst$estimate),
    p_value = tst$p.value
  )
}

# EDA long table for ig_gradient vs each outcome with boundary flag.
d18_eda_corr_long_ceil_01 <- d18_ggir_data |>
  select(any_of(c("exam", x_col_eda_ceil_01, sel_dep_var_ceil_01))) |>
  transmute(
    exam = normalize_exam_factor_01(exam),
    ig_gradient = as.numeric(.data[[x_col_eda_ceil_01]]),
    across(any_of(sel_dep_var_ceil_01), as.numeric)
  ) |>
  pivot_longer(
    cols = any_of(sel_dep_var_ceil_01),
    names_to = "outcome",
    values_to = "y"
  ) |>
  mutate(
    boundary_flag = case_when(
      outcome == "mmt8" ~ as.integer(y == 80),
      outcome %in% c("mmt10", "fi2", "fi3") ~ as.integer(y == 100),
      outcome == "haq" ~ as.integer(y == 0),
      TRUE ~ 0L
    )
  ) |>
  filter(!is.na(exam), !is.na(ig_gradient), !is.na(y))

exam_levels_eda_ceil_01 <- exam_levels_ref_01

# Kendall summary for all values.
res_eda_kendall_all_ceil_01 <- d18_eda_corr_long_ceil_01 |>
  group_by(outcome, exam) |>
  reframe(kendall_test_safe_eda_01(ig_gradient, y)) |>
  ungroup() |>
  rename(
    n_all = n,
    tau_all = tau,
    p_all = p_value
  )

# Kendall summary for no-ceiling/floor values.
res_eda_kendall_nb_ceil_01 <- d18_eda_corr_long_ceil_01 |>
  filter(boundary_flag == 0L) |>
  group_by(outcome, exam) |>
  reframe(kendall_test_safe_eda_01(ig_gradient, y)) |>
  ungroup() |>
  rename(
    n_no_cei_flo = n,
    tau_no_cei_flo = tau,
    p_no_cei_flo = p_value
  )

# Panel-level labels containing both all-cohort and no-ceiling/floor Kendall stats.
res_eda_kendall_tab_ceil_01 <- full_join(
  res_eda_kendall_all_ceil_01,
  res_eda_kendall_nb_ceil_01,
  by = c("outcome", "exam")
  ) |>
  mutate(
    exam = factor(exam, levels = exam_levels_eda_ceil_01),
    label = paste0(
      "ALL: n=", n_all,
      " | tau=", if_else(is.na(tau_all), "NA", as.character(round(tau_all, 3))),
      " | p=", if_else(is.na(p_all), "NA", as.character(scales::pvalue(p_all, accuracy = 0.001))),
      "\nNO_CEI/FLO: n=", n_no_cei_flo,
      " | tau=", if_else(is.na(tau_no_cei_flo), "NA", as.character(round(tau_no_cei_flo, 3))),
      " | p=", if_else(is.na(p_no_cei_flo), "NA", as.character(scales::pvalue(p_no_cei_flo, accuracy = 0.001)))
    )
  )

# Plot one multi-figure per outcome: facets by exam with both cohorts in each panel.
plot_kendall_exam_multi_ceil_01 <- function(outcome_i) {
  dat_all <- d18_eda_corr_long_ceil_01 |>
    filter(outcome == outcome_i)

  dat_no_cei_flo <- d18_eda_corr_long_ceil_01 |>
    filter(outcome == outcome_i, boundary_flag == 0L)

  dat_plot <- dat_all |>
    filter(outcome == outcome_i)

  lab_plot <- res_eda_kendall_tab_ceil_01 |>
    filter(outcome == outcome_i)

  ggplot(dat_plot, aes(x = ig_gradient, y = y)) +
    geom_point(alpha = 0.45, size = 3.0, color = "grey50") +
    geom_smooth(
      aes(color = "Trend: all values", linetype = "Trend: all values"),
      method = "lm",
      se = FALSE,
      linewidth = 1.2
    ) +
    geom_point(
      data = dat_no_cei_flo,
      aes(x = ig_gradient, y = y),
      inherit.aes = FALSE,
      alpha = 0.8,
      size = 3.0,
      color = "#2b8cbe"
    ) +
    geom_smooth(
      data = dat_no_cei_flo,
      aes(
        x = ig_gradient,
        y = y,
        color = "Trend: without ceiling/floor",
        linetype = "Trend: without ceiling/floor"
      ),
      inherit.aes = FALSE,
      method = "lm",
      se = FALSE,
      linewidth = 1.2
    ) +
    geom_text(
      data = lab_plot,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.02,
      vjust = 1.15,
      size = 3.6,
      lineheight = 1.0,
      fontface = "bold"
    ) +
    facet_wrap(~exam, scales = "free_y", ncol = 2) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0("Kendall correlation: ", x_lab_eda_ceil_01, " vs ", outcome_i),
      subtitle = "Grey = all values | Blue/Orange = without ceiling/floor",
      x = x_lab_eda_ceil_01,
      y = outcome_i
    ) +
    scale_color_manual(
      values = c(
        "Trend: all values" = "grey30",
        "Trend: without ceiling/floor" = "#d95f0e"
      ),
      name = "Trend"
    ) +
    scale_linetype_manual(
      values = c(
        "Trend: all values" = "solid",
        "Trend: without ceiling/floor" = "solid"
      ),
      name = "Trend"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      strip.text = element_text(face = "bold", size = 16),
      plot.title = element_text(face = "bold", size = 22),
      plot.subtitle = element_text(face = "bold", size = 15),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 13, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11, face = "bold")
    )
}

# Multi-figure list for all selected outcomes.
res_eda_kendall_plot_list_ceil_01 <- sel_dep_var_ceil_01 |>
  set_names() |>
  map(plot_kendall_exam_multi_ceil_01)

# Save one image per ig_gradient-outcome combination (safe export).
dir.create(output_fig_dir_ceil_01, recursive = TRUE, showWarnings = FALSE)
run_date_eda_ceil_01 <- format(Sys.Date(), "%y%m%d")

save_eda_plot_safe_01 <- function(plot_obj, outcome_nm) {
  fp <- file.path(
    output_fig_dir_ceil_01,
    paste0(run_date_eda_ceil_01, "_eda_kendall_", outcome_nm, "_01.png")
  )

  # Remove previous file to avoid silently keeping outdated image.
  removed_old <- TRUE
  if (file.exists(fp)) {
    removed_old <- isTRUE(file.remove(fp))
  }

  err_msg <- NA_character_
  ok <- tryCatch(
    {
      ggplot2::ggsave(
        filename = fp,
        plot = plot_obj,
        width = 14,
        height = 10,
        units = "in",
        dpi = 300,
        bg = "white",
        limitsize = FALSE
      )
      TRUE
    },
    error = function(e) {
      err_msg <<- conditionMessage(e)
      FALSE
    }
  )

  tibble(
    outcome = outcome_nm,
    file_path = fp,
    removed_old = removed_old,
    export_ok = ok,
    file_exists_after = file.exists(fp),
    file_mtime = if_else(file.exists(fp), as.character(file.info(fp)$mtime), NA_character_),
    error_msg = err_msg
  )
}

res_eda_kendall_export_ceil_01 <- imap_dfr(
  res_eda_kendall_plot_list_ceil_01,
  ~ save_eda_plot_safe_01(.x, .y)
)

print(res_eda_kendall_export_ceil_01)

if (any(!res_eda_kendall_export_ceil_01$export_ok)) {
  stop(
    "EDA Kendall export failed for: ",
    paste(res_eda_kendall_export_ceil_01$outcome[!res_eda_kendall_export_ceil_01$export_ok], collapse = ", ")
  )
}

# **********************************************************************
## PhVAS dual-axis EDA: clinical parameter + ig_gradient by exam ----
# **********************************************************************
# Resolve PhVAS column in the current dataset.
phvas_col_eda_ceil_01 <- intersect(c("ph_vas", "PhVAS", "phvas"), names(d18_ggir_data))
phvas_col_eda_ceil_01 <- if (length(phvas_col_eda_ceil_01) > 0) phvas_col_eda_ceil_01[[1]] else NA_character_

if (!is.na(phvas_col_eda_ceil_01)) {
  # Long table for dual-axis plots: x = PhVAS, y1 = clinical outcome, y2 = ig_gradient.
  d18_eda_phvas_dual_long_01 <- d18_ggir_data |>
    select(any_of(c("exam", phvas_col_eda_ceil_01, x_col_eda_ceil_01, sel_dep_var_ceil_01))) |>
    transmute(
      exam = normalize_exam_factor_01(exam),
      exam_order_factor = normalize_exam_factor_01(exam),
      ph_vas = as.numeric(.data[[phvas_col_eda_ceil_01]]),
      ig_gradient = as.numeric(.data[[x_col_eda_ceil_01]]),
      across(any_of(sel_dep_var_ceil_01), as.numeric)
    ) |>
    pivot_longer(
      cols = any_of(sel_dep_var_ceil_01),
      names_to = "outcome",
      values_to = "y_clin"
    ) |>
    filter(
      !is.na(ph_vas),
      !is.na(ig_gradient),
      !is.na(y_clin),
      !is.na(exam_order_factor)
    )

  # One multi-figure per outcome: dual y-axis and facets by exam_order_factor.
  plot_eda_phvas_dual_multi_01 <- function(outcome_i) {
    dat_plot <- d18_eda_phvas_dual_long_01 |>
      filter(outcome == outcome_i)

    if (nrow(dat_plot) == 0) {
      return(NULL)
    }

    y_min <- suppressWarnings(min(dat_plot$y_clin, na.rm = TRUE))
    y_max <- suppressWarnings(max(dat_plot$y_clin, na.rm = TRUE))
    g_min <- suppressWarnings(min(dat_plot$ig_gradient, na.rm = TRUE))
    g_max <- suppressWarnings(max(dat_plot$ig_gradient, na.rm = TRUE))

    y_rng <- ifelse(is.finite(y_max - y_min) && (y_max - y_min) > 0, y_max - y_min, 1)
    g_rng <- ifelse(is.finite(g_max - g_min) && (g_max - g_min) > 0, g_max - g_min, 1)

    dat_plot <- dat_plot |>
      mutate(
        ig_scaled = ((ig_gradient - g_min) / g_rng) * y_rng + y_min
      )

    pred_lab <- x_lab_eda_ceil_01

    ggplot(dat_plot, aes(x = ph_vas)) +
      geom_point(aes(y = y_clin, color = "Clinical"), alpha = 0.7, size = 2.4) +
      geom_smooth(aes(y = y_clin, color = "Clinical"), method = "lm", se = FALSE, linewidth = 1.1) +
      geom_point(aes(y = ig_scaled, color = .env$pred_lab), alpha = 0.7, size = 2.4) +
      geom_smooth(aes(y = ig_scaled, color = .env$pred_lab), method = "lm", se = FALSE, linewidth = 1.1) +
      facet_wrap(~exam_order_factor, ncol = 2) +
      scale_color_manual(
        values = c("Clinical" = "#2b8cbe", setNames("#d95f0e", pred_lab)),
        name = NULL
      ) +
      scale_y_continuous(
        name = outcome_i,
        sec.axis = sec_axis(
          ~ ((. - y_min) / y_rng) * g_rng + g_min,
          name = pred_lab
        )
      ) +
      labs(
        title = paste0("PhVAS vs ", outcome_i, " + ", pred_lab),
        subtitle = "Facet: exam_order_factor",
        x = "ph_vas"
      ) +
      theme_bw(base_size = 16) +
      theme(
        strip.text = element_text(face = "bold", size = 14),
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold")
      )
  }

  res_eda_phvas_dual_plot_list_01 <- sel_dep_var_ceil_01 |>
    set_names() |>
    map(plot_eda_phvas_dual_multi_01) |>
    purrr::compact()

  save_eda_phvas_dual_plot_safe_01 <- function(plot_obj, outcome_nm) {
    fp <- file.path(
      output_fig_dir_ceil_01,
      paste0(run_date_eda_ceil_01, "_eda_phvas_dual_", outcome_nm, "_01.png")
    )

    removed_old <- TRUE
    if (file.exists(fp)) {
      removed_old <- isTRUE(file.remove(fp))
    }

    err_msg <- NA_character_
    ok <- tryCatch(
      {
        ggplot2::ggsave(
          filename = fp,
          plot = plot_obj,
          width = 14,
          height = 10,
          units = "in",
          dpi = 300,
          bg = "white",
          limitsize = FALSE
        )
        TRUE
      },
      error = function(e) {
        err_msg <<- conditionMessage(e)
        FALSE
      }
    )

    tibble(
      outcome = outcome_nm,
      file_path = fp,
      removed_old = removed_old,
      export_ok = ok,
      file_exists_after = file.exists(fp),
      file_mtime = if_else(file.exists(fp), as.character(file.info(fp)$mtime), NA_character_),
      error_msg = err_msg
    )
  }

  res_eda_phvas_dual_export_01 <- imap_dfr(
    res_eda_phvas_dual_plot_list_01,
    ~ save_eda_phvas_dual_plot_safe_01(.x, .y)
  )

  print(res_eda_phvas_dual_export_01)

  if (any(!res_eda_phvas_dual_export_01$export_ok)) {
    stop(
      "EDA PhVAS dual-axis export failed for: ",
      paste(res_eda_phvas_dual_export_01$outcome[!res_eda_phvas_dual_export_01$export_ok], collapse = ", ")
    )
  }

  # Kendall tau correlations by outcome and exam:
  # 1) PhVAS vs ig_gradient
  # 2) PhVAS vs clinical outcome
  calc_kendall_safe_01 <- function(x, y) {
    dat <- tibble(x = as.numeric(x), y = as.numeric(y)) |>
      filter(is.finite(x), is.finite(y))

    if (nrow(dat) < 3 || dplyr::n_distinct(dat$x) < 2 || dplyr::n_distinct(dat$y) < 2) {
      return(tibble(n = nrow(dat), tau = NA_real_, p_value = NA_real_))
    }

    tst <- suppressWarnings(cor.test(dat$x, dat$y, method = "kendall", exact = FALSE))
    tibble(
      n = nrow(dat),
      tau = unname(tst$estimate),
      p_value = tst$p.value
    )
  }

  # Pair: PhVAS vs ig_gradient (same for all outcomes but exported per outcome-exam to match figures).
  res_eda_phvas_dual_kendall_ig_01 <- d18_eda_phvas_dual_long_01 |>
    group_by(outcome, exam_order_factor) |>
    group_modify(~ calc_kendall_safe_01(.x$ph_vas, .x$ig_gradient)) |>
    ungroup() |>
    mutate(
      pair = paste0("phvas_vs_", x_lab_eda_ceil_01),
      x_var = "ph_vas",
      y_var = x_lab_eda_ceil_01,
      .before = 1
    )

  # Pair: PhVAS vs clinical outcome.
  res_eda_phvas_dual_kendall_clin_01 <- d18_eda_phvas_dual_long_01 |>
    group_by(outcome, exam_order_factor) |>
    group_modify(~ calc_kendall_safe_01(.x$ph_vas, .x$y_clin)) |>
    ungroup() |>
    mutate(
      pair = "phvas_vs_clinical",
      x_var = "ph_vas",
      y_var = "y_clin",
      .before = 1
    )

  # Combined long table for export.
  res_eda_phvas_dual_kendall_01 <- bind_rows(
    res_eda_phvas_dual_kendall_ig_01,
    res_eda_phvas_dual_kendall_clin_01
  ) |>
    mutate(
      exam_order_factor = as.character(exam_order_factor),
      p_label = if_else(is.na(p_value), "NA", as.character(scales::pvalue(p_value, accuracy = 0.001)))
    )

  # Optional wide summary for quick reading per exam.
  res_eda_phvas_dual_kendall_wide_01 <- res_eda_phvas_dual_kendall_01 |>
    select(outcome, exam_order_factor, pair, n, tau, p_value) |>
    pivot_wider(
      names_from = pair,
      values_from = c(n, tau, p_value),
      names_sep = "__"
    )

  # Export both plot manifest and Kendall tables.
  rio::export(
    list(
      phvas_dual_plot_manifest = janitor::clean_names(res_eda_phvas_dual_export_01),
      phvas_dual_kendall_long = janitor::clean_names(res_eda_phvas_dual_kendall_01),
      phvas_dual_kendall_wide = janitor::clean_names(res_eda_phvas_dual_kendall_wide_01)
    ),
    file.path(
      output_tab_dir_ceil_01,
      paste0(run_date_eda_ceil_01, "_eda_phvas_dual_kendall_01.xlsx")
    )
  )
} else {
  warning("PhVAS column not found for dual-axis EDA plots (expected one of: ph_vas, PhVAS, phvas).")
}

# **********************************************************************
## Exam-trend dual-axis EDA: ig_gradient + clinical outcome (<=M30) ----
# **********************************************************************
# Resolve ID and exam columns for the trend analysis block.
id_col_eda_trend_01 <- intersect(c("patient_id", "immet_id", "iimet_id", "projekt_id", "bbm_id", "sample"), names(d18_ggir_data))
id_col_eda_trend_01 <- if (length(id_col_eda_trend_01) > 0) id_col_eda_trend_01[[1]] else NA_character_

exam_col_eda_trend_01 <- intersect(c("exam_order", "exam", "poradie_vysetrenia"), names(d18_ggir_data))
exam_col_eda_trend_01 <- if (length(exam_col_eda_trend_01) > 0) exam_col_eda_trend_01[[1]] else NA_character_

if (is.na(id_col_eda_trend_01) || is.na(exam_col_eda_trend_01)) {
  warning("EDA trend dual-axis block skipped: missing ID or exam column.")
} else {
  # Long EDA table: one row per patient-visit-outcome with exams limited to <= M30.
  d18_eda_exam_dual_long_01 <- d18_ggir_data |>
    select(any_of(c(id_col_eda_trend_01, exam_col_eda_trend_01, x_col_eda_ceil_01, sel_dep_var_ceil_01))) |>
    transmute(
      patient_id = as.character(.data[[id_col_eda_trend_01]]),
      exam_raw = as.character(.data[[exam_col_eda_trend_01]]),
      exam_month = readr::parse_number(exam_raw),
      exam_label = paste0("M", exam_month),
      ig_gradient = as.numeric(.data[[x_col_eda_ceil_01]]),
      across(any_of(sel_dep_var_ceil_01), as.numeric)
    ) |>
    filter(!is.na(patient_id), !is.na(exam_month), exam_month <= 30) |>
    pivot_longer(
      cols = any_of(sel_dep_var_ceil_01),
      names_to = "outcome",
      values_to = "y"
    ) |>
    filter(!is.na(ig_gradient), !is.na(y)) |>
    mutate(
      exam_label = factor(exam_label, levels = exam_levels_ref_01)
    )

  # Outcome-specific bounds for fixed Y1 axis (ceiling + 5; HAQ floor anchored at 0).
  mmt8_max_obs_eda_01 <- d18_eda_exam_dual_long_01 |>
    filter(outcome == "mmt8") |>
    summarise(v = suppressWarnings(max(y, na.rm = TRUE)), .groups = "drop") |>
    pull(v)
  mmt8_ceiling_eda_01 <- if (length(mmt8_max_obs_eda_01) == 1 && is.finite(mmt8_max_obs_eda_01) && mmt8_max_obs_eda_01 > 100) 800 else 80

  res_eda_exam_dual_bounds_01 <- tibble(
    outcome = sel_dep_var_ceil_01,
    boundary_type = if_else(outcome == "haq", "floor", "ceiling"),
    boundary_value = case_when(
      outcome == "mmt8" ~ as.numeric(mmt8_ceiling_eda_01),
      outcome %in% c("mmt10", "fi2", "fi3") ~ 100,
      outcome == "haq" ~ 0,
      TRUE ~ NA_real_
    ),
    y1_lower = if_else(outcome == "haq", 0, 0),
    y1_upper = case_when(
      outcome == "mmt8" ~ as.numeric(mmt8_ceiling_eda_01 + 5),
      outcome %in% c("mmt10", "fi2", "fi3") ~ 105,
      outcome == "haq" ~ 1.5,
      TRUE ~ NA_real_
    )
  )

  d18_eda_exam_dual_long_01 <- d18_eda_exam_dual_long_01 |>
    left_join(select(res_eda_exam_dual_bounds_01, outcome, boundary_value), by = "outcome") |>
    mutate(
      boundary_flag = as.integer(dplyr::near(y, boundary_value))
    )

  # Build one dual-axis trend plot per outcome for boundary-patient cohorts.
  plot_eda_exam_dual_m30_01 <- function(outcome_i) {
    dat_i <- d18_eda_exam_dual_long_01 |>
      filter(outcome == outcome_i)

    boundary_patients <- dat_i |>
      group_by(patient_id) |>
      summarise(has_boundary = any(boundary_flag == 1L, na.rm = TRUE), .groups = "drop") |>
      filter(has_boundary) |>
      pull(patient_id)

    dat_i <- dat_i |>
      filter(patient_id %in% boundary_patients)

    if (nrow(dat_i) == 0) {
      return(NULL)
    }

    dat_med <- dat_i |>
      group_by(exam_month, exam_label) |>
      summarise(
        y_median = median(y, na.rm = TRUE),
        ig_median = median(ig_gradient, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) |>
      arrange(exam_month)

    if (nrow(dat_med) == 0) {
      return(NULL)
    }

    bounds_i <- res_eda_exam_dual_bounds_01 |>
      filter(outcome == outcome_i) |>
      slice_head(n = 1)

    y1_high <- bounds_i$y1_upper[[1]]
    if (!is.finite(y1_high)) {
      y1_high <- suppressWarnings(max(dat_med$y_median, na.rm = TRUE))
    }

    range_y <- range(dat_med$y_median, na.rm = TRUE)
    range_ig <- range(dat_med$ig_median, na.rm = TRUE)
    if (!all(is.finite(c(range_y, range_ig)))) {
      return(NULL)
    }

    if ((range_ig[2] - range_ig[1]) <= 0 || (range_y[2] - range_y[1]) <= 0) {
      return(NULL)
    }

    # Use the same range-to-range scaling method as in the user prototype.
    a <- (range_y[2] - range_y[1]) / (range_ig[2] - range_ig[1])
    b <- range_y[1] - a * range_ig[1]

    dat_med <- dat_med |>
      mutate(ig_scaled = a * ig_median + b)

    y1_low_plot <- if (outcome_i == "haq") 0 else NA_real_

    p_main <- ggplot(dat_med, aes(x = exam_month)) +
      geom_point(aes(y = y_median, color = "Clinical"), size = 3) +
      geom_line(aes(y = y_median, color = "Clinical"), linewidth = 2.4) +
      geom_point(aes(y = ig_scaled, color = x_lab_eda_ceil_01), shape = 1, size = 3, stroke = 1.2) +
      geom_line(aes(y = ig_scaled, color = x_lab_eda_ceil_01, linetype = x_lab_eda_ceil_01), linewidth = 2.4) +
      scale_x_continuous(
        breaks = dat_med$exam_month,
        labels = paste0("M", dat_med$exam_month)
      ) +
      scale_y_continuous(
        name = paste0("Median ", outcome_i),
        sec.axis = sec_axis(
          ~ (. - b) / a,
          name = paste0("Median ", x_lab_eda_ceil_01)
        )
      ) +
      scale_color_manual(
        name = "Measure",
        values = c("Clinical" = "black", setNames("darkred", x_lab_eda_ceil_01))
      ) +
      scale_linetype_manual(
        name = "Measure",
        values = setNames("dashed", x_lab_eda_ceil_01)
      ) +
      labs(
        title = paste0("Boundary-patient trend (<=M30): ", outcome_i, " + ", x_lab_eda_ceil_01),
        subtitle = "Only patients with at least one ceiling/floor visit in the selected outcome",
        x = "Exam"
      ) +
      coord_cartesian(ylim = c(y1_low_plot, y1_high)) +
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
      distinct(patient_id, exam_month, exam_label) |>
      count(exam_month, exam_label, name = "n_patients") |>
      arrange(exam_month)

    exam_levels_tab <- tab_counts |>
      pull(exam_label) |>
      as.character()

    tab_plot_df <- bind_rows(
      tab_counts |>
        transmute(
          exam_label = factor(as.character(exam_label), levels = exam_levels_tab),
          row = "Exam",
          value = as.character(exam_label)
        ),
      tab_counts |>
        transmute(
          exam_label = factor(as.character(exam_label), levels = exam_levels_tab),
          row = "N patients",
          value = as.character(n_patients)
        )
    ) |>
      mutate(row = factor(row, levels = c("N patients", "Exam")))

    p_table <- ggplot(tab_plot_df, aes(x = exam_label, y = row)) +
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

  res_eda_exam_dual_plot_list_01 <- sel_dep_var_ceil_01 |>
    set_names() |>
    map(plot_eda_exam_dual_m30_01) |>
    purrr::compact()

  save_eda_exam_dual_plot_safe_01 <- function(plot_obj, outcome_nm) {
    fp <- file.path(
      output_fig_dir_ceil_01,
      paste0(run_date_eda_ceil_01, "_eda_exam_dual_m30_", outcome_nm, "_01.png")
    )

    removed_old <- TRUE
    if (file.exists(fp)) {
      removed_old <- isTRUE(file.remove(fp))
    }

    err_msg <- NA_character_
    ok <- tryCatch(
      {
        ggplot2::ggsave(
          filename = fp,
          plot = plot_obj,
          width = 24,
          height = 18,
          units = "in",
          dpi = 300,
          bg = "white",
          limitsize = FALSE
        )
        TRUE
      },
      error = function(e) {
        err_msg <<- conditionMessage(e)
        FALSE
      }
    )

    tibble(
      outcome = outcome_nm,
      file_path = fp,
      removed_old = removed_old,
      export_ok = ok,
      file_exists_after = file.exists(fp),
      file_mtime = if_else(file.exists(fp), as.character(file.info(fp)$mtime), NA_character_),
      error_msg = err_msg
    )
  }

  res_eda_exam_dual_export_01 <- imap_dfr(
    res_eda_exam_dual_plot_list_01,
    ~ save_eda_exam_dual_plot_safe_01(.x, .y)
  )

  print(res_eda_exam_dual_export_01)

  if (any(!res_eda_exam_dual_export_01$export_ok)) {
    stop(
      "EDA exam dual-axis export failed for: ",
      paste(res_eda_exam_dual_export_01$outcome[!res_eda_exam_dual_export_01$export_ok], collapse = ", ")
    )
  }
}

# **********************************************************************
# 4. ig_gradient -> clinical scale conversion (hurdle + OOF BA-like) ----
# **********************************************************************

# Select the first existing column name from a candidate list.
pick_first_col_01 <- function(data, candidates, allow_missing = FALSE) {
  hit <- intersect(candidates, names(data))
  if (length(hit) == 0) {
    if (allow_missing) {
      return(NA_character_)
    }
    stop("No required column found among: ", paste(candidates, collapse = ", "))
  }
  hit[[1]]
}

# Clamp numeric values into a fixed interval.
clamp_01 <- function(x, low, high) {
  pmin(pmax(x, low), high)
}

# Save ggplot object safely without stopping the pipeline.
save_plot_safe_01 <- function(plot_obj, file_path, width = 8, height = 5, dpi = 300) {
  tryCatch(
    {
      ggplot2::ggsave(filename = file_path, plot = plot_obj, width = width, height = height, dpi = dpi)
      TRUE
    },
    error = function(e) FALSE
  )
}

# Convert exam labels (e.g., M0, M6) into numeric month values.
parse_exam_time_01 <- function(exam_order_chr) {
  lab <- stringr::str_to_lower(as.character(exam_order_chr))
  lab <- stringr::str_trim(lab)
  out <- suppressWarnings(as.numeric(stringr::str_extract(lab, "\\d+")))
  out
}

# Keep only covariates that exist and have non-missing variation.
get_covars_use_01 <- function(train_df, covars) {
  covars_present <- intersect(covars, names(train_df))
  covars_present[map_lgl(covars_present, ~ {
    vv <- train_df[[.x]]
    any(!is.na(vv)) && dplyr::n_distinct(vv[!is.na(vv)]) > 1
  })]
}

# Build grouped CV folds so the same patient never appears in train and test.
make_group_folds_01 <- function(ids, k = 5, seed = 20260303) {
  ids_unique <- unique(ids[!is.na(ids)])
  n_ids <- length(ids_unique)
  if (n_ids == 0) {
    stop("No valid patient IDs for grouped CV.")
  }
  k_use <- min(k, n_ids)
  set.seed(seed)
  ids_shuffled <- sample(ids_unique, size = n_ids, replace = FALSE)
  tibble(
    patient_id = ids_shuffled,
    fold = rep(seq_len(k_use), length.out = n_ids)
  )
}

# Compute a correlation safely (returns NA for degenerate inputs).
safe_cor_01 <- function(x, y, method) {
  keep <- complete.cases(x, y)
  if (sum(keep) < 3 || sd(x[keep]) == 0 || sd(y[keep]) == 0) {
    return(NA_real_)
  }
  suppressWarnings(cor(x[keep], y[keep], method = method))
}

# Compute ROC AUC from binary labels and continuous scores.
calc_auc_01 <- function(y_true, score) {
  keep <- !is.na(y_true) & !is.na(score)
  y <- as.integer(y_true[keep])
  s <- as.numeric(score[keep])

  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  if (n_pos == 0 || n_neg == 0) {
    return(NA_real_)
  }

  ranks <- rank(s, ties.method = "average")
  (sum(ranks[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

# Compute PR-AUC (average precision) from binary labels and scores.
calc_pr_auc_01 <- function(y_true, score) {
  keep <- !is.na(y_true) & !is.na(score)
  y <- as.integer(y_true[keep])
  s <- as.numeric(score[keep])

  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  if (n_pos == 0 || n_neg == 0) {
    return(NA_real_)
  }

  ord <- order(s, decreasing = TRUE)
  y <- y[ord]

  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  precision <- tp / (tp + fp)
  idx_pos <- which(y == 1L)

  mean(precision[idx_pos])
}

# Compute global regression performance and linear calibration metrics.
calc_reg_metrics_01 <- function(df) {
  n <- nrow(df)
  mae <- mean(abs(df$y - df$y_hat), na.rm = TRUE)
  rmse <- sqrt(mean((df$y - df$y_hat)^2, na.rm = TRUE))

  sse <- sum((df$y - df$y_hat)^2, na.rm = TRUE)
  sst <- sum((df$y - mean(df$y, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- ifelse(sst == 0, NA_real_, 1 - sse / sst)

  if (n >= 3 && sd(df$y_hat, na.rm = TRUE) > 0) {
    fit_cal <- lm(y ~ y_hat, data = df)
    cal_intercept <- unname(coef(fit_cal)[1])
    cal_slope <- unname(coef(fit_cal)[2])
  } else {
    cal_intercept <- NA_real_
    cal_slope <- NA_real_
  }

  tibble(
    n = n,
    mae = mae,
    rmse = rmse,
    r2 = r2,
    spearman = safe_cor_01(df$y, df$y_hat, "spearman"),
    kendall = safe_cor_01(df$y, df$y_hat, "kendall"),
    calibration_intercept = cal_intercept,
    calibration_slope = cal_slope
  )
}

# Compute BA-like bias, LoA, and proportional-bias trend statistics.
calc_ba_stats_01 <- function(df) {
  if (nrow(df) < 2) {
    return(
      tibble(
        n = nrow(df),
        bias = NA_real_,
        sd_diff = NA_real_,
        loa_lower = NA_real_,
        loa_upper = NA_real_,
        trend_slope = NA_real_,
        trend_p = NA_real_
      )
    )
  }

  bias <- mean(df$ba_diff, na.rm = TRUE)
  sd_diff <- sd(df$ba_diff, na.rm = TRUE)
  loa_lower <- bias - 1.96 * sd_diff
  loa_upper <- bias + 1.96 * sd_diff

  if (nrow(df) >= 3 && sd(df$ba_mean, na.rm = TRUE) > 0) {
    fit_trend <- lm(ba_diff ~ ba_mean, data = df)
    trend_slope <- unname(coef(fit_trend)[2])
    trend_p <- summary(fit_trend)$coefficients[2, 4]
  } else {
    trend_slope <- NA_real_
    trend_p <- NA_real_
  }

  tibble(
    n = nrow(df),
    bias = bias,
    sd_diff = sd_diff,
    loa_lower = loa_lower,
    loa_upper = loa_upper,
    trend_slope = trend_slope,
    trend_p = trend_p
  )
}

# Compute MAE/RMSE for a selected prediction column.
calc_err_vs_pred_01 <- function(df, pred_col) {
  tibble(
    mae = mean(abs(df$y - df[[pred_col]]), na.rm = TRUE),
    rmse = sqrt(mean((df$y - df[[pred_col]])^2, na.rm = TRUE))
  )
}

# Compute correlation estimate and p-value for a chosen method.
calc_cor_test_01 <- function(df, predictor_col, method) {
  dat <- df |>
    select(all_of(c("y", predictor_col))) |>
    filter(!is.na(y), !is.na(.data[[predictor_col]]))

  if (nrow(dat) < 5 || dplyr::n_distinct(dat[[predictor_col]]) < 2 || dplyr::n_distinct(dat$y) < 2) {
    return(tibble(estimate = NA_real_, p_value = NA_real_))
  }

  tst <- tryCatch(
    suppressWarnings(cor.test(dat[[predictor_col]], dat$y, method = method, exact = FALSE)),
    error = function(e) NULL
  )

  if (is.null(tst)) {
    return(tibble(estimate = NA_real_, p_value = NA_real_))
  }

  tibble(estimate = unname(tst$estimate), p_value = tst$p.value)
}

# Fit non-ceiling association model with spline predictor and optional covariates.
fit_assoc_ns_model_01 <- function(df, predictor_col, covars_pref = c("age_m0_cov", "exam_order_factor")) {
  covars_use <- get_covars_use_01(df, covars_pref)
  rhs <- c(paste0("splines::ns(", predictor_col, ", df = 3)"), covars_use)
  form <- as.formula(paste0("y ~ ", paste(rhs, collapse = " + ")))

  fit <- tryCatch(lm(form, data = df), error = function(e) NULL)
  list(model = fit, covars_use = covars_use, formula = form)
}

# Compute marginal delta between P90 and P10 of a predictor from a fitted model.
calc_delta_p90_p10_01 <- function(fit_model, df, predictor_col) {
  if (is.null(fit_model)) {
    return(NA_real_)
  }

  p10 <- as.numeric(quantile(df[[predictor_col]], 0.10, na.rm = TRUE))
  p90 <- as.numeric(quantile(df[[predictor_col]], 0.90, na.rm = TRUE))
  if (!is.finite(p10) || !is.finite(p90)) {
    return(NA_real_)
  }

  dat_low <- df
  dat_high <- df
  dat_low[[predictor_col]] <- p10
  dat_high[[predictor_col]] <- p90

  pred_low <- tryCatch(predict(fit_model, newdata = dat_low), error = function(e) rep(NA_real_, nrow(dat_low)))
  pred_high <- tryCatch(predict(fit_model, newdata = dat_high), error = function(e) rep(NA_real_, nrow(dat_high)))
  mean(pred_high - pred_low, na.rm = TRUE)
}

# Bootstrap CI for P90-P10 effect size from the spline association model.
bootstrap_delta_p90_p10_01 <- function(df, predictor_col, covars_pref, n_boot = 200, seed = 20260303) {
  if (nrow(df) < 30) {
    return(tibble(ci_low = NA_real_, ci_high = NA_real_))
  }

  set.seed(seed)
  deltas <- map_dbl(seq_len(n_boot), \(b) {
    idx <- sample.int(nrow(df), size = nrow(df), replace = TRUE)
    dat_b <- df[idx, , drop = FALSE]
    fit_b <- fit_assoc_ns_model_01(dat_b, predictor_col = predictor_col, covars_pref = covars_pref)
    calc_delta_p90_p10_01(fit_b$model, dat_b, predictor_col = predictor_col)
  })
  deltas <- deltas[is.finite(deltas)]
  if (length(deltas) < max(20, n_boot * 0.2)) {
    return(tibble(ci_low = NA_real_, ci_high = NA_real_))
  }

  tibble(
    ci_low = as.numeric(quantile(deltas, 0.025, na.rm = TRUE)),
    ci_high = as.numeric(quantile(deltas, 0.975, na.rm = TRUE))
  )
}

# Compute non-ceiling direct association table for one predictor direction.
calc_assoc_non_ceiling_01 <- function(df_non_boundary, predictor_col, predictor_direction) {
  outcomes_local <- unique(df_non_boundary$outcome)

  map_dfr(outcomes_local, function(outcome_i) {
    dat <- df_non_boundary |>
      filter(outcome == outcome_i) |>
      filter(!is.na(.data[[predictor_col]]), !is.na(y))

    if (nrow(dat) < 10 || dplyr::n_distinct(dat[[predictor_col]]) < 2 || dplyr::n_distinct(dat$y) < 2) {
      return(tibble(
        outcome = outcome_i,
        predictor_direction = predictor_direction,
        n = nrow(dat),
        tau_kendall = NA_real_,
        tau_p = NA_real_,
        rho_spearman = NA_real_,
        rho_p = NA_real_,
        p_model = NA_real_,
        delta_p90_p10 = NA_real_,
        delta_ci_low = NA_real_,
        delta_ci_high = NA_real_,
        covars_used = NA_character_
      ))
    }

    kend <- calc_cor_test_01(dat, predictor_col, method = "kendall")
    spea <- calc_cor_test_01(dat, predictor_col, method = "spearman")

    fit_full <- fit_assoc_ns_model_01(dat, predictor_col = predictor_col, covars_pref = c("age_m0_cov", "exam_order_factor"))
    delta <- calc_delta_p90_p10_01(fit_full$model, dat, predictor_col = predictor_col)
    delta_ci <- bootstrap_delta_p90_p10_01(
      dat,
      predictor_col = predictor_col,
      covars_pref = c("age_m0_cov", "exam_order_factor"),
      n_boot = 200,
      seed = 20260303
    )

    p_model <- NA_real_
    if (!is.null(fit_full$model)) {
      if (length(fit_full$covars_use) > 0) {
        fit_red <- tryCatch(
          lm(as.formula(paste0("y ~ ", paste(fit_full$covars_use, collapse = " + "))), data = dat),
          error = function(e) NULL
        )
      } else {
        fit_red <- tryCatch(lm(y ~ 1, data = dat), error = function(e) NULL)
      }

      if (!is.null(fit_red)) {
        cmp <- tryCatch(anova(fit_red, fit_full$model), error = function(e) NULL)
        if (!is.null(cmp) && nrow(cmp) >= 2) {
          p_model <- cmp$`Pr(>F)`[2]
        }
      }
    }

    tibble(
      outcome = outcome_i,
      predictor_direction = predictor_direction,
      n = nrow(dat),
      tau_kendall = kend$estimate,
      tau_p = kend$p_value,
      rho_spearman = spea$estimate,
      rho_p = spea$p_value,
      p_model = p_model,
      delta_p90_p10 = delta,
      delta_ci_low = delta_ci$ci_low,
      delta_ci_high = delta_ci$ci_high,
      covars_used = paste(fit_full$covars_use, collapse = ", ")
    )
  })
}

# Estimate M0/M6 interaction effects for direct non-ceiling association.
calc_interaction_m0_m6_01 <- function(df_non_boundary, predictor_col, predictor_direction) {
  outcomes_local <- unique(df_non_boundary$outcome)

  map_dfr(outcomes_local, function(outcome_i) {
    dat <- df_non_boundary |>
      filter(outcome == outcome_i, exam_order %in% c("M0", "M6")) |>
      filter(!is.na(.data[[predictor_col]]), !is.na(y)) |>
      mutate(exam_m0m6 = factor(exam_order, levels = c("M0", "M6")))

    if (nrow(dat) < 20 || dplyr::n_distinct(dat$exam_m0m6) < 2 || dplyr::n_distinct(dat[[predictor_col]]) < 2) {
      return(tibble(
        outcome = outcome_i,
        predictor_direction = predictor_direction,
        n = nrow(dat),
        effect_m0 = NA_real_,
        effect_m6 = NA_real_,
        p_interaction = NA_real_
      ))
    }

    covars_use <- get_covars_use_01(dat, c("age_m0_cov"))
    rhs_main <- c(predictor_col, "exam_m0m6", covars_use)
    form_main <- as.formula(paste0("y ~ ", paste(rhs_main, collapse = " + ")))
    form_int <- as.formula(paste0("y ~ ", paste(c(rhs_main, paste0(predictor_col, ":exam_m0m6")), collapse = " + ")))

    fit_main <- tryCatch(lm(form_main, data = dat), error = function(e) NULL)
    fit_int <- tryCatch(lm(form_int, data = dat), error = function(e) NULL)
    if (is.null(fit_main) || is.null(fit_int)) {
      return(tibble(
        outcome = outcome_i,
        predictor_direction = predictor_direction,
        n = nrow(dat),
        effect_m0 = NA_real_,
        effect_m6 = NA_real_,
        p_interaction = NA_real_
      ))
    }

    cf <- coef(fit_int)
    eff_m0 <- ifelse(predictor_col %in% names(cf), unname(cf[[predictor_col]]), NA_real_)
    int_name <- names(cf)[stringr::str_detect(names(cf), paste0("^", predictor_col, ":exam_m0m6M6$|^exam_m0m6M6:", predictor_col, "$"))]
    int_coef <- ifelse(length(int_name) == 1, unname(cf[[int_name]]), NA_real_)
    eff_m6 <- eff_m0 + int_coef

    cmp <- tryCatch(anova(fit_main, fit_int), error = function(e) NULL)
    p_int <- if (!is.null(cmp) && nrow(cmp) >= 2) cmp$`Pr(>F)`[2] else NA_real_

    tibble(
      outcome = outcome_i,
      predictor_direction = predictor_direction,
      n = nrow(dat),
      effect_m0 = eff_m0,
      effect_m6 = eff_m6,
      p_interaction = p_int
    )
  })
}

# Fit boundary (ceiling/floor) probability model with linear/spline candidate forms.
fit_binary_link_01 <- function(train_df,
                               flag_col = "boundary_flag_model",
                               covars = c("exam_time_month", "age_m0_cov")) {
  default_prob <- mean(train_df[[flag_col]], na.rm = TRUE)
  if (!is.finite(default_prob)) {
    default_prob <- 0
  }

  covars_use <- get_covars_use_01(train_df, covars)
  rhs_linear <- c("ig_gradient", covars_use)
  rhs_spline <- c("splines::ns(ig_gradient, df = 3)", covars_use)

  if (nrow(train_df) < 10 ||
    dplyr::n_distinct(train_df$ig_gradient) < 2 ||
    dplyr::n_distinct(train_df[[flag_col]]) < 2) {
    return(list(model = NULL, model_label = "intercept", default_prob = default_prob, covars = covars_use))
  }

  model_lin <- tryCatch(
    glm(
      reformulate(rhs_linear, response = flag_col),
      family = binomial(),
      data = train_df
    ),
    error = function(e) NULL
  )

  model_spline <- NULL
  if (dplyr::n_distinct(train_df$ig_gradient) >= 5) {
    model_spline <- tryCatch(
      glm(
        as.formula(paste0(flag_col, " ~ ", paste(rhs_spline, collapse = " + "))),
        family = binomial(),
        data = train_df
      ),
      error = function(e) NULL
    )
  }

  models <- purrr::compact(list(linear = model_lin, spline = model_spline))
  if (length(models) == 0) {
    return(list(model = NULL, model_label = "intercept", default_prob = default_prob, covars = covars_use))
  }

  aic_vec <- map_dbl(models, ~ tryCatch(AIC(.x), error = function(e) Inf))
  best_name <- names(which.min(aic_vec))[1]

  list(
    model = models[[best_name]],
    model_label = best_name,
    default_prob = default_prob,
    covars = covars_use
  )
}

# Predict boundary probability from a fitted boundary model.
pred_binary_link_01 <- function(fit_obj, new_df) {
  if (is.null(fit_obj$model)) {
    return(rep(fit_obj$default_prob, nrow(new_df)))
  }

  pred <- as.numeric(predict(fit_obj$model, newdata = new_df, type = "response"))
  pred[!is.finite(pred)] <- fit_obj$default_prob
  clamp_01(pred, 1e-6, 1 - 1e-6)
}

# Fit continuous outcome model for the non-boundary part.
fit_cont_link_01 <- function(train_df,
                             covars = c("exam_time_month", "age_m0_cov")) {
  default_value <- mean(train_df$y_model, na.rm = TRUE)
  if (!is.finite(default_value)) {
    default_value <- 0
  }

  covars_use <- get_covars_use_01(train_df, covars)
  rhs_linear <- c("ig_gradient", covars_use)
  rhs_spline <- c("splines::ns(ig_gradient, df = 3)", covars_use)

  if (nrow(train_df) < 5 ||
    dplyr::n_distinct(train_df$ig_gradient) < 2 ||
    dplyr::n_distinct(train_df$y_model) < 2) {
    return(list(model = NULL, model_label = "intercept", default_value = default_value, covars = covars_use))
  }

  model_lin <- tryCatch(
    lm(reformulate(rhs_linear, response = "y_model"), data = train_df),
    error = function(e) NULL
  )

  model_spline <- NULL
  if (dplyr::n_distinct(train_df$ig_gradient) >= 5) {
    model_spline <- tryCatch(
      lm(as.formula(paste0("y_model ~ ", paste(rhs_spline, collapse = " + "))), data = train_df),
      error = function(e) NULL
    )
  }

  models <- purrr::compact(list(linear = model_lin, spline = model_spline))
  if (length(models) == 0) {
    return(list(model = NULL, model_label = "intercept", default_value = default_value, covars = covars_use))
  }

  aic_vec <- map_dbl(models, ~ tryCatch(AIC(.x), error = function(e) Inf))
  best_name <- names(which.min(aic_vec))[1]

  list(
    model = models[[best_name]],
    model_label = best_name,
    default_value = default_value,
    covars = covars_use
  )
}

# Predict continuous non-boundary values from a fitted model.
pred_cont_link_01 <- function(fit_obj, new_df) {
  if (is.null(fit_obj$model)) {
    return(rep(fit_obj$default_value, nrow(new_df)))
  }

  pred <- as.numeric(predict(fit_obj$model, newdata = new_df))
  pred[!is.finite(pred)] <- fit_obj$default_value
  pred
}

# Run full hurdle modeling in grouped CV and return out-of-fold predictions.
run_hurdle_cv_01 <- function(df_outcome,
                             covars = c("exam_time_month", "age_m0_cov")) {
  outcome_name <- unique(df_outcome$outcome)[1]
  boundary_type <- unique(df_outcome$boundary_type)[1]
  min_y <- unique(df_outcome$min_y)[1]
  max_y <- unique(df_outcome$max_y)[1]
  is_haq <- identical(outcome_name, "haq")

  map_dfr(sort(unique(df_outcome$fold)), function(fold_id) {
    train <- df_outcome |> filter(fold != fold_id)
    test <- df_outcome |> filter(fold == fold_id)
    if (nrow(test) == 0 || nrow(train) == 0) {
      return(tibble())
    }

    mean_train_global <- mean(train$y, na.rm = TRUE)
    if (!is.finite(mean_train_global)) {
      mean_train_global <- 0
    }
    mean_by_exam <- train |>
      group_by(exam_order) |>
      summarise(y_hat_base_exam = mean(y, na.rm = TRUE), .groups = "drop")

    covars_in_data <- covars[covars %in% names(train)]
    for (cv in covars_in_data) {
      med <- median(train[[cv]], na.rm = TRUE)
      if (!is.finite(med)) {
        med <- 0
      }
      train[[cv]][is.na(train[[cv]])] <- med
      test[[cv]][is.na(test[[cv]])] <- med
    }

    if (is_haq) {
      train <- train |> mutate(boundary_flag_model = as.integer(y == min_y))
      test <- test |> mutate(boundary_flag_model = as.integer(y == min_y))
      train_cont <- train |>
        filter(y > min_y) |>
        mutate(y_model = y)
    } else {
      train <- train |> mutate(boundary_flag_model = as.integer(y == max_y))
      test <- test |> mutate(boundary_flag_model = as.integer(y == max_y))
      train_cont <- train |>
        filter(y < max_y) |>
        mutate(y_model = y)
    }

    if (nrow(train_cont) == 0) {
      train_cont <- train |> mutate(y_model = y)
    }

    fit_prob <- fit_binary_link_01(train, "boundary_flag_model", covars = covars_in_data)
    fit_cont <- fit_cont_link_01(train_cont, covars = covars_in_data)

    p_boundary <- pred_binary_link_01(fit_prob, test)
    y_hat_cont <- pred_cont_link_01(fit_cont, test)
    y_hat_cont <- clamp_01(y_hat_cont, min_y, max_y)

    if (is_haq) {
      y_hat <- (1 - p_boundary) * y_hat_cont
    } else {
      y_hat <- p_boundary * max_y + (1 - p_boundary) * y_hat_cont
    }
    y_hat <- clamp_01(y_hat, min_y, max_y)

    test_baseline <- test |>
      left_join(mean_by_exam, by = "exam_order") |>
      mutate(
        y_hat_base_global = mean_train_global,
        y_hat_base_exam = if_else(is.na(y_hat_base_exam), mean_train_global, y_hat_base_exam)
      )

    test_baseline |>
      transmute(
        patient_id,
        exam_order,
        exam_order_factor = as.character(exam_order_factor),
        exam_time_month,
        age_m0_cov,
        fold,
        outcome,
        boundary_type,
        min_y,
        max_y,
        ig_gradient,
        ig_grad_pos,
        y,
        y_hat,
        y_hat_base_global,
        y_hat_base_exam,
        p_boundary,
        boundary_flag,
        model_prob = fit_prob$model_label,
        model_cont = fit_cont$model_label
      )
  })
}

# Resolve patient identifier column used for grouped CV and joins.
id_col_01 <- pick_first_col_01(
  d18_ggir_data,
  c("patient_id", "immet_id", "projekt_id", "bbm_id", "sample")
)
# Resolve exam/time column used to derive visit timing.
exam_col_01 <- pick_first_col_01(
  d18_ggir_data,
  c("exam_order", "exam", "poradie_vysetrenia"),
  allow_missing = TRUE
)
# Resolve ig_gradient predictor column.
x_col_01 <- pick_first_col_01(d18_ggir_data, sel_indep_var_ceil_01)
x_lab_01 <- stringr::str_remove(x_col_01, "_enmo$")
# Resolve age-at-baseline covariate column (age_m0 preferred).
age_m0_col_01 <- pick_first_col_01(
  d18_ggir_data,
  c("age_m0", "age", "vek"),
  allow_missing = TRUE
)

# Keep only outcomes present in the current dataset.
outcomes_01 <- intersect(sel_dep_var_ceil_01, names(d18_ggir_data))
if (length(outcomes_01) == 0) {
  stop("No outcome from sel_dep_var_ceil_01 found in d18_ggir_data.")
}

# List optional covariates and detect which are available.
optional_covars_01 <- c("age", "age_m0", "sex", "pohlavi", "bcm", "wear_time", "valid_days")
present_optional_covars_01 <- intersect(optional_covars_01, names(d18_ggir_data))

# Build wide modeling table with ID, time, predictor, age_m0, outcomes, and optional covariates.
d_convert_01 <- d18_ggir_data |>
  mutate(
    patient_id = as.character(.data[[id_col_01]]),
    patient_id = if (stringr::str_detect(id_col_01, "sample")) {
      stringr::str_remove(patient_id, "_[Mm]\\d+$")
    } else {
      patient_id
    },
    exam_order = if (!is.na(exam_col_01)) as.character(.data[[exam_col_01]]) else NA_character_,
    exam_order_std = stringr::str_to_upper(stringr::str_trim(exam_order)),
    exam_time_month = parse_exam_time_01(exam_order),
    exam_order_factor = normalize_exam_factor_01(exam_order),
    ig_gradient = as.numeric(.data[[x_col_01]]),
    ig_grad_pos = -ig_gradient,
    age_raw_cov = if (!is.na(age_m0_col_01)) as.numeric(.data[[age_m0_col_01]]) else NA_real_
  ) |>
  group_by(patient_id) |>
  mutate(
    age_m0_cov = {
      age_m0_vals <- age_raw_cov[exam_order_std == "M0" & !is.na(age_raw_cov)]
      if (length(age_m0_vals) > 0) age_m0_vals[[1]] else NA_real_
    }
  ) |>
  ungroup() |>
  select(
    patient_id,
    exam_order,
    exam_order_factor,
    exam_time_month,
    ig_gradient,
    ig_grad_pos,
    age_m0_cov,
    all_of(outcomes_01),
    any_of(present_optional_covars_01)
  )

# Convert to long format: one row per patient-visit-outcome.
d_convert_long_01 <- d_convert_01 |>
  pivot_longer(
    cols = all_of(outcomes_01),
    names_to = "outcome",
    values_to = "y"
  ) |>
  mutate(y = as.numeric(y))

# QC table with selected columns, classes, and data completeness rates.
res_qc_types_01 <- tibble(
  id_column_used = id_col_01,
  exam_column_used = exam_col_01,
  predictor_column_used = x_col_01,
  age_m0_column_used = age_m0_col_01,
  patient_id_class = class(d_convert_long_01$patient_id)[1],
  exam_order_class = class(d_convert_long_01$exam_order)[1],
  exam_time_month_class = class(d_convert_long_01$exam_time_month)[1],
  ig_gradient_class = class(d_convert_long_01$ig_gradient)[1],
  age_m0_cov_class = class(d_convert_long_01$age_m0_cov)[1],
  y_class = class(d_convert_long_01$y)[1],
  exam_time_non_missing_pct = mean(!is.na(d_convert_long_01$exam_time_month)),
  age_m0_non_missing_pct = mean(!is.na(d_convert_long_01$age_m0_cov)),
  n_rows = nrow(d_convert_long_01),
  n_patients = dplyr::n_distinct(d_convert_long_01$patient_id)
)

# QC table indicating availability of optional covariates.
res_qc_optional_covars_01 <- tibble(
  optional_covariate = optional_covars_01,
  present_in_data = optional_covariate %in% present_optional_covars_01
)

# QC table with missingness profile by outcome.
res_qc_missing_01 <- d_convert_long_01 |>
  group_by(outcome) |>
  summarise(
    n_total = n(),
    n_missing_ig = sum(is.na(ig_gradient)),
    n_missing_exam_time = sum(is.na(exam_time_month)),
    n_missing_age_m0 = sum(is.na(age_m0_cov)),
    n_missing_y = sum(is.na(y)),
    n_valid = sum(!is.na(ig_gradient) & !is.na(y)),
    .groups = "drop"
  )

# Filter to rows with valid ID, predictor, and outcome for modeling.
d_convert_model_01 <- d_convert_long_01 |>
  filter(!is.na(patient_id), !is.na(ig_gradient), !is.na(y))

# QC table reporting ig_gradient outliers (IQR rule) per outcome.
res_qc_outlier_ig_01 <- d_convert_model_01 |>
  group_by(outcome) |>
  summarise(
    q1 = quantile(ig_gradient, 0.25, na.rm = TRUE),
    q3 = quantile(ig_gradient, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower_fence = q1 - 1.5 * iqr,
    upper_fence = q3 + 1.5 * iqr,
    outlier_n = sum(ig_gradient < lower_fence | ig_gradient > upper_fence, na.rm = TRUE),
    outlier_pct = outlier_n / n(),
    .groups = "drop"
  )

# Define default outcome bounds and boundary type (ceiling/floor).
res_bounds_01 <- tibble(
  outcome = c("mmt8", "mmt10", "fi2", "fi3", "haq"),
  min_y = c(0, 0, 0, 0, 0),
  max_y = c(80, 100, 100, 100, 3),
  boundary_type = c("ceiling", "ceiling", "ceiling", "ceiling", "floor")
) |>
  filter(outcome %in% outcomes_01)

# Adjust mmt8 maximum if observed scale indicates 0-800 coding.
if ("mmt8" %in% res_bounds_01$outcome) {
  mmt8_max_obs_01 <- d_convert_model_01 |>
    filter(outcome == "mmt8") |>
    summarise(mx = max(y, na.rm = TRUE)) |>
    pull(mx)

  if (is.finite(mmt8_max_obs_01) && mmt8_max_obs_01 > 100) {
    res_bounds_01 <- res_bounds_01 |>
      mutate(max_y = if_else(outcome == "mmt8", 800, max_y))
  }
}

# Adjust haq maximum if observed values exceed the default upper bound.
if ("haq" %in% res_bounds_01$outcome) {
  haq_max_obs_01 <- d_convert_model_01 |>
    filter(outcome == "haq") |>
    summarise(mx = max(y, na.rm = TRUE)) |>
    pull(mx)

  if (is.finite(haq_max_obs_01) && haq_max_obs_01 > 3) {
    res_bounds_01 <- res_bounds_01 |>
      mutate(max_y = if_else(outcome == "haq", haq_max_obs_01, max_y))
  }
}

# Attach bounds and derive boundary flags used in hurdle modeling.
d_convert_model_01 <- d_convert_model_01 |>
  left_join(res_bounds_01, by = "outcome") |>
  mutate(
    ceiling_y = as.integer(boundary_type == "ceiling" & y == max_y),
    floor_y = as.integer(boundary_type == "floor" & y == min_y),
    boundary_flag = as.integer(
      (boundary_type == "ceiling" & y == max_y) |
        (boundary_type == "floor" & y == min_y)
    ),
    ceiling_haq = as.integer(outcome == "haq" & y == max_y)
  )

# Dataset used for direct association analyses outside ceiling/floor values.
d_non_boundary_01 <- d_convert_model_01 |>
  mutate(
    non_boundary_flag = case_when(
      outcome == "haq" ~ y > min_y,
      TRUE ~ y < max_y
    )
  ) |>
  filter(non_boundary_flag)

# ig_gradient informativeness summary for unique patient-visit rows.
d_ig_unique_01 <- d_convert_01 |>
  distinct(patient_id, exam_order, .keep_all = TRUE)

# Overall informativeness statistics of ig_gradient.
res_ig_info_overall_01 <- d_ig_unique_01 |>
  summarise(
    n = n(),
    mean = mean(ig_gradient, na.rm = TRUE),
    sd = sd(ig_gradient, na.rm = TRUE),
    iqr = IQR(ig_gradient, na.rm = TRUE),
    min = min(ig_gradient, na.rm = TRUE),
    max = max(ig_gradient, na.rm = TRUE)
  )

# Informativeness statistics by visit time.
res_ig_info_by_exam_01 <- d_ig_unique_01 |>
  group_by(exam_order) |>
  summarise(
    n = n(),
    mean = mean(ig_gradient, na.rm = TRUE),
    sd = sd(ig_gradient, na.rm = TRUE),
    iqr = IQR(ig_gradient, na.rm = TRUE),
    min = min(ig_gradient, na.rm = TRUE),
    max = max(ig_gradient, na.rm = TRUE),
    .groups = "drop"
  )

# Correlation of ig_gradient with other activity-related metrics (orientation check).
activity_cols_01 <- names(d18_ggir_data)[
  stringr::str_detect(names(d18_ggir_data), "enmo|mvpa|^ig_")
]
activity_cols_01 <- setdiff(activity_cols_01, c(x_col_01, "ig_gradient"))
activity_cols_01 <- activity_cols_01[map_lgl(activity_cols_01, ~ is.numeric(d18_ggir_data[[.x]]))]

res_ig_corr_activ_01 <- map_dfr(activity_cols_01, function(vv) {
  dat <- d18_ggir_data |>
    transmute(ig_gradient = as.numeric(.data[[x_col_01]]), other = as.numeric(.data[[vv]])) |>
    filter(!is.na(ig_gradient), !is.na(other))

  if (nrow(dat) < 10 || dplyr::n_distinct(dat$ig_gradient) < 2 || dplyr::n_distinct(dat$other) < 2) {
    return(tibble(variable = vv, spearman = NA_real_, p_value = NA_real_, n = nrow(dat)))
  }

  tst <- suppressWarnings(cor.test(dat$ig_gradient, dat$other, method = "spearman", exact = FALSE))
  tibble(variable = vv, spearman = unname(tst$estimate), p_value = tst$p.value, n = nrow(dat))
})

# Short interpretive QC conclusion on whether ig_gradient has enough spread.
res_ig_info_conclusion_01 <- res_ig_info_overall_01 |>
  mutate(
    informativeness_level = case_when(
      sd < 0.10 | iqr < 0.10 ~ "low",
      sd < 0.20 | iqr < 0.20 ~ "moderate",
      TRUE ~ "good"
    ),
    conclusion = case_when(
      informativeness_level == "low" ~ "ig_gradient spread is narrow; strong replacement performance is unlikely.",
      informativeness_level == "moderate" ~ "ig_gradient spread is moderate; partial replacement may be possible in selected outcomes.",
      TRUE ~ "ig_gradient spread is adequate; weak replacement likely reflects model/outcome mismatch rather than range restriction."
    )
  )

# Create grouped patient folds for leakage-safe cross-validation.
fold_map_01 <- make_group_folds_01(d_convert_model_01$patient_id, k = 5, seed = 20260303)

# Add fold assignment to modeling rows.
d_convert_cv_01 <- d_convert_model_01 |>
  left_join(fold_map_01, by = "patient_id")

# Generate out-of-fold hurdle predictions for each outcome.
res_conversion_oof_01 <- d_convert_cv_01 |>
  group_by(outcome) |>
  group_split(.keep = TRUE) |>
  map_dfr(run_hurdle_cv_01)

# Global OOF regression metrics by outcome.
res_metrics_global_01 <- res_conversion_oof_01 |>
  group_by(outcome) |>
  group_modify(~ calc_reg_metrics_01(.x)) |>
  ungroup()

# OOF boundary classification metrics (AUC/PR-AUC) by outcome.
res_metrics_classif_01 <- res_conversion_oof_01 |>
  group_by(outcome, boundary_type) |>
  summarise(
    n = n(),
    prevalence = mean(boundary_flag, na.rm = TRUE),
    auc = calc_auc_01(boundary_flag, p_boundary),
    pr_auc = calc_pr_auc_01(boundary_flag, p_boundary),
    .groups = "drop"
  )

# OOF regression metrics stratified by boundary vs non-boundary subsets.
res_metrics_strata_01 <- res_conversion_oof_01 |>
  mutate(
    stratum = case_when(
      outcome == "haq" & boundary_flag == 1L ~ "haq_eq_0",
      outcome == "haq" & boundary_flag == 0L ~ "haq_gt_0",
      boundary_flag == 1L ~ "ceiling_eq_max",
      TRUE ~ "below_ceiling"
    )
  ) |>
  group_by(outcome, stratum) |>
  group_modify(~ calc_reg_metrics_01(.x)) |>
  ungroup()

# Reliability table: predicted vs observed boundary probabilities by decile.
res_reliability_01 <- res_conversion_oof_01 |>
  mutate(prob_bin = ntile(p_boundary, 10)) |>
  group_by(outcome, prob_bin) |>
  summarise(
    n = n(),
    mean_pred_prob = mean(p_boundary, na.rm = TRUE),
    observed_rate = mean(boundary_flag, na.rm = TRUE),
    .groups = "drop"
  )

# Named list of reliability plots, one per outcome.
res_reliability_plot_01 <- res_reliability_01 |>
  group_nest(outcome) |>
  mutate(
    plot = map2(data, outcome, \(dat, outcome_nm) {
      ggplot(dat, aes(x = mean_pred_prob, y = observed_rate)) +
        geom_point(size = 2) +
        geom_line() +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
        coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        labs(
          title = paste("Reliability:", outcome_nm),
          x = "Predicted boundary probability",
          y = "Observed boundary frequency"
        )
    })
  ) |>
  select(outcome, plot) |>
  tibble::deframe()

# BA-like working table with diff/mean and clinically relevant strata.
res_ba_like_01 <- res_conversion_oof_01 |>
  mutate(
    ba_diff = y - y_hat,
    ba_mean = (y + y_hat) / 2,
    stratum = case_when(
      outcome == "haq" & boundary_flag == 1L ~ "haq_eq_0",
      outcome == "haq" & boundary_flag == 0L ~ "haq_gt_0",
      boundary_flag == 1L ~ "ceiling_eq_max",
      TRUE ~ "below_ceiling"
    )
  )

# BA-like summary statistics by outcome.
res_ba_like_stats_01 <- res_ba_like_01 |>
  group_by(outcome) |>
  group_modify(~ calc_ba_stats_01(.x)) |>
  ungroup()

# BA-like summary statistics by outcome and boundary stratum.
res_ba_like_stats_strata_01 <- res_ba_like_01 |>
  group_by(outcome, stratum) |>
  group_modify(~ calc_ba_stats_01(.x)) |>
  ungroup()

# Named list of BA-like plots, one per outcome.
res_ba_like_plot_01 <- res_ba_like_01 |>
  group_nest(outcome) |>
  mutate(
    plot = map2(data, outcome, \(dat, outcome_nm) {
      stat <- calc_ba_stats_01(dat)
      ggplot(dat, aes(x = ba_mean, y = ba_diff)) +
        geom_point(alpha = 0.65) +
        geom_hline(yintercept = stat$bias, linetype = "solid") +
        geom_hline(yintercept = stat$loa_lower, linetype = "dashed") +
        geom_hline(yintercept = stat$loa_upper, linetype = "dashed") +
        geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
        labs(
          title = paste("BA-like (OOF):", outcome_nm),
          x = "Mean of observed and predicted",
          y = "Observed - predicted"
        )
    })
  ) |>
  select(outcome, plot) |>
  tibble::deframe()

# Direct non-ceiling association table for both ig_gradient directions.
res_assoc_non_ceiling_01 <- bind_rows(
  calc_assoc_non_ceiling_01(d_non_boundary_01, predictor_col = "ig_gradient", predictor_direction = "ig_gradient"),
  calc_assoc_non_ceiling_01(d_non_boundary_01, predictor_col = "ig_grad_pos", predictor_direction = "ig_grad_pos")
)

# Interaction table (M0 vs M6) for both ig_gradient directions.
res_interaction_m0m6_01 <- bind_rows(
  calc_interaction_m0_m6_01(d_non_boundary_01, predictor_col = "ig_gradient", predictor_direction = "ig_gradient"),
  calc_interaction_m0_m6_01(d_non_boundary_01, predictor_col = "ig_grad_pos", predictor_direction = "ig_grad_pos")
)

# Direction consistency summary against expected clinical sign.
res_ig_direction_summary_01 <- res_assoc_non_ceiling_01 |>
  mutate(
    expected_sign = if_else(outcome == "haq", -1, 1),
    observed_sign = sign(delta_p90_p10),
    sign_consistent = observed_sign == expected_sign
  ) |>
  group_by(predictor_direction) |>
  summarise(
    n_outcomes = n(),
    n_consistent = sum(sign_consistent, na.rm = TRUE),
    consistency_pct = n_consistent / n_outcomes,
    mean_abs_delta = mean(abs(delta_p90_p10), na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(consistency_pct), desc(mean_abs_delta))

# Recommended direction for interpretation based on consistency and effect size.
res_ig_direction_choice_01 <- res_ig_direction_summary_01 |>
  slice(1) |>
  mutate(
    recommendation = case_when(
      predictor_direction == "ig_grad_pos" ~ "Use ig_grad_pos for intuitive interpretation (higher value = clinically better in most outcomes).",
      TRUE ~ "Keep original ig_gradient direction; sign consistency is not improved by inversion."
    )
  )

# Benchmark table: compare model error vs simple baseline predictions.
res_benchmark_gain_01 <- bind_rows(
  res_conversion_oof_01 |> mutate(subset = "global"),
  res_conversion_oof_01 |> filter(boundary_flag == 0) |> mutate(subset = "non_boundary")
) |>
  group_by(outcome, subset) |>
  summarise(
    n = n(),
    mae_model = mean(abs(y - y_hat), na.rm = TRUE),
    rmse_model = sqrt(mean((y - y_hat)^2, na.rm = TRUE)),
    mae_base_global = mean(abs(y - y_hat_base_global), na.rm = TRUE),
    rmse_base_global = sqrt(mean((y - y_hat_base_global)^2, na.rm = TRUE)),
    mae_base_exam = mean(abs(y - y_hat_base_exam), na.rm = TRUE),
    rmse_base_exam = sqrt(mean((y - y_hat_base_exam)^2, na.rm = TRUE)),
    gain_mae_vs_global = mae_base_global - mae_model,
    gain_rmse_vs_global = rmse_base_global - rmse_model,
    gain_mae_vs_exam = mae_base_exam - mae_model,
    gain_rmse_vs_exam = rmse_base_exam - rmse_model,
    .groups = "drop"
  )

# Export-ready diagnostic figures focused on non-ceiling behavior and boundary separation.
diag_plot_date_01 <- format(Sys.Date(), "%y%m%d")
res_diag_fig_manifest_01 <- map_dfr(outcomes_01, function(outcome_i) {
  dat_nb <- d_non_boundary_01 |>
    filter(outcome == outcome_i)

  dat_boundary <- d_convert_model_01 |>
    filter(outcome == outcome_i) |>
    mutate(boundary_group = if_else(boundary_flag == 1L, "boundary", "non_boundary"))

  files_saved <- c()

  if (nrow(dat_nb) >= 10) {
    p_scatter_exam <- ggplot(
      dat_nb |> filter(exam_order %in% c("M0", "M6")),
      aes(x = ig_gradient, y = y, color = exam_order)
    ) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "lm", formula = y ~ splines::ns(x, df = 3), se = TRUE) +
      labs(
        title = paste0(outcome_i, ": non-boundary scatter by exam"),
        x = "ig_gradient",
        y = outcome_i
      )
    fp1 <- file.path(output_fig_dir_ceil_01, paste0(diag_plot_date_01, "_", outcome_i, "_nonboundary_scatter_exam_01.png"))
    if (save_plot_safe_01(p_scatter_exam, fp1)) files_saved <- c(files_saved, fp1)

    dat_q <- dat_nb |>
      mutate(y_quartile = factor(ntile(y, 4), labels = c("Q1", "Q2", "Q3", "Q4")))
    p_quartile_density <- ggplot(dat_q, aes(x = ig_gradient, fill = y_quartile)) +
      geom_density(alpha = 0.35) +
      labs(
        title = paste0(outcome_i, ": ig_gradient density by outcome quartiles (non-boundary)"),
        x = "ig_gradient",
        y = "Density",
        fill = "Outcome quartile"
      )
    fp2 <- file.path(output_fig_dir_ceil_01, paste0(diag_plot_date_01, "_", outcome_i, "_nonboundary_quartile_density_01.png"))
    if (save_plot_safe_01(p_quartile_density, fp2)) files_saved <- c(files_saved, fp2)

    p_scatter_pos <- ggplot(dat_nb, aes(x = ig_grad_pos, y = y)) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "lm", formula = y ~ splines::ns(x, df = 3), se = TRUE, color = "steelblue") +
      labs(
        title = paste0(outcome_i, ": non-boundary scatter (ig_grad_pos)"),
        x = "ig_grad_pos",
        y = outcome_i
      )
    fp3 <- file.path(output_fig_dir_ceil_01, paste0(diag_plot_date_01, "_", outcome_i, "_nonboundary_scatter_pos_01.png"))
    if (save_plot_safe_01(p_scatter_pos, fp3)) files_saved <- c(files_saved, fp3)
  }

  if (nrow(dat_boundary) >= 10 && dplyr::n_distinct(dat_boundary$boundary_group) > 1) {
    p_boundary_density <- ggplot(dat_boundary, aes(x = ig_gradient, fill = boundary_group)) +
      geom_density(alpha = 0.35) +
      labs(
        title = paste0(outcome_i, ": ig_gradient density boundary vs non-boundary"),
        x = "ig_gradient",
        y = "Density",
        fill = "Group"
      )
    fp4 <- file.path(output_fig_dir_ceil_01, paste0(diag_plot_date_01, "_", outcome_i, "_boundary_density_01.png"))
    if (save_plot_safe_01(p_boundary_density, fp4)) files_saved <- c(files_saved, fp4)
  }

  tibble(
    outcome = outcome_i,
    n_figures_saved = length(files_saved),
    figure_files = paste(files_saved, collapse = " | ")
  )
})

# Optional sensitivity analysis: fixed-effects vs mixed-effects (random intercept by patient).
res_mixed_sensitivity_01 <- {
  predictor_choice <- if (nrow(res_ig_direction_choice_01) > 0) res_ig_direction_choice_01$predictor_direction[1] else "ig_gradient"
  predictor_col <- ifelse(predictor_choice == "ig_grad_pos", "ig_grad_pos", "ig_gradient")
  lme4_available <- requireNamespace("lme4", quietly = TRUE)

  map_dfr(outcomes_01, function(outcome_i) {
    dat <- d_non_boundary_01 |>
      filter(outcome == outcome_i) |>
      filter(!is.na(.data[[predictor_col]]), !is.na(y))

    if (nrow(dat) < 20 || dplyr::n_distinct(dat$patient_id) < 5 || dplyr::n_distinct(dat[[predictor_col]]) < 2) {
      return(tibble(
        outcome = outcome_i,
        predictor_direction = predictor_choice,
        n = nrow(dat),
        fixed_estimate = NA_real_,
        fixed_ci_low = NA_real_,
        fixed_ci_high = NA_real_,
        mixed_estimate = NA_real_,
        mixed_ci_low = NA_real_,
        mixed_ci_high = NA_real_,
        mixed_model_used = FALSE,
        note = "Insufficient data for sensitivity model."
      ))
    }

    rhs <- c(predictor_col)
    if ("exam_order_factor" %in% names(dat) && dplyr::n_distinct(dat$exam_order_factor) > 1) rhs <- c(rhs, "exam_order_factor")
    if ("age_m0_cov" %in% names(dat) && dplyr::n_distinct(dat$age_m0_cov[!is.na(dat$age_m0_cov)]) > 1) rhs <- c(rhs, "age_m0_cov")
    form_fix <- as.formula(paste0("y ~ ", paste(rhs, collapse = " + ")))

    fit_fix <- tryCatch(lm(form_fix, data = dat), error = function(e) NULL)
    if (is.null(fit_fix)) {
      return(tibble(
        outcome = outcome_i,
        predictor_direction = predictor_choice,
        n = nrow(dat),
        fixed_estimate = NA_real_,
        fixed_ci_low = NA_real_,
        fixed_ci_high = NA_real_,
        mixed_estimate = NA_real_,
        mixed_ci_low = NA_real_,
        mixed_ci_high = NA_real_,
        mixed_model_used = FALSE,
        note = "Fixed model failed."
      ))
    }

    cf_fix <- coef(summary(fit_fix))
    fix_est <- ifelse(predictor_col %in% rownames(cf_fix), cf_fix[predictor_col, "Estimate"], NA_real_)
    ci_fix <- tryCatch(confint(fit_fix, parm = predictor_col, level = 0.95), error = function(e) c(NA_real_, NA_real_))

    if (!lme4_available) {
      return(tibble(
        outcome = outcome_i,
        predictor_direction = predictor_choice,
        n = nrow(dat),
        fixed_estimate = fix_est,
        fixed_ci_low = ci_fix[1],
        fixed_ci_high = ci_fix[2],
        mixed_estimate = NA_real_,
        mixed_ci_low = NA_real_,
        mixed_ci_high = NA_real_,
        mixed_model_used = FALSE,
        note = "lme4 not available."
      ))
    }

    form_mix <- as.formula(paste0("y ~ ", paste(rhs, collapse = " + "), " + (1|patient_id)"))
    fit_mix <- tryCatch(lme4::lmer(form_mix, data = dat, REML = FALSE), error = function(e) NULL)
    if (is.null(fit_mix)) {
      return(tibble(
        outcome = outcome_i,
        predictor_direction = predictor_choice,
        n = nrow(dat),
        fixed_estimate = fix_est,
        fixed_ci_low = ci_fix[1],
        fixed_ci_high = ci_fix[2],
        mixed_estimate = NA_real_,
        mixed_ci_low = NA_real_,
        mixed_ci_high = NA_real_,
        mixed_model_used = FALSE,
        note = "Mixed model failed."
      ))
    }

    fixef_mix <- lme4::fixef(fit_mix)
    mix_est <- ifelse(predictor_col %in% names(fixef_mix), unname(fixef_mix[[predictor_col]]), NA_real_)
    ci_mix <- tryCatch(suppressMessages(confint(fit_mix, parm = predictor_col, method = "Wald")), error = function(e) c(NA_real_, NA_real_))
    if (is.matrix(ci_mix)) ci_mix <- as.numeric(ci_mix[1, ])

    tibble(
      outcome = outcome_i,
      predictor_direction = predictor_choice,
      n = nrow(dat),
      fixed_estimate = fix_est,
      fixed_ci_low = ci_fix[1],
      fixed_ci_high = ci_fix[2],
      mixed_estimate = mix_est,
      mixed_ci_low = ci_mix[1],
      mixed_ci_high = ci_mix[2],
      mixed_model_used = TRUE,
      note = "OK"
    )
  })
}

# Extract non-boundary performance subset used in decision summary.
res_non_boundary_01 <- res_metrics_strata_01 |>
  filter(
    (outcome == "haq" & stratum == "haq_gt_0") |
      (outcome != "haq" & stratum == "below_ceiling")
  ) |>
  select(
    outcome,
    n_non_boundary = n,
    mae_non_boundary = mae,
    rmse_non_boundary = rmse
  )

# Extract boundary-only performance subset used in decision summary.
res_boundary_only_01 <- res_metrics_strata_01 |>
  filter(
    (outcome == "haq" & stratum == "haq_eq_0") |
      (outcome != "haq" & stratum == "ceiling_eq_max")
  ) |>
  select(
    outcome,
    n_boundary = n,
    mae_boundary = mae,
    rmse_boundary = rmse
  )
# **********************************************************************
## Final table ----
# **********************************************************************
# Final decision table combining accuracy, boundary, and BA-like criteria.
res_conversion_decision_01 <- res_bounds_01 |>
  left_join(res_non_boundary_01, by = "outcome") |>
  left_join(res_boundary_only_01, by = "outcome") |>
  left_join(
    res_metrics_classif_01 |>
      select(outcome, auc_boundary = auc, pr_auc_boundary = pr_auc),
    by = "outcome"
  ) |>
  left_join(
    res_ba_like_stats_01 |>
      select(outcome, bias, trend_slope, trend_p),
    by = "outcome"
  ) |>
  mutate(
    range_y = max_y - min_y,
    non_boundary_ok = !is.na(mae_non_boundary) & mae_non_boundary <= 0.10 * range_y,
    boundary_prob_ok = !is.na(auc_boundary) & auc_boundary >= 0.70,
    ba_bias_ok = !is.na(bias) & abs(bias) <= 0.10 * range_y,
    ba_trend_ok = is.na(trend_p) | trend_p >= 0.05,
    usability = case_when(
      non_boundary_ok & boundary_prob_ok & ba_bias_ok & ba_trend_ok ~ "usable",
      non_boundary_ok & (boundary_prob_ok | ba_bias_ok) ~ "partially_usable",
      TRUE ~ "limited"
    ),
    boundary_comment = case_when(
      is.na(auc_boundary) ~ "Boundary probability not estimable.",
      auc_boundary >= 0.80 ~ "Boundary probability estimated well.",
      auc_boundary >= 0.70 ~ "Boundary probability estimated moderately.",
      TRUE ~ "Boundary probability estimated weakly."
    )
  ) |>
  select(
    outcome, boundary_type, min_y, max_y, range_y,
    n_non_boundary, mae_non_boundary, rmse_non_boundary,
    n_boundary, mae_boundary, rmse_boundary,
    auc_boundary, pr_auc_boundary,
    bias, trend_slope, trend_p,
    non_boundary_ok, boundary_prob_ok, ba_bias_ok, ba_trend_ok,
    usability, boundary_comment
  )

# **********************************************************************
## Export ----
# **********************************************************************
rio::export(
  list(
    qc_types = janitor::clean_names(res_qc_types_01),
    qc_missing = janitor::clean_names(res_qc_missing_01),
    qc_ig_outlier = janitor::clean_names(res_qc_outlier_ig_01),
    qc_optional_covars = janitor::clean_names(res_qc_optional_covars_01),
    bounds = janitor::clean_names(res_bounds_01),
    oof_predictions = janitor::clean_names(res_conversion_oof_01),
    metrics_global = janitor::clean_names(res_metrics_global_01),
    metrics_classif = janitor::clean_names(res_metrics_classif_01),
    metrics_strata = janitor::clean_names(res_metrics_strata_01),
    assoc_non_ceiling = janitor::clean_names(res_assoc_non_ceiling_01),
    interaction_m0_m6 = janitor::clean_names(res_interaction_m0m6_01),
    ig_direction = janitor::clean_names(res_ig_direction_summary_01),
    ig_direction_choice = janitor::clean_names(res_ig_direction_choice_01),
    benchmark_gain = janitor::clean_names(res_benchmark_gain_01),
    reliability = janitor::clean_names(res_reliability_01),
    ba_like = janitor::clean_names(res_ba_like_stats_01),
    ba_like_strata = janitor::clean_names(res_ba_like_stats_strata_01),
    ig_info_overall = janitor::clean_names(res_ig_info_overall_01),
    ig_info_by_exam = janitor::clean_names(res_ig_info_by_exam_01),
    ig_corr_activ = janitor::clean_names(res_ig_corr_activ_01),
    ig_info_conclusion = janitor::clean_names(res_ig_info_conclusion_01),
    diag_fig_manifest = janitor::clean_names(res_diag_fig_manifest_01),
    mixed_sensitivity = janitor::clean_names(res_mixed_sensitivity_01),
    decision = janitor::clean_names(res_conversion_decision_01)
  ),
  file.path(
    output_tab_dir_ceil_01,
    paste0(format(Sys.Date(), "%y%m%d"), "_ig_gradient_hurdle_conversion_01.xlsx")
  )
)

# **********************************************************************
# 5. Disease activity in ceiling/floor ----
# **********************************************************************
# Resolve anchor activity column (PhVAS preferred; anchor_activity fallback).
anchor_col_01 <- pick_first_col_01(
  d18_ggir_data,
  c("ph_vas", "PhVAS", "phvas", "anchor_activity"),
  allow_missing = TRUE
)

# Build anchor table on visit level and align IDs/exam coding with conversion dataset.
d_anchor_activity_01 <- d18_ggir_data |>
  mutate(
    patient_id = as.character(.data[[id_col_01]]),
    patient_id = if (stringr::str_detect(id_col_01, "sample")) {
      stringr::str_remove(patient_id, "_[Mm]\\d+$")
    } else {
      patient_id
    },
    exam_order = if (!is.na(exam_col_01)) as.character(.data[[exam_col_01]]) else NA_character_,
    anchor_activity = if (!is.na(anchor_col_01)) as.numeric(.data[[anchor_col_01]]) else NA_real_
  ) |>
  group_by(patient_id, exam_order) |>
  summarise(
    anchor_activity = {
      vals <- anchor_activity[!is.na(anchor_activity)]
      if (length(vals) > 0) vals[[1]] else NA_real_
    },
    .groups = "drop"
  )
# **********************************************************************
## No ceiling/floor cohort ----
# **********************************************************************
# Join anchor to long conversion data and keep only no-ceiling/floor rows.
d_disease_activity_01 <- d_convert_model_01 |>
  left_join(d_anchor_activity_01, by = c("patient_id", "exam_order")) |>
  filter(
    boundary_flag == 0L,
    !is.na(anchor_activity),
    !is.na(y),
    !is.na(ig_gradient),
    !is.na(patient_id)
  ) |>
  mutate(patient_id = factor(patient_id))

# Display labels for model comparison plots/tables.
disease_model_labels_01 <- c(
  model_a = "Model 1: clin",
  model_b = paste0("Model 2: clin + ", x_lab_01),
  model_c = paste0("Model 3: ", x_lab_01)
)

make_models_named_abc_01 <- function(fit_obj) {
  out <- list(
    fit_obj$model_a,
    fit_obj$model_b,
    fit_obj$model_c
  )
  names(out) <- unname(disease_model_labels_01[c("model_a", "model_b", "model_c")])
  out
}

# Ensure all displayed model labels use the currently selected activity parameter.
normalize_model_display_labels_01 <- function(x_chr) {
  out <- as.character(x_chr)
  out <- stringr::str_replace_all(out, fixed("Model 2: clin + ig_gradient"), disease_model_labels_01[["model_b"]])
  out <- stringr::str_replace_all(out, fixed("Model 3: ig_gradient"), disease_model_labels_01[["model_c"]])
  out <- stringr::str_replace_all(out, fixed("ig_gradient"), x_lab_01)
  out
}

# Recode compare_performance auto-names (model1/model2/...) to custom labels.
relabel_compare_performance_01 <- function(perf_obj, model_labels) {
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
  perf_obj[[nm_col]] <- normalize_model_display_labels_01(new_vals)

  # Some plotting methods use additional character columns for labels.
  chr_cols <- names(perf_obj)[vapply(perf_obj, is.character, logical(1))]
  for (cc in chr_cols) {
    perf_obj[[cc]] <- normalize_model_display_labels_01(perf_obj[[cc]])
  }

  perf_obj
}

# Safely fit one mixed model.
fit_lmer_safe_01 <- function(formula_obj, data) {
  tryCatch(
    lme4::lmer(
      formula_obj,
      data = data,
      REML = FALSE
    ),
    error = function(e) NULL
  )
}

# Build a fixed-effects table with estimate, p-value, and 95% CI.
# Uses broom.mixed when available, with a safe fallback from lme4 output.
tidy_lmer_fixed_01 <- function(fit_obj) {
  if (is.null(fit_obj)) {
    return(
      tibble(
        term = NA_character_,
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        note = "Model fit unavailable."
      )
    )
  }

  tidy_tbl <- if (requireNamespace("broom.mixed", quietly = TRUE)) {
    tryCatch(
      broom.mixed::tidy(
        fit_obj,
        effects = "fixed",
        conf.int = TRUE,
        conf.level = 0.95
      ),
      error = function(e) NULL
    )
  } else {
    NULL
  }

  if (!is.null(tidy_tbl) && nrow(tidy_tbl) > 0) {
    nm <- names(tidy_tbl)
    se_col <- intersect(c("std.error", "std_error"), nm)
    stat_col <- intersect(c("statistic", "t.value", "z.value"), nm)
    p_col <- intersect(c("p.value", "p_value"), nm)
    low_col <- intersect(c("conf.low", "conf_low"), nm)
    high_col <- intersect(c("conf.high", "conf_high"), nm)

    se_col <- if (length(se_col) > 0) se_col[[1]] else NA_character_
    stat_col <- if (length(stat_col) > 0) stat_col[[1]] else NA_character_
    p_col <- if (length(p_col) > 0) p_col[[1]] else NA_character_
    low_col <- if (length(low_col) > 0) low_col[[1]] else NA_character_
    high_col <- if (length(high_col) > 0) high_col[[1]] else NA_character_

    out <- tidy_tbl |>
      transmute(
        term = as.character(term),
        estimate = as.numeric(estimate),
        std_error = if (!is.na(se_col)) as.numeric(.data[[se_col]]) else NA_real_,
        statistic = if (!is.na(stat_col)) as.numeric(.data[[stat_col]]) else NA_real_,
        p_value = if (!is.na(p_col)) as.numeric(.data[[p_col]]) else NA_real_,
        conf_low = if (!is.na(low_col)) as.numeric(.data[[low_col]]) else NA_real_,
        conf_high = if (!is.na(high_col)) as.numeric(.data[[high_col]]) else NA_real_
      ) |>
      mutate(
        p_value = coalesce(
          p_value,
          2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
        ),
        note = "OK"
      )

    return(out)
  }

  coef_mat <- coef(summary(fit_obj))
  coef_tbl <- as_tibble(coef_mat, rownames = "term")

  est_col <- intersect(c("Estimate", "estimate"), names(coef_tbl))
  se_col <- intersect(c("Std. Error", "std.error", "std_error"), names(coef_tbl))
  stat_col <- intersect(c("t value", "z value", "statistic"), names(coef_tbl))
  p_col <- intersect(c("Pr(>|t|)", "Pr(>|z|)", "p.value", "p_value"), names(coef_tbl))

  est_col <- if (length(est_col) > 0) est_col[[1]] else NA_character_
  se_col <- if (length(se_col) > 0) se_col[[1]] else NA_character_
  stat_col <- if (length(stat_col) > 0) stat_col[[1]] else NA_character_
  p_col <- if (length(p_col) > 0) p_col[[1]] else NA_character_

  out <- coef_tbl |>
    transmute(
      term = as.character(term),
      estimate = if (!is.na(est_col)) as.numeric(.data[[est_col]]) else NA_real_,
      std_error = if (!is.na(se_col)) as.numeric(.data[[se_col]]) else NA_real_,
      statistic = if (!is.na(stat_col)) as.numeric(.data[[stat_col]]) else NA_real_,
      p_value = if (!is.na(p_col)) as.numeric(.data[[p_col]]) else NA_real_
    ) |>
    mutate(
      p_value = coalesce(
        p_value,
        2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
      )
    )

  parm_fix <- names(lme4::fixef(fit_obj))
  ci_fix <- tryCatch(
    suppressMessages(confint(fit_obj, parm = parm_fix, method = "Wald")),
    error = function(e) NULL
  )

  if (!is.null(ci_fix)) {
    ci_tbl <- as_tibble(as.data.frame(ci_fix), rownames = "term")
    names(ci_tbl)[2:3] <- c("conf_low", "conf_high")
    out <- out |>
      left_join(ci_tbl, by = "term")
  } else {
    out <- out |>
      mutate(conf_low = NA_real_, conf_high = NA_real_)
  }

  out |>
    mutate(note = "OK")
}

# Build one table summarizing fixed effects of Model 2 across outcomes.
build_model2_summary_table_01 <- function(res_models) {
  imap_dfr(
    res_models,
    function(x, nm) {
      fit_b <- if (!is.null(x$fit)) x$fit$model_b else NULL
      tidy_lmer_fixed_01(fit_b) |>
        mutate(
          outcome = nm,
          model = "model_b",
          model_label = disease_model_labels_01[["model_b"]],
          n = x$n,
          n_patients = x$n_patients,
          .before = 1
        )
    }
  )
}

# Build one table describing model formulas and selected covariates by outcome.
build_model_formula_table_01 <- function(res_models) {
  imap_dfr(
    res_models,
    function(x, nm) {
      fit_obj <- x$fit
      tibble(
        outcome = nm,
        n = x$n,
        n_patients = x$n_patients,
        formula_a = if (!is.null(fit_obj)) paste(fit_obj$formula_a, collapse = " ") else NA_character_,
        formula_b = if (!is.null(fit_obj)) paste(fit_obj$formula_b, collapse = " ") else NA_character_,
        formula_c = if (!is.null(fit_obj)) paste(fit_obj$formula_c, collapse = " ") else NA_character_,
        covariates = if (!is.null(fit_obj) && length(fit_obj$covars_use) > 0) {
          paste(fit_obj$covars_use, collapse = ", ")
        } else {
          NA_character_
        },
        note = x$note
      )
    }
  )
}

# Build post-hoc table for Model 3 exam contrasts in All cohort.
# Contrast estimate is computed as exam_to - exam_from.
build_model3_posthoc_pairs_01 <- function(res_models, data_all, pairs_tbl) {
  imap_dfr(
    res_models,
    function(x, nm) {
      fit_c <- if (!is.null(x$fit)) x$fit$model_c else NULL

      dat_outcome <- data_all |>
        filter(outcome == nm, !is.na(exam_order_factor))

      fit_frame <- if (!is.null(fit_c)) {
        tryCatch(as_tibble(stats::model.frame(fit_c)), error = function(e) tibble())
      } else {
        tibble()
      }

      exam_present <- if ("exam_order_factor" %in% names(fit_frame)) {
        fit_frame |>
          distinct(exam_order_factor) |>
          mutate(exam_order_factor = as.character(exam_order_factor)) |>
          pull(exam_order_factor)
      } else {
        dat_outcome |>
          distinct(exam_order_factor) |>
          mutate(exam_order_factor = as.character(exam_order_factor)) |>
          pull(exam_order_factor)
      }

      pairs_use <- pairs_tbl |>
        filter(exam_from %in% exam_present, exam_to %in% exam_present)

      if (nrow(pairs_use) == 0) {
        return(tibble())
      }

      boundary_exam <- dat_outcome |>
        group_by(exam_order_factor) |>
        summarise(
          has_boundary = any(boundary_flag == 1L, na.rm = TRUE),
          n_boundary = sum(boundary_flag == 1L, na.rm = TRUE),
          .groups = "drop"
        ) |>
        mutate(exam_order_factor = as.character(exam_order_factor))

      pair_meta <- pairs_use |>
        left_join(
          boundary_exam |>
            rename(
              exam_from = exam_order_factor,
              has_boundary_from = has_boundary,
              n_boundary_from = n_boundary
            ),
          by = "exam_from"
        ) |>
        left_join(
          boundary_exam |>
            rename(
              exam_to = exam_order_factor,
              has_boundary_to = has_boundary,
              n_boundary_to = n_boundary
            ),
          by = "exam_to"
        ) |>
        mutate(
          pair_has_boundary = coalesce(has_boundary_from, FALSE) | coalesce(has_boundary_to, FALSE)
        )

      if (is.null(fit_c)) {
        return(
          pair_meta |>
            transmute(
              outcome = nm,
              model = "model_c",
              model_label = disease_model_labels_01[["model_c"]],
              pair,
              exam_from,
              exam_to,
              estimate = NA_real_,
              std_error = NA_real_,
              conf_low = NA_real_,
              conf_high = NA_real_,
              p_value = NA_real_,
              has_boundary_from,
              has_boundary_to,
              pair_has_boundary,
              n_boundary_from,
              n_boundary_to,
              n = x$n,
              n_patients = x$n_patients,
              note = "Model 3 unavailable."
            )
        )
      }

      if (!("exam_order_factor" %in% names(fit_frame))) {
        return(
          pair_meta |>
            transmute(
              outcome = nm,
              model = "model_c",
              model_label = disease_model_labels_01[["model_c"]],
              pair,
              exam_from,
              exam_to,
              estimate = NA_real_,
              std_error = NA_real_,
              conf_low = NA_real_,
              conf_high = NA_real_,
              p_value = NA_real_,
              has_boundary_from,
              has_boundary_to,
              pair_has_boundary,
              n_boundary_from,
              n_boundary_to,
              n = x$n,
              n_patients = x$n_patients,
              note = "Model 3 has no exam term."
            )
        )
      }

      exam_levels_fit <- levels(factor(fit_frame$exam_order_factor))
      if (length(exam_levels_fit) < 2) {
        return(
          pair_meta |>
            transmute(
              outcome = nm,
              model = "model_c",
              model_label = disease_model_labels_01[["model_c"]],
              pair,
              exam_from,
              exam_to,
              estimate = NA_real_,
              std_error = NA_real_,
              conf_low = NA_real_,
              conf_high = NA_real_,
              p_value = NA_real_,
              has_boundary_from,
              has_boundary_to,
              pair_has_boundary,
              n_boundary_from,
              n_boundary_to,
              n = x$n,
              n_patients = x$n_patients,
              note = "Model 3 has fewer than two exam levels."
            )
        )
      }

      fixef_vec <- lme4::fixef(fit_c)
      vc_fix <- tryCatch(as.matrix(stats::vcov(fit_c)), error = function(e) NULL)

      if (is.null(vc_fix)) {
        return(
          pair_meta |>
            transmute(
              outcome = nm,
              model = "model_c",
              model_label = disease_model_labels_01[["model_c"]],
              pair,
              exam_from,
              exam_to,
              estimate = NA_real_,
              std_error = NA_real_,
              conf_low = NA_real_,
              conf_high = NA_real_,
              p_value = NA_real_,
              has_boundary_from,
              has_boundary_to,
              pair_has_boundary,
              n_boundary_from,
              n_boundary_to,
              n = x$n,
              n_patients = x$n_patients,
              note = "Model 3 vcov unavailable."
            )
        )
      }

      get_exam_coef_name_01 <- function(exam_level, base_level) {
        if (is.na(exam_level) || is.na(base_level) || identical(exam_level, base_level)) {
          return(NA_character_)
        }
        paste0("exam_order_factor", exam_level)
      }

      base_level <- exam_levels_fit[[1]]

      contrast_tbl <- pmap_dfr(
        list(pair_meta$pair, pair_meta$exam_from, pair_meta$exam_to),
        \(pair_i, exam_from_i, exam_to_i) {
          coef_from <- get_exam_coef_name_01(exam_from_i, base_level)
          coef_to <- get_exam_coef_name_01(exam_to_i, base_level)

          coef_ok <- function(coef_name) {
            is.na(coef_name) || coef_name %in% names(fixef_vec)
          }

          if (!coef_ok(coef_from) || !coef_ok(coef_to)) {
            return(
              tibble(
                pair = pair_i,
                estimate = NA_real_,
                std_error = NA_real_,
                conf_low = NA_real_,
                conf_high = NA_real_,
                p_value = NA_real_,
                note = "Missing exam coefficient in Model 3."
              )
            )
          }

          est_from <- if (is.na(coef_from)) 0 else unname(fixef_vec[[coef_from]])
          est_to <- if (is.na(coef_to)) 0 else unname(fixef_vec[[coef_to]])
          est <- est_to - est_from

          var_from <- if (is.na(coef_from)) 0 else vc_fix[coef_from, coef_from]
          var_to <- if (is.na(coef_to)) 0 else vc_fix[coef_to, coef_to]
          cov_from_to <- if (is.na(coef_from) || is.na(coef_to)) {
            0
          } else {
            vc_fix[coef_from, coef_to]
          }

          var_diff <- as.numeric(var_to + var_from - 2 * cov_from_to)
          var_diff <- if (is.finite(var_diff)) max(var_diff, 0) else NA_real_
          se <- if (is.finite(var_diff)) sqrt(var_diff) else NA_real_
          z_stat <- if (is.finite(se) && se > 0) est / se else NA_real_
          p_val <- if (is.finite(z_stat)) 2 * stats::pnorm(abs(z_stat), lower.tail = FALSE) else NA_real_
          conf_low <- if (is.finite(se)) est - 1.96 * se else NA_real_
          conf_high <- if (is.finite(se)) est + 1.96 * se else NA_real_

          tibble(
            pair = pair_i,
            estimate = as.numeric(est),
            std_error = as.numeric(se),
            conf_low = as.numeric(conf_low),
            conf_high = as.numeric(conf_high),
            p_value = as.numeric(p_val),
            note = "OK"
          )
        }
      )

      pair_meta |>
        left_join(contrast_tbl, by = "pair") |>
        transmute(
          outcome = nm,
          model = "model_c",
          model_label = disease_model_labels_01[["model_c"]],
          pair,
          exam_from,
          exam_to,
          estimate,
          std_error,
          conf_low,
          conf_high,
          p_value,
          has_boundary_from,
          has_boundary_to,
          pair_has_boundary,
          n_boundary_from,
          n_boundary_to,
          n = x$n,
          n_patients = x$n_patients,
          note = coalesce(note, "Post-hoc contrast unavailable.")
        )
    }
  )
}

# Fit disease-activity models A/B/C for one clinical outcome.
fit_disease_models_one_outcome_01 <- function(dat_outcome) {
  covars_pref <- c("age_m0_cov", "sex", "pohlavi", "bcm", "wear_time", "valid_days")
  covars_use <- get_covars_use_01(dat_outcome, covars_pref)

  exam_term <- if ("exam_order_factor" %in% names(dat_outcome) &&
    dplyr::n_distinct(dat_outcome$exam_order_factor) > 1) {
    "exam_order_factor"
  } else {
    NULL
  }

  rhs_a <- c("y", exam_term, covars_use)

  rhs_b <- c("y", "ig_gradient", exam_term, covars_use)
  rhs_c <- c("ig_gradient", exam_term, covars_use)

  form_a <- as.formula(paste0("anchor_activity ~ ", paste(c(rhs_a, "(1|patient_id)"), collapse = " + ")))
  form_b <- as.formula(paste0("anchor_activity ~ ", paste(c(rhs_b, "(1|patient_id)"), collapse = " + ")))
  form_c <- as.formula(paste0("anchor_activity ~ ", paste(c(rhs_c, "(1|patient_id)"), collapse = " + ")))

  fit_a <- fit_lmer_safe_01(form_a, dat_outcome)
  fit_b <- fit_lmer_safe_01(form_b, dat_outcome)
  fit_c <- fit_lmer_safe_01(form_c, dat_outcome)

  list(
    model_a = fit_a,
    model_b = fit_b,
    model_c = fit_c,
    formula_a = format(form_a),
    formula_b = format(form_b),
    formula_c = format(form_c),
    covars_use = covars_use
  )
}

# Run A/B/C models for each clinical outcome in no-ceiling/floor cohort.
res_disease_models_01 <- outcomes_01 |>
  set_names() |>
  map(function(outcome_i) {
    dat <- d_disease_activity_01 |>
      filter(outcome == outcome_i)

    if (nrow(dat) < 20 || dplyr::n_distinct(dat$patient_id) < 5) {
      return(
        list(
          outcome = outcome_i,
          n = nrow(dat),
          n_patients = dplyr::n_distinct(dat$patient_id),
          fit = NULL,
          compare = NULL,
          note = "Insufficient no-ceiling/floor data for mixed modeling."
        )
      )
    }

    fit_obj <- fit_disease_models_one_outcome_01(dat)

    models_named <- make_models_named_abc_01(fit_obj)
    models_named <- models_named[!map_lgl(models_named, is.null)]

    display_labels <- names(models_named)
    perf_tbl <- NULL
    if (length(models_named) >= 2 && requireNamespace("performance", quietly = TRUE)) {
      perf_tbl <- do.call(performance::compare_performance, models_named) |>
        relabel_compare_performance_01(model_labels = display_labels) |>
        as_tibble() |>
        mutate(outcome = outcome_i, .before = 1)
    }

    list(
      outcome = outcome_i,
      n = nrow(dat),
      n_patients = dplyr::n_distinct(dat$patient_id),
      fit = fit_obj,
      compare = perf_tbl,
      note = if (is.null(perf_tbl)) "Model comparison unavailable (missing package or <2 valid models)." else "OK"
    )
  })

# **********************************************************************
### Tables ----
# **********************************************************************
# Flat comparison table for all outcomes (A/B/C via compare_performance).
res_disease_model_compare_01 <- map_dfr(res_disease_models_01, "compare")

# Lightweight run summary for diagnostics.
res_disease_model_status_01 <- imap_dfr(
  res_disease_models_01,
  function(x, nm) {
    tibble(
      outcome = nm,
      n = x$n,
      n_patients = x$n_patients,
      has_model_a = !is.null(x$fit$model_a),
      has_model_b = !is.null(x$fit$model_b),
      has_model_c = !is.null(x$fit$model_c),
      model_a_label = disease_model_labels_01[["model_a"]],
      model_b_label = disease_model_labels_01[["model_b"]],
      model_c_label = disease_model_labels_01[["model_c"]],
      has_compare = !is.null(x$compare),
      note = x$note
    )
  }
)

# Model formulas and selected covariates by outcome.
res_disease_model_formula_01 <- build_model_formula_table_01(res_disease_models_01)

# Fixed-effects summary for Model 2 with p-value and 95% CI.
res_disease_model2_tidy_01 <- build_model2_summary_table_01(res_disease_models_01)

# **********************************************************************
### Figures ----
# **********************************************************************
# Build one performance plot per outcome via see::plot(compare_performance(...)).
build_disease_perf_plot_01 <- function(model_entry, outcome_nm) {
  if (
    is.null(model_entry$fit) ||
      !requireNamespace("performance", quietly = TRUE) ||
      !requireNamespace("see", quietly = TRUE)
  ) {
    return(NULL)
  }

  models_named <- make_models_named_abc_01(model_entry$fit)
  models_named <- models_named[!map_lgl(models_named, is.null)]

  if (length(models_named) < 2) {
    return(NULL)
  }

  cmp_obj <- do.call(performance::compare_performance, models_named) |>
    relabel_compare_performance_01(model_labels = names(models_named))

  p_obj <- tryCatch(
    plot(cmp_obj),
    error = function(e) NULL
  )

  if (is.null(p_obj) || !inherits(p_obj, "ggplot")) {
    return(NULL)
  }

  p_obj +
    ggplot2::labs(title = paste0("Outcome: ", outcome_nm)) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.title = ggplot2::element_text(size = 11, face = "bold"),
      axis.text = ggplot2::element_text(size = 9, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 12, hjust = 1),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 8, face = "bold")
    )
}

# Build one sjPlot prediction figure per outcome from Model 3
# (PhVAS ~ ig_gradient + exam + covariates + (1|patient_id)).
build_disease_ig_exam_plot_01 <- function(model_entry, outcome_nm, cohort_subtitle) {
  if (
    is.null(model_entry$fit) ||
      is.null(model_entry$fit$model_c) ||
      !requireNamespace("sjPlot", quietly = TRUE)
  ) {
    return(NULL)
  }

  fit_c <- model_entry$fit$model_c
  terms_use <- c("ig_gradient")

  p_obj <- tryCatch(
    sjPlot::plot_model(fit_c, type = "pred", terms = terms_use),
    error = function(e) NULL
  )

  if (is.null(p_obj) || !inherits(p_obj, "ggplot")) {
    return(NULL)
  }

  raw_dat <- tryCatch(
    as_tibble(stats::model.frame(fit_c)),
    error = function(e) NULL
  )

  has_raw <- !is.null(raw_dat) &&
    all(c("ig_gradient", "anchor_activity") %in% names(raw_dat))
  raw_layer <- if (has_raw) {
    ggplot2::geom_point(
      data = raw_dat,
      mapping = ggplot2::aes(x = ig_gradient, y = anchor_activity),
      inherit.aes = FALSE,
      color = "grey45",
      alpha = 0.35,
      size = 2.1
    )
  } else {
    NULL
  }

  p_obj +
    raw_layer +
    ggplot2::labs(
      title = paste0("Outcome: ", outcome_nm),
      subtitle = cohort_subtitle,
      x = x_lab_01,
      y = "PhVAS"
    ) +
    ggplot2::theme_bw(base_size = 20) +
    ggplot2::guides(
      color = ggplot2::guide_legend(title = NULL),
      fill = ggplot2::guide_legend(title = NULL),
      linetype = ggplot2::guide_legend(title = NULL)
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.title = ggplot2::element_text(size = 26, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 18, face = "bold"),
      axis.title = ggplot2::element_text(size = 22, face = "bold"),
      axis.text = ggplot2::element_text(size = 18, face = "bold"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 17, face = "bold")
    )
}

res_disease_perf_plot_list_01 <- imap(res_disease_models_01, build_disease_perf_plot_01) |>
  purrr::compact()

# Combine outcome-level performance plots into one multi-figure.
if (length(res_disease_perf_plot_list_01) > 0) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    res_disease_perf_multi_plot_01 <- patchwork::wrap_plots(
      res_disease_perf_plot_list_01,
      ncol = 2
    ) +
      patchwork::plot_annotation(
        title = "Disease Activity Models: compare_performance (Model 1/2/3)",
        subtitle = "No ceiling/floor cohort",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 18, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = 12, face = "bold")
        )
      ) &
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 11, face = "bold"),
        axis.text = ggplot2::element_text(size = 9, face = "bold"),
        strip.text = ggplot2::element_text(size = 11, face = "bold"),
        legend.title = ggplot2::element_text(size = 10, face = "bold"),
        legend.text = ggplot2::element_text(size = 8, face = "bold")
      )
  } else {
    # Fallback: faceted bar chart of Performance_Score from compare_performance.
    model_nm_col_01 <- intersect(c("Name", "Model", "name", "model"), names(res_disease_model_compare_01))
    model_nm_col_01 <- if (length(model_nm_col_01) > 0) model_nm_col_01[[1]] else NA_character_

    if (!is.na(model_nm_col_01) && "Performance_Score" %in% names(res_disease_model_compare_01)) {
      res_disease_perf_multi_plot_01 <- res_disease_model_compare_01 |>
        mutate(model_name = as.character(.data[[model_nm_col_01]])) |>
        ggplot(aes(x = model_name, y = Performance_Score, fill = model_name)) +
        geom_col(width = 0.7, show.legend = FALSE) +
        facet_wrap(~ outcome, scales = "free_y") +
        labs(
          title = "Disease Activity Models: compare_performance (Model 1/2/3)",
          subtitle = "No ceiling/floor cohort",
          x = "Model",
          y = "Performance score"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          plot.title = element_text(size = 18, face = "bold"),
          plot.subtitle = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 9, face = "bold"),
          strip.text = element_text(size = 11, face = "bold")
        )
    } else {
      res_disease_perf_multi_plot_01 <- NULL
    }
  }
} else {
  res_disease_perf_multi_plot_01 <- NULL
}

# Save one multi-figure with model performance comparison.
res_disease_perf_plot_file_01 <- file.path(
  output_fig_dir_ceil_01,
  paste0(format(Sys.Date(), "%y%m%d"), "_disease_activity_compare_performance_01.png")
)

res_disease_perf_plot_saved_01 <- if (!is.null(res_disease_perf_multi_plot_01)) {
  save_plot_safe_01(
    plot_obj = res_disease_perf_multi_plot_01,
    file_path = res_disease_perf_plot_file_01,
    width = 20,
    height = 12,
    dpi = 300
  )
} else {
  FALSE
}

res_disease_perf_plot_manifest_01 <- tibble(
  file_path = res_disease_perf_plot_file_01,
  saved = res_disease_perf_plot_saved_01,
  n_outcomes_plotted = length(res_disease_perf_plot_list_01)
)

# Build and export one multi-figure with Model 3 prediction curves
# (x = ig_gradient, y = PhVAS, groups = exam) for no-ceiling/floor cohort.
res_disease_ig_exam_plot_list_01 <- imap(
  res_disease_models_01,
  \(x, nm) build_disease_ig_exam_plot_01(x, nm, "No ceiling/floor cohort")
) |>
  purrr::compact()

if (length(res_disease_ig_exam_plot_list_01) > 0 && requireNamespace("patchwork", quietly = TRUE)) {
  res_disease_ig_exam_multi_plot_01 <- patchwork::wrap_plots(
    res_disease_ig_exam_plot_list_01,
    ncol = 2
  ) +
    patchwork::plot_annotation(
      title = paste0("PhVAS vs ", x_lab_01, " (Model 3)"),
      subtitle = "No ceiling/floor cohort",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 34, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 24, face = "bold")
      )
    ) &
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 22, face = "bold"),
      axis.text = ggplot2::element_text(size = 18, face = "bold"),
      strip.text = ggplot2::element_text(size = 24, face = "bold"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 17, face = "bold")
    )
} else {
  res_disease_ig_exam_multi_plot_01 <- NULL
}

res_disease_ig_exam_plot_file_01 <- file.path(
  output_fig_dir_ceil_01,
  paste0(format(Sys.Date(), "%y%m%d"), "_disease_activity_phvas_ig_gradient_exam_no_cei_flo_01.png")
)

res_disease_ig_exam_plot_saved_01 <- if (!is.null(res_disease_ig_exam_multi_plot_01)) {
  save_plot_safe_01(
    plot_obj = res_disease_ig_exam_multi_plot_01,
    file_path = res_disease_ig_exam_plot_file_01,
    width = 16,
    height = 10,
    dpi = 300
  )
} else {
  FALSE
}

res_disease_ig_exam_plot_manifest_01 <- tibble(
  file_path = res_disease_ig_exam_plot_file_01,
  saved = res_disease_ig_exam_plot_saved_01,
  n_outcomes_plotted = length(res_disease_ig_exam_plot_list_01)
)

# **********************************************************************
## All cohort ----
# **********************************************************************
# Join anchor to long conversion data and keep all rows (including boundary rows).
d_disease_activity_all_01 <- d_convert_model_01 |>
  left_join(d_anchor_activity_01, by = c("patient_id", "exam_order")) |>
  filter(
    !is.na(anchor_activity),
    !is.na(y),
    !is.na(ig_gradient),
    !is.na(patient_id)
  ) |>
  mutate(patient_id = factor(patient_id))

# Run A/B/C models for each clinical outcome in full cohort.
res_disease_models_all_01 <- outcomes_01 |>
  set_names() |>
  map(function(outcome_i) {
    dat <- d_disease_activity_all_01 |>
      filter(outcome == outcome_i)

    if (nrow(dat) < 20 || dplyr::n_distinct(dat$patient_id) < 5) {
      return(
        list(
          outcome = outcome_i,
          n = nrow(dat),
          n_patients = dplyr::n_distinct(dat$patient_id),
          fit = NULL,
          compare = NULL,
          note = "Insufficient all-cohort data for mixed modeling."
        )
      )
    }

    fit_obj <- fit_disease_models_one_outcome_01(dat)

    models_named <- make_models_named_abc_01(fit_obj)
    models_named <- models_named[!map_lgl(models_named, is.null)]

    display_labels <- names(models_named)
    perf_tbl <- NULL
    if (length(models_named) >= 2 && requireNamespace("performance", quietly = TRUE)) {
      perf_tbl <- do.call(performance::compare_performance, models_named) |>
        relabel_compare_performance_01(model_labels = display_labels) |>
        as_tibble() |>
        mutate(outcome = outcome_i, .before = 1)
    }

    list(
      outcome = outcome_i,
      n = nrow(dat),
      n_patients = dplyr::n_distinct(dat$patient_id),
      fit = fit_obj,
      compare = perf_tbl,
      note = if (is.null(perf_tbl)) "Model comparison unavailable (missing package or <2 valid models)." else "OK"
    )
  })
# **********************************************************************
### Tables ----
# **********************************************************************
# Flat comparison table for all outcomes in full cohort.
res_disease_model_compare_all_01 <- map_dfr(res_disease_models_all_01, "compare")

# Lightweight run summary for diagnostics in full cohort.
res_disease_model_status_all_01 <- imap_dfr(
  res_disease_models_all_01,
  function(x, nm) {
    tibble(
      outcome = nm,
      n = x$n,
      n_patients = x$n_patients,
      has_model_a = !is.null(x$fit$model_a),
      has_model_b = !is.null(x$fit$model_b),
      has_model_c = !is.null(x$fit$model_c),
      model_a_label = disease_model_labels_01[["model_a"]],
      model_b_label = disease_model_labels_01[["model_b"]],
      model_c_label = disease_model_labels_01[["model_c"]],
      has_compare = !is.null(x$compare),
      note = x$note
    )
  }
)

# Model formulas and selected covariates by outcome (all cohort).
res_disease_model_formula_all_01 <- build_model_formula_table_01(res_disease_models_all_01)

# Fixed-effects summary for Model 2 with p-value and 95% CI (all cohort).
res_disease_model2_tidy_all_01 <- build_model2_summary_table_01(res_disease_models_all_01)

# **********************************************************************
### Post hoc ----
# **********************************************************************
# Requested Model 3 exam-pair contrasts (computed as exam_to - exam_from).
res_posthoc_pairs_model3_all_01 <- tibble(
  pair = c("M0-M3", "M0-M6", "M0-M18"),
  exam_from = c("M0", "M0", "M0"),
  exam_to = c("M3", "M6", "M18")
)

# Post-hoc comparison table for Model 3 in all cohort.
res_disease_model3_posthoc_all_01 <- build_model3_posthoc_pairs_01(
  res_models = res_disease_models_all_01,
  data_all = d_disease_activity_all_01,
  pairs_tbl = res_posthoc_pairs_model3_all_01
)

# **********************************************************************
### Figures ----
# **********************************************************************
# Build and combine performance plots for full cohort.
res_disease_perf_plot_list_all_01 <- imap(res_disease_models_all_01, build_disease_perf_plot_01) |>
  purrr::compact()

if (length(res_disease_perf_plot_list_all_01) > 0) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    res_disease_perf_multi_plot_all_01 <- patchwork::wrap_plots(
      res_disease_perf_plot_list_all_01,
      ncol = 2
    ) +
      patchwork::plot_annotation(
        title = "Disease Activity Models: compare_performance (Model 1/2/3)",
        subtitle = "All cohort",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 18, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = 12, face = "bold")
        )
      ) &
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 11, face = "bold"),
        axis.text = ggplot2::element_text(size = 9, face = "bold"),
        strip.text = ggplot2::element_text(size = 11, face = "bold"),
        legend.title = ggplot2::element_text(size = 10, face = "bold"),
        legend.text = ggplot2::element_text(size = 8, face = "bold")
      )
  } else {
    model_nm_col_all_01 <- intersect(c("Name", "Model", "name", "model"), names(res_disease_model_compare_all_01))
    model_nm_col_all_01 <- if (length(model_nm_col_all_01) > 0) model_nm_col_all_01[[1]] else NA_character_

    if (!is.na(model_nm_col_all_01) && "Performance_Score" %in% names(res_disease_model_compare_all_01)) {
      res_disease_perf_multi_plot_all_01 <- res_disease_model_compare_all_01 |>
        mutate(model_name = as.character(.data[[model_nm_col_all_01]])) |>
        ggplot(aes(x = model_name, y = Performance_Score, fill = model_name)) +
        geom_col(width = 0.7, show.legend = FALSE) +
        facet_wrap(~ outcome, scales = "free_y") +
        labs(
          title = "Disease Activity Models: compare_performance (Model 1/2/3)",
          subtitle = "All cohort",
          x = "Model",
          y = "Performance score"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          plot.title = element_text(size = 18, face = "bold"),
          plot.subtitle = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 9, face = "bold"),
          strip.text = element_text(size = 11, face = "bold")
        )
    } else {
      res_disease_perf_multi_plot_all_01 <- NULL
    }
  }
} else {
  res_disease_perf_multi_plot_all_01 <- NULL
}

# Save one multi-figure with model performance comparison for full cohort.
res_disease_perf_plot_file_all_01 <- file.path(
  output_fig_dir_ceil_01,
  paste0(format(Sys.Date(), "%y%m%d"), "_disease_activity_compare_performance_all_cohort_01.png")
)

res_disease_perf_plot_saved_all_01 <- if (!is.null(res_disease_perf_multi_plot_all_01)) {
  save_plot_safe_01(
    plot_obj = res_disease_perf_multi_plot_all_01,
    file_path = res_disease_perf_plot_file_all_01,
    width = 20,
    height = 12,
    dpi = 300
  )
} else {
  FALSE
}

res_disease_perf_plot_manifest_all_01 <- tibble(
  file_path = res_disease_perf_plot_file_all_01,
  saved = res_disease_perf_plot_saved_all_01,
  n_outcomes_plotted = length(res_disease_perf_plot_list_all_01)
)

# Build and export one multi-figure with Model 3 prediction curves
# (x = ig_gradient, y = PhVAS, groups = exam) for all cohort.
res_disease_ig_exam_plot_list_all_01 <- imap(
  res_disease_models_all_01,
  \(x, nm) build_disease_ig_exam_plot_01(x, nm, "All cohort")
) |>
  purrr::compact()

if (length(res_disease_ig_exam_plot_list_all_01) > 0 && requireNamespace("patchwork", quietly = TRUE)) {
  res_disease_ig_exam_multi_plot_all_01 <- patchwork::wrap_plots(
    res_disease_ig_exam_plot_list_all_01,
    ncol = 2
  ) +
    patchwork::plot_annotation(
      title = paste0("PhVAS vs ", x_lab_01, " (Model 3)"),
      subtitle = "All cohort",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 34, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 24, face = "bold")
      )
    ) &
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 22, face = "bold"),
      axis.text = ggplot2::element_text(size = 18, face = "bold"),
      strip.text = ggplot2::element_text(size = 24, face = "bold"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 17, face = "bold")
    )
} else {
  res_disease_ig_exam_multi_plot_all_01 <- NULL
}

res_disease_ig_exam_plot_file_all_01 <- file.path(
  output_fig_dir_ceil_01,
  paste0(format(Sys.Date(), "%y%m%d"), "_disease_activity_phvas_ig_gradient_exam_all_cohort_01.png")
)

res_disease_ig_exam_plot_saved_all_01 <- if (!is.null(res_disease_ig_exam_multi_plot_all_01)) {
  save_plot_safe_01(
    plot_obj = res_disease_ig_exam_multi_plot_all_01,
    file_path = res_disease_ig_exam_plot_file_all_01,
    width = 16,
    height = 10,
    dpi = 300
  )
} else {
  FALSE
}

res_disease_ig_exam_plot_manifest_all_01 <- tibble(
  file_path = res_disease_ig_exam_plot_file_all_01,
  saved = res_disease_ig_exam_plot_saved_all_01,
  n_outcomes_plotted = length(res_disease_ig_exam_plot_list_all_01)
)

# **********************************************************************
## Ceiling/floor cohort ----
# **********************************************************************
### Post hoc ----
# **********************************************************************
# Requested Model 3 exam-pair contrasts (computed as exam_to - exam_from).
res_posthoc_pairs_model3_cei_flo_01 <- tibble(
  pair = c("M0-M3", "M0-M6", "M0-M18"),
  exam_from = c("M0", "M0", "M0"),
  exam_to = c("M3", "M6", "M18")
)

# Keep only patients who had at least one ceiling/floor value for each outcome.
res_boundary_patients_by_outcome_cei_flo_01 <- d_convert_model_01 |>
  filter(boundary_flag == 1L, !is.na(patient_id)) |>
  transmute(outcome, patient_id = as.character(patient_id)) |>
  distinct()

# Full pool for post-hoc metadata (keeps boundary rows to report boundary counts in pairs).
d_disease_activity_cei_flo_pool_01 <- d_convert_model_01 |>
  left_join(d_anchor_activity_01, by = c("patient_id", "exam_order")) |>
  filter(
    !is.na(anchor_activity),
    !is.na(y),
    !is.na(ig_gradient),
    !is.na(patient_id)
  ) |>
  mutate(
    patient_id_chr = as.character(patient_id),
    patient_id = factor(patient_id)
  ) |>
  semi_join(
    res_boundary_patients_by_outcome_cei_flo_01,
    by = c("outcome", "patient_id_chr" = "patient_id")
  ) |>
  select(-patient_id_chr)

# Fit models only on no-ceiling/floor rows of boundary-patient cohort.
d_disease_activity_cei_flo_nonboundary_01 <- d_disease_activity_cei_flo_pool_01 |>
  filter(boundary_flag == 0L) |>
  mutate(patient_id = factor(patient_id))

res_disease_models_cei_flo_01 <- outcomes_01 |>
  set_names() |>
  map(function(outcome_i) {
    dat <- d_disease_activity_cei_flo_nonboundary_01 |>
      filter(outcome == outcome_i)

    if (nrow(dat) < 20 || dplyr::n_distinct(dat$patient_id) < 5) {
      return(
        list(
          outcome = outcome_i,
          n = nrow(dat),
          n_patients = dplyr::n_distinct(dat$patient_id),
          fit = NULL,
          compare = NULL,
          note = "Insufficient ceiling/floor-cohort non-boundary data for mixed modeling."
        )
      )
    }

    fit_obj <- fit_disease_models_one_outcome_01(dat)

    models_named <- make_models_named_abc_01(fit_obj)
    models_named <- models_named[!map_lgl(models_named, is.null)]

    display_labels <- names(models_named)
    perf_tbl <- NULL
    if (length(models_named) >= 2 && requireNamespace("performance", quietly = TRUE)) {
      perf_tbl <- do.call(performance::compare_performance, models_named) |>
        relabel_compare_performance_01(model_labels = display_labels) |>
        as_tibble() |>
        mutate(outcome = outcome_i, .before = 1)
    }

    list(
      outcome = outcome_i,
      n = nrow(dat),
      n_patients = dplyr::n_distinct(dat$patient_id),
      fit = fit_obj,
      compare = perf_tbl,
      note = if (is.null(perf_tbl)) "Model comparison unavailable (missing package or <2 valid models)." else "OK"
    )
  })

# Post-hoc comparison table for Model 3 in ceiling/floor cohort.
res_disease_model3_posthoc_cei_flo_01 <- build_model3_posthoc_pairs_01(
  res_models = res_disease_models_cei_flo_01,
  data_all = d_disease_activity_cei_flo_pool_01,
  pairs_tbl = res_posthoc_pairs_model3_cei_flo_01
)

# **********************************************************************
### Export ----
# **********************************************************************
# Export disease-activity model tables for no-ceiling/floor, all, and ceiling/floor cohorts.
rio::export(
  list(
    no_cei_flo_model_compare = janitor::clean_names(res_disease_model_compare_01),
    no_cei_flo_model_status = janitor::clean_names(res_disease_model_status_01),
    no_cei_flo_model_formulas = janitor::clean_names(res_disease_model_formula_01),
    no_cei_flo_model2_fixed_effects = janitor::clean_names(res_disease_model2_tidy_01),
    ceiling_flo_cohort_model3_posthoc = janitor::clean_names(res_disease_model3_posthoc_cei_flo_01),
    no_cei_flo_perf_plot_manifest = janitor::clean_names(res_disease_perf_plot_manifest_01),
    no_cei_flo_ig_plot_manifest = janitor::clean_names(res_disease_ig_exam_plot_manifest_01),
    all_cohort_model_compare = janitor::clean_names(res_disease_model_compare_all_01),
    all_cohort_model_status = janitor::clean_names(res_disease_model_status_all_01),
    all_cohort_model_formulas = janitor::clean_names(res_disease_model_formula_all_01),
    all_cohort_model2_fixed_effects = janitor::clean_names(res_disease_model2_tidy_all_01),
    all_cohort_model3_posthoc = janitor::clean_names(res_disease_model3_posthoc_all_01),
    all_cohort_perf_plot_manifest = janitor::clean_names(res_disease_perf_plot_manifest_all_01),
    all_cohort_ig_plot_manifest = janitor::clean_names(res_disease_ig_exam_plot_manifest_all_01)
  ),
  file.path(
    output_tab_dir_ceil_01,
    paste0(format(Sys.Date(), "%y%m%d"), "_disease_activity_models_01.xlsx")
  )
)

# Export Model 3 post-hoc table also as a standalone file for direct use.
rio::export(
  janitor::clean_names(res_disease_model3_posthoc_all_01),
  file.path(
    output_tab_dir_ceil_01,
    paste0(format(Sys.Date(), "%y%m%d"), "_res_disease_model3_posthoc_all_01.xlsx")
  )
)

# Export Model 3 post-hoc table for ceiling/floor cohort as standalone file.
rio::export(
  janitor::clean_names(res_disease_model3_posthoc_cei_flo_01),
  file.path(
    output_tab_dir_ceil_01,
    paste0(format(Sys.Date(), "%y%m%d"), "_res_disease_model3_posthoc_cei_flo_01.xlsx")
  )
)

# **********************************************************************
## Ceiling/floor cohort: ceiling-only PhVAS models ----
# **********************************************************************
##  PhVAS ~ ig_gradient + covariates 
# Boundary lookup with fallback when current session has incomplete res_bounds_01.
boundary_lookup_01 <- if (
  exists("res_bounds_01") &&
    is.data.frame(res_bounds_01) &&
    all(c("outcome", "boundary_type") %in% names(res_bounds_01))
) {
  res_bounds_01 |>
    select(outcome, boundary_type) |>
    distinct()
} else {
  tibble(
    outcome = sel_dep_var_ceil_01,
    boundary_type = case_when(
      outcome %in% c("mmt8", "mmt10", "fi2", "fi3") ~ "ceiling",
      outcome == "haq" ~ "floor",
      TRUE ~ NA_character_
    )
  )
}

# Prepare boundary metadata and keep only ceiling rows for each clinical outcome.
d_disease_ceiling_phvas_01 <- d_convert_model_01 |>
  left_join(d_anchor_activity_01, by = c("patient_id", "exam_order")) |>
  left_join(boundary_lookup_01, by = "outcome", suffix = c("", "_lookup")) |>
  mutate(boundary_type = coalesce(boundary_type, boundary_type_lookup)) |>
  select(-any_of("boundary_type_lookup")) |>
  filter(
    boundary_type == "ceiling",
    boundary_flag == 1L,
    !is.na(anchor_activity),
    !is.na(ig_gradient),
    !is.na(patient_id)
  ) |>
  mutate(patient_id = factor(patient_id))

# Fit one PhVAS model in ceiling-only cohort for a chosen clinical parameter.
run_phvas_ceiling_model_01 <- function(outcome_i) {
  boundary_type_i <- boundary_lookup_01 |>
    filter(outcome == outcome_i) |>
    pull(boundary_type)
  boundary_type_i <- if (length(boundary_type_i) > 0) boundary_type_i[[1]] else NA_character_

  dat <- d_disease_ceiling_phvas_01 |>
    filter(outcome == outcome_i)

  if (!identical(boundary_type_i, "ceiling")) {
    return(
      list(
        fit = NULL,
        summary = tibble(
          outcome = outcome_i,
          boundary_type = boundary_type_i,
          n = 0L,
          n_patients = 0L,
          ig_beta = NA_real_,
          ig_se = NA_real_,
          ig_p_value = NA_real_,
          ig_ci_low = NA_real_,
          ig_ci_high = NA_real_,
          aic = NA_real_,
          bic = NA_real_,
          r2_marginal = NA_real_,
          r2_conditional = NA_real_,
          formula = NA_character_,
          covariates = NA_character_,
          note = "Outcome has no ceiling definition."
        ),
        coef = tibble(),
        note = "Outcome has no ceiling definition."
      )
    )
  }

  if (nrow(dat) < 20 || dplyr::n_distinct(dat$patient_id) < 5 || dplyr::n_distinct(dat$ig_gradient) < 2) {
    return(
      list(
        fit = NULL,
        summary = tibble(
          outcome = outcome_i,
          boundary_type = boundary_type_i,
          n = nrow(dat),
          n_patients = dplyr::n_distinct(dat$patient_id),
          ig_beta = NA_real_,
          ig_se = NA_real_,
          ig_p_value = NA_real_,
          ig_ci_low = NA_real_,
          ig_ci_high = NA_real_,
          aic = NA_real_,
          bic = NA_real_,
          r2_marginal = NA_real_,
          r2_conditional = NA_real_,
          formula = NA_character_,
          covariates = NA_character_,
          note = "Insufficient ceiling-only data for mixed modeling."
        ),
        coef = tibble(),
        note = "Insufficient ceiling-only data for mixed modeling."
      )
    )
  }

  covars_pref <- c("age_m0_cov", "sex", "pohlavi", "bcm", "wear_time", "valid_days")
  covars_use <- get_covars_use_01(dat, covars_pref)
  exam_term <- if ("exam_order_factor" %in% names(dat) &&
    dplyr::n_distinct(dat$exam_order_factor) > 1) {
    "exam_order_factor"
  } else {
    NULL
  }

  rhs <- c("ig_gradient", exam_term, covars_use)
  form <- as.formula(paste0("anchor_activity ~ ", paste(c(rhs, "(1|patient_id)"), collapse = " + ")))
  fit <- fit_lmer_safe_01(form, dat)

  if (is.null(fit)) {
    return(
      list(
        fit = NULL,
        summary = tibble(
          outcome = outcome_i,
          boundary_type = boundary_type_i,
          n = nrow(dat),
          n_patients = dplyr::n_distinct(dat$patient_id),
          ig_beta = NA_real_,
          ig_se = NA_real_,
          ig_p_value = NA_real_,
          ig_ci_low = NA_real_,
          ig_ci_high = NA_real_,
          aic = NA_real_,
          bic = NA_real_,
          r2_marginal = NA_real_,
          r2_conditional = NA_real_,
          formula = paste(deparse(form), collapse = " "),
          covariates = if (length(covars_use) > 0) paste(covars_use, collapse = ", ") else NA_character_,
          note = "Model fit failed."
        ),
        coef = tibble(),
        note = "Model fit failed."
      )
    )
  }

  coef_mat <- coef(summary(fit))
  coef_tbl_raw <- as_tibble(coef_mat, rownames = "term")
  est_col <- intersect(c("Estimate", "estimate"), names(coef_tbl_raw))
  se_col <- intersect(c("Std. Error", "std.error", "std_error"), names(coef_tbl_raw))
  stat_col <- intersect(c("t value", "z value", "statistic"), names(coef_tbl_raw))
  p_col <- intersect(c("Pr(>|t|)", "Pr(>|z|)", "p.value", "p_value"), names(coef_tbl_raw))
  est_name <- if (length(est_col) > 0) est_col[[1]] else NA_character_
  se_name <- if (length(se_col) > 0) se_col[[1]] else NA_character_
  stat_name <- if (length(stat_col) > 0) stat_col[[1]] else NA_character_
  p_name <- if (length(p_col) > 0) p_col[[1]] else NA_character_

  coef_tbl <- coef_tbl_raw |>
    transmute(term) |>
    mutate(
      estimate = if (!is.na(est_name)) as.numeric(coef_tbl_raw[[est_name]]) else NA_real_,
      std_error = if (!is.na(se_name)) as.numeric(coef_tbl_raw[[se_name]]) else NA_real_,
      statistic = if (!is.na(stat_name)) as.numeric(coef_tbl_raw[[stat_name]]) else NA_real_,
      p_value = if (!is.na(p_name)) as.numeric(coef_tbl_raw[[p_name]]) else NA_real_,
      # Fallback when p-value is not provided by model summary (e.g., plain lme4::lmer).
      p_value = dplyr::coalesce(
        p_value,
        2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
      )
    )

  parm_fix <- names(lme4::fixef(fit))
  ci_fix <- tryCatch(
    suppressMessages(confint(fit, parm = parm_fix, method = "Wald")),
    error = function(e) NULL
  )

  if (!is.null(ci_fix)) {
    ci_tbl <- as_tibble(as.data.frame(ci_fix), rownames = "term")
    names(ci_tbl)[2:3] <- c("ci_low", "ci_high")
    coef_tbl <- coef_tbl |>
      left_join(ci_tbl, by = "term")
  } else {
    coef_tbl <- coef_tbl |>
      mutate(ci_low = NA_real_, ci_high = NA_real_)
  }

  ig_row <- coef_tbl |>
    filter(term == "ig_gradient") |>
    slice_head(n = 1)

  r2_obj <- tryCatch(
    performance::r2_nakagawa(fit),
    error = function(e) NULL
  )

  r2_marginal <- if (!is.null(r2_obj) && "R2_marginal" %in% names(r2_obj)) unname(r2_obj$R2_marginal) else NA_real_
  r2_conditional <- if (!is.null(r2_obj) && "R2_conditional" %in% names(r2_obj)) unname(r2_obj$R2_conditional) else NA_real_

  summary_tbl <- tibble(
    outcome = outcome_i,
    boundary_type = boundary_type_i,
    n = nrow(dat),
    n_patients = dplyr::n_distinct(dat$patient_id),
    ig_beta = if (nrow(ig_row) == 1) ig_row$estimate else NA_real_,
    ig_se = if (nrow(ig_row) == 1) ig_row$std_error else NA_real_,
    ig_p_value = if (nrow(ig_row) == 1) ig_row$p_value else NA_real_,
    ig_ci_low = if (nrow(ig_row) == 1) ig_row$ci_low else NA_real_,
    ig_ci_high = if (nrow(ig_row) == 1) ig_row$ci_high else NA_real_,
    aic = AIC(fit),
    bic = BIC(fit),
    r2_marginal = r2_marginal,
    r2_conditional = r2_conditional,
    formula = paste(deparse(form), collapse = " "),
    covariates = if (length(covars_use) > 0) paste(covars_use, collapse = ", ") else NA_character_,
    note = "OK"
  )

  list(
    fit = fit,
    summary = summary_tbl,
    coef = coef_tbl,
    note = "OK"
  )
}

# **********************************************************************
### mmt8 ceiling cohort ----
# **********************************************************************
res_phvas_ceiling_mmt8_01 <- run_phvas_ceiling_model_01("mmt8")

# **********************************************************************
### mmt10 ceiling cohor ----
# **********************************************************************
res_phvas_ceiling_mmt10_01 <- run_phvas_ceiling_model_01("mmt10")

# **********************************************************************
### fi2 ceiling cohort ----
# **********************************************************************
res_phvas_ceiling_fi2_01 <- run_phvas_ceiling_model_01("fi2")

# **********************************************************************
### fi3 ceiling cohort ----
# **********************************************************************
res_phvas_ceiling_fi3_01 <- run_phvas_ceiling_model_01("fi3")

# **********************************************************************
### haq ceiling cohort ----
# **********************************************************************
res_phvas_ceiling_haq_01 <- run_phvas_ceiling_model_01("haq")

# **********************************************************************
### Collect and Export ----
# **********************************************************************
# Collect per-parameter model objects.
res_phvas_ceiling_models_01 <- list(
  mmt8 = res_phvas_ceiling_mmt8_01,
  mmt10 = res_phvas_ceiling_mmt10_01,
  fi2 = res_phvas_ceiling_fi2_01,
  fi3 = res_phvas_ceiling_fi3_01,
  haq = res_phvas_ceiling_haq_01
)

# Summary table across all selected clinical parameters.
res_phvas_ceiling_summary_01 <- map_dfr(res_phvas_ceiling_models_01, "summary")

# Fixed-effects coefficients across all selected clinical parameters.
res_phvas_ceiling_coef_01 <- imap_dfr(
  res_phvas_ceiling_models_01,
  function(x, nm) {
    if (nrow(x$coef) == 0) {
      return(
        tibble(
          outcome = nm,
          term = NA_character_,
          estimate = NA_real_,
          std_error = NA_real_,
          statistic = NA_real_,
          p_value = NA_real_,
          ci_low = NA_real_,
          ci_high = NA_real_,
          note = x$note
        )
      )
    }

    x$coef |>
      mutate(
        outcome = nm,
        note = x$note,
        .before = 1
      )
  }
)

# Status table for quick diagnostics.
res_phvas_ceiling_status_01 <- imap_dfr(
  res_phvas_ceiling_models_01,
  function(x, nm) {
    tibble(
      outcome = nm,
      has_fit = !is.null(x$fit),
      n_terms = nrow(x$coef),
      note = x$note
    )
  }
)

# Export PhVAS ceiling-only model results.
rio::export(
  list(
    phvas_ceiling_status = janitor::clean_names(res_phvas_ceiling_status_01),
    phvas_ceiling_summary = janitor::clean_names(res_phvas_ceiling_summary_01),
    phvas_ceiling_coefficients = janitor::clean_names(res_phvas_ceiling_coef_01)
  ),
  file.path(
    output_tab_dir_ceil_01,
    paste0(format(Sys.Date(), "%y%m%d"), "_phvas_ig_gradient_ceiling_models_01.xlsx")
  )
)


# **********************************************************************
# 6. Prediction perfomance - model comparison ----
# **********************************************************************

# **********************************************************************
## No ceiling/floor cohort - without any celing/floor  ----
# **********************************************************************
# Prediction performance compares models on the original PhVAS scale.
# Standardization is not applied because all compared models share the same response
# and we report absolute error metrics (MAE/RMSE) in clinically interpretable units.

# Build a common base dataset for prediction-performance analyses.
d_predperf_base_01 <- d_convert_model_01 |>
  left_join(d_anchor_activity_01, by = c("patient_id", "exam_order")) |>
  filter(
    !is.na(outcome),
    !is.na(anchor_activity),
    !is.na(ig_gradient),
    !is.na(patient_id)
  ) |>
  (\(dat) {
    if ("exam_order_factor" %in% names(dat)) {
      dat |>
        mutate(
          patient_id = factor(patient_id),
          exam_order_factor = normalize_exam_factor_01(exam_order_factor)
        )
    } else {
      dat |>
        mutate(
          patient_id = factor(patient_id),
          exam_order_factor = normalize_exam_factor_01(exam_order)
        )
    }
  })()

outcomes_predperf_01 <- sort(unique(d_predperf_base_01$outcome))

mode_value_01 <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[[which.max(tabulate(match(x, ux)))]]
}

predict_lmer_fixed_01 <- function(fit_obj, newdata) {
  tryCatch(
    as.numeric(stats::predict(fit_obj, newdata = newdata, re.form = NA, allow.new.levels = TRUE)),
    error = function(e) rep(NA_real_, nrow(newdata))
  )
}

calc_pred_metrics_vec_01 <- function(obs, pred) {
  dat <- tibble(observed = as.numeric(obs), predicted = as.numeric(pred)) |>
    filter(is.finite(observed), is.finite(predicted))

  if (nrow(dat) < 3) {
    return(
      tibble(
        n = nrow(dat),
        mae = NA_real_,
        rmse = NA_real_,
        r2 = NA_real_,
        pearson = NA_real_,
        spearman = NA_real_,
        bias = NA_real_,
        calib_intercept = NA_real_,
        calib_slope = NA_real_
      )
    )
  }

  resid <- dat$observed - dat$predicted
  sse <- sum(resid^2, na.rm = TRUE)
  sst <- sum((dat$observed - mean(dat$observed, na.rm = TRUE))^2, na.rm = TRUE)
  r2_val <- if (is.finite(sst) && sst > 0) 1 - sse / sst else NA_real_

  calib_fit <- tryCatch(
    stats::lm(observed ~ predicted, data = dat),
    error = function(e) NULL
  )

  calib_int <- if (!is.null(calib_fit)) unname(stats::coef(calib_fit)[["(Intercept)"]]) else NA_real_
  calib_slope <- if (!is.null(calib_fit)) unname(stats::coef(calib_fit)[["predicted"]]) else NA_real_

  tibble(
    n = nrow(dat),
    mae = mean(abs(resid), na.rm = TRUE),
    rmse = sqrt(mean(resid^2, na.rm = TRUE)),
    r2 = r2_val,
    pearson = suppressWarnings(stats::cor(dat$observed, dat$predicted, method = "pearson")),
    spearman = suppressWarnings(stats::cor(dat$observed, dat$predicted, method = "spearman")),
    bias = mean(resid, na.rm = TRUE),
    calib_intercept = calib_int,
    calib_slope = calib_slope
  )
}

fit_predperf_model_c_only_01 <- function(dat_outcome) {
  covars_pref <- c("age_m0_cov", "sex", "pohlavi", "bcm", "wear_time", "valid_days")
  covars_use <- get_covars_use_01(dat_outcome, covars_pref)

  exam_term <- if ("exam_order_factor" %in% names(dat_outcome) &&
    dplyr::n_distinct(dat_outcome$exam_order_factor) > 1) {
    "exam_order_factor"
  } else {
    NULL
  }

  rhs_c <- c("ig_gradient", exam_term, covars_use)
  form_c <- as.formula(paste0("anchor_activity ~ ", paste(c(rhs_c, "(1|patient_id)"), collapse = " + ")))
  fit_c <- fit_lmer_safe_01(form_c, dat_outcome)

  list(
    model_c = fit_c,
    formula_c = format(form_c),
    covars_use = covars_use
  )
}

fit_predperf_models_one_outcome_01 <- function(dat_outcome, mode = c("abc", "c_only")) {
  mode <- match.arg(mode)

  if (mode == "abc") {
    fit_obj <- fit_disease_models_one_outcome_01(dat_outcome)
    models_named <- make_models_named_abc_01(fit_obj)
  } else {
    fit_obj <- fit_predperf_model_c_only_01(dat_outcome)
    models_named <- setNames(list(fit_obj$model_c), unname(disease_model_labels_01[["model_c"]]))
  }

  models_named <- models_named[!map_lgl(models_named, is.null)]

  list(
    fit_obj = fit_obj,
    models_named = models_named
  )
}

build_predperf_pred_long_01 <- function(models_named, dat_outcome) {
  imap_dfr(
    models_named,
    function(fit_obj, model_label) {
      if (is.null(fit_obj)) return(tibble())

      pred <- predict_lmer_fixed_01(fit_obj, dat_outcome)
      tibble(
        model_label = as.character(model_label),
        observed = as.numeric(dat_outcome$anchor_activity),
        predicted = as.numeric(pred),
        ig_gradient = as.numeric(dat_outcome$ig_gradient),
        exam_order_factor = if ("exam_order_factor" %in% names(dat_outcome)) as.character(dat_outcome$exam_order_factor) else NA_character_,
        patient_id = if ("patient_id" %in% names(dat_outcome)) as.character(dat_outcome$patient_id) else NA_character_
      )
    }
  ) |>
    filter(is.finite(observed), is.finite(predicted))
}

build_predperf_accuracy_plot_01 <- function(pred_long, outcome_nm, cohort_subtitle) {
  if (nrow(pred_long) == 0) return(NULL)

  ggplot(pred_long, aes(x = observed, y = predicted, color = model_label, fill = model_label)) +
    geom_abline(linetype = "dashed", color = "grey40") +
    geom_point(alpha = 0.35, size = 1.8) +
    geom_smooth(method = "lm", alpha = 0.12, linewidth = 1.1) +
    labs(
      title = paste0("Outcome: ", outcome_nm),
      subtitle = cohort_subtitle,
      x = "PhVAS (observed)",
      y = "PhVAS (predicted)",
      color = NULL,
      fill = NULL
    ) +
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(size = 24, face = "bold"),
      plot.subtitle = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 16, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 15, face = "bold")
    )
}

# Bootstrap 95% CI for RMSE and MAE from observed-predicted pairs.
bootstrap_mae_rmse_ci_01 <- function(pred_df, n_boot = 500, conf_level = 0.95, seed = NULL) {
  dat <- pred_df |>
    transmute(
      observed = as.numeric(observed),
      predicted = as.numeric(predicted),
      patient_id = as.character(patient_id)
    ) |>
    filter(is.finite(observed), is.finite(predicted))

  if (nrow(dat) < 5) {
    return(
      tibble(
        metric = c("RMSE", "MAE"),
        estimate = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        n = nrow(dat),
        n_boot = n_boot,
        bootstrap_unit = NA_character_
      )
    )
  }

  metric_fun <- function(obs, pred) {
    tibble(
      metric = c("RMSE", "MAE"),
      value = c(
        sqrt(mean((obs - pred)^2, na.rm = TRUE)),
        mean(abs(obs - pred), na.rm = TRUE)
      )
    )
  }

  if (!is.null(seed)) set.seed(seed)
  alpha <- (1 - conf_level) / 2

  has_patient_boot <- "patient_id" %in% names(dat) &&
    any(!is.na(dat$patient_id)) &&
    dplyr::n_distinct(dat$patient_id[!is.na(dat$patient_id)]) >= 3

  boot_tbl <- if (has_patient_boot) {
    ids <- unique(dat$patient_id[!is.na(dat$patient_id)])
    map_dfr(
      seq_len(n_boot),
      \(b) {
        ids_b <- sample(ids, size = length(ids), replace = TRUE)
        dat_b <- map_dfr(ids_b, \(id_i) dat |> filter(patient_id == id_i))
        metric_fun(dat_b$observed, dat_b$predicted) |>
          mutate(iter = b)
      }
    )
  } else {
    n_i <- nrow(dat)
    map_dfr(
      seq_len(n_boot),
      \(b) {
        idx <- sample.int(n_i, size = n_i, replace = TRUE)
        dat_b <- dat[idx, , drop = FALSE]
        metric_fun(dat_b$observed, dat_b$predicted) |>
          mutate(iter = b)
      }
    )
  }

  ci_tbl <- boot_tbl |>
    group_by(metric) |>
    summarise(
      ci_low = stats::quantile(value, probs = alpha, na.rm = TRUE),
      ci_high = stats::quantile(value, probs = 1 - alpha, na.rm = TRUE),
      .groups = "drop"
    )

  metric_fun(dat$observed, dat$predicted) |>
    rename(estimate = value) |>
    left_join(ci_tbl, by = "metric") |>
    mutate(
      n = nrow(dat),
      n_boot = n_boot,
      bootstrap_unit = if_else(has_patient_boot, "patient_id", "row")
    )
}

# One outcome-level uncertainty plot: models compared for RMSE and MAE with bootstrap 95% CI.
build_predperf_uncertainty_plot_01 <- function(metric_ci_tbl, outcome_nm, cohort_subtitle) {
  if (nrow(metric_ci_tbl) == 0) return(NULL)

  model_levels <- unname(disease_model_labels_01[c("model_a", "model_b", "model_c")])

  plot_dat <- metric_ci_tbl |>
    mutate(
      model_label = factor(model_label, levels = model_levels),
      metric = factor(metric, levels = c("RMSE", "MAE"))
    ) |>
    filter(!is.na(model_label), !is.na(metric))

  if (nrow(plot_dat) == 0) return(NULL)

  ggplot(plot_dat, aes(x = model_label, y = estimate, color = model_label)) +
    geom_pointrange(
      aes(ymin = ci_low, ymax = ci_high),
      linewidth = 0.8,
      fatten = 2.0
    ) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(
      title = paste0("Outcome: ", outcome_nm),
      subtitle = cohort_subtitle,
      x = NULL,
      y = "Metric value\n(bootstrap 95% CI)",
      color = NULL
    ) +
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(size = 24, face = "bold"),
      plot.subtitle = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(angle = 12, hjust = 1),
      strip.text = element_text(size = 18, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 15, face = "bold")
    )
}

make_grid_for_fit_01 <- function(fit_obj, dat_outcome, n_points = 80) {
  if (is.null(fit_obj) || nrow(dat_outcome) == 0) return(tibble())

  ig_rng <- range(dat_outcome$ig_gradient, na.rm = TRUE)
  if (!all(is.finite(ig_rng))) return(tibble())

  has_exam <- "exam_order_factor" %in% names(dat_outcome) &&
    dplyr::n_distinct(dat_outcome$exam_order_factor) > 1

  base_grid <- if (has_exam) {
    expand_grid(
      ig_gradient = seq(ig_rng[1], ig_rng[2], length.out = n_points),
      exam_order_factor = factor(levels(droplevels(dat_outcome$exam_order_factor)), levels = levels(droplevels(dat_outcome$exam_order_factor)))
    )
  } else {
    tibble(
      ig_gradient = seq(ig_rng[1], ig_rng[2], length.out = n_points),
      exam_order_factor = if ("exam_order_factor" %in% names(dat_outcome)) mode_value_01(dat_outcome$exam_order_factor) else NA
    )
  }

  if ("y" %in% names(dat_outcome)) {
    base_grid <- base_grid |>
      mutate(y = median(dat_outcome$y, na.rm = TRUE))
  }

  covar_candidates <- c("age_m0_cov", "sex", "pohlavi", "bcm", "wear_time", "valid_days")
  covars_present <- intersect(covar_candidates, names(dat_outcome))

  for (nm in covars_present) {
    if (is.numeric(dat_outcome[[nm]])) {
      base_grid[[nm]] <- median(dat_outcome[[nm]], na.rm = TRUE)
    } else {
      mv <- mode_value_01(dat_outcome[[nm]])
      base_grid[[nm]] <- if (is.factor(dat_outcome[[nm]])) {
        factor(mv, levels = levels(dat_outcome[[nm]]))
      } else {
        mv
      }
    }
  }

  if ("patient_id" %in% names(dat_outcome)) {
    pid_ref <- levels(factor(dat_outcome$patient_id))
    pid_ref <- if (length(pid_ref) > 0) pid_ref[[1]] else as.character(dat_outcome$patient_id[[1]])
    base_grid <- base_grid |>
      mutate(patient_id = factor(pid_ref, levels = levels(factor(dat_outcome$patient_id))))
  }

  pred <- predict_lmer_fixed_01(fit_obj, base_grid)
  base_grid |>
    mutate(predicted = pred) |>
    filter(is.finite(predicted))
}

build_predperf_grid_plot_01 <- function(models_named, dat_outcome, outcome_nm, cohort_subtitle) {
  if (length(models_named) == 0 || nrow(dat_outcome) == 0) return(NULL)

  grid_long <- imap_dfr(
    models_named,
    function(fit_obj, model_label) {
      make_grid_for_fit_01(fit_obj, dat_outcome) |>
        mutate(model_label = as.character(model_label))
    }
  )

  if (nrow(grid_long) == 0) return(NULL)

  raw_dat <- dat_outcome |>
    transmute(
      ig_gradient = as.numeric(ig_gradient),
      anchor_activity = as.numeric(anchor_activity),
      exam_order_factor = as.character(exam_order_factor)
    ) |>
    filter(is.finite(ig_gradient), is.finite(anchor_activity))

  p <- ggplot() +
    geom_point(
      data = raw_dat,
      aes(x = ig_gradient, y = anchor_activity),
      color = "grey45",
      alpha = 0.30,
      size = 1.8
    ) +
    geom_line(
      data = grid_long,
      aes(x = ig_gradient, y = predicted, color = model_label),
      linewidth = 1.2
    ) +
    labs(
      title = paste0("Outcome: ", outcome_nm),
      subtitle = cohort_subtitle,
      x = x_lab_01,
      y = "PhVAS (predicted)",
      color = NULL
    ) +
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(size = 24, face = "bold"),
      plot.subtitle = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 16, face = "bold"),
      strip.text = element_text(size = 16, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 15, face = "bold")
    )

  if ("exam_order_factor" %in% names(grid_long) && dplyr::n_distinct(grid_long$exam_order_factor) > 1) {
    p <- p + facet_wrap(~ exam_order_factor)
  }

  p
}

combine_plot_list_01 <- function(plot_list, title_txt, subtitle_txt) {
  if (length(plot_list) == 0 || !requireNamespace("patchwork", quietly = TRUE)) return(NULL)

  patchwork::wrap_plots(plot_list, ncol = 2) +
    patchwork::plot_annotation(
      title = title_txt,
      subtitle = subtitle_txt,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 30, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 22, face = "bold")
      )
    )
}

run_predperf_cohort_01 <- function(data_in, cohort_label, fit_mode = c("abc", "c_only")) {
  fit_mode <- match.arg(fit_mode)
  set.seed(20260305 + sum(as.integer(charToRaw(cohort_label))))

  res_by_outcome <- outcomes_predperf_01 |>
    set_names() |>
    map(function(outcome_i) {
      dat <- data_in |>
        filter(outcome == outcome_i) |>
        filter(!is.na(anchor_activity), !is.na(ig_gradient), !is.na(patient_id))

      if (fit_mode == "abc") {
        dat <- dat |>
          filter(!is.na(y))
      }

      if (nrow(dat) < 20 || dplyr::n_distinct(dat$patient_id) < 5) {
        return(
          list(
            outcome = outcome_i,
            n = nrow(dat),
            n_patients = dplyr::n_distinct(dat$patient_id),
            models_named = list(),
            pred_long = tibble(),
            pred_metrics = tibble(),
            metric_ci_boot = tibble(),
            compare = NULL,
            accuracy_plot = NULL,
            grid_plot = NULL,
            uncertainty_plot = NULL,
            note = "Insufficient data for mixed-model prediction performance."
          )
        )
      }

      fit_res <- fit_predperf_models_one_outcome_01(dat, mode = fit_mode)
      models_named <- fit_res$models_named
      pred_long <- build_predperf_pred_long_01(models_named, dat)

      pred_metrics <- if (nrow(pred_long) > 0) {
        pred_long |>
          group_by(model_label) |>
          group_modify(~ calc_pred_metrics_vec_01(.x$observed, .x$predicted)) |>
          ungroup() |>
          mutate(
            outcome = outcome_i,
            cohort = cohort_label,
            n_outcome = nrow(dat),
            n_patients = dplyr::n_distinct(dat$patient_id),
            .before = 1
          )
      } else {
        tibble()
      }

      metric_ci_boot <- if (nrow(pred_long) > 0) {
        pred_long |>
          group_by(model_label) |>
          group_modify(~ bootstrap_mae_rmse_ci_01(.x, n_boot = 500, conf_level = 0.95)) |>
          ungroup() |>
          mutate(
            outcome = outcome_i,
            cohort = cohort_label,
            n_outcome = nrow(dat),
            n_patients = dplyr::n_distinct(dat$patient_id),
            .before = 1
          )
      } else {
        tibble()
      }

      compare_tbl <- NULL
      if (length(models_named) >= 2 && requireNamespace("performance", quietly = TRUE)) {
        compare_tbl <- do.call(performance::compare_performance, models_named) |>
          relabel_compare_performance_01(model_labels = names(models_named)) |>
          as_tibble() |>
          mutate(
            outcome = outcome_i,
            cohort = cohort_label,
            .before = 1
          )
      }

      list(
        outcome = outcome_i,
        n = nrow(dat),
        n_patients = dplyr::n_distinct(dat$patient_id),
        models_named = models_named,
        pred_long = pred_long,
        pred_metrics = pred_metrics,
        metric_ci_boot = metric_ci_boot,
        compare = compare_tbl,
        accuracy_plot = build_predperf_accuracy_plot_01(pred_long, outcome_i, cohort_label),
        grid_plot = build_predperf_grid_plot_01(models_named, dat, outcome_i, cohort_label),
        uncertainty_plot = build_predperf_uncertainty_plot_01(metric_ci_boot, outcome_i, cohort_label),
        note = if (length(models_named) == 0) "No valid fitted model." else "OK"
      )
    })

  status_tbl <- imap_dfr(
    res_by_outcome,
    function(x, nm) {
      tibble(
        outcome = nm,
        cohort = cohort_label,
        n = x$n,
        n_patients = x$n_patients,
        n_models = length(x$models_named),
        model_labels = paste(names(x$models_named), collapse = " | "),
        has_compare = !is.null(x$compare),
        note = x$note
      )
    }
  )

  metrics_tbl <- map_dfr(res_by_outcome, "pred_metrics")
  metric_ci_tbl <- map_dfr(res_by_outcome, "metric_ci_boot")
  compare_tbl <- map_dfr(res_by_outcome, "compare")
  accuracy_plot_list <- map(res_by_outcome, "accuracy_plot") |> purrr::compact()
  grid_plot_list <- map(res_by_outcome, "grid_plot") |> purrr::compact()
  uncertainty_plot_list <- map(res_by_outcome, "uncertainty_plot") |> purrr::compact()

  list(
    by_outcome = res_by_outcome,
    status = status_tbl,
    pred_metrics = metrics_tbl,
    metric_ci_boot = metric_ci_tbl,
    compare = compare_tbl,
    accuracy_plot_list = accuracy_plot_list,
    grid_plot_list = grid_plot_list,
    uncertainty_plot_list = uncertainty_plot_list
  )
}

# Keep only patients who never reached outcome-specific ceiling/floor.
d_predperf_no_cei_flo_patients_01 <- d_predperf_base_01 |>
  group_by(outcome, patient_id) |>
  mutate(patient_has_boundary = any(boundary_flag == 1L, na.rm = TRUE)) |>
  ungroup() |>
  filter(!patient_has_boundary) |>
  select(-patient_has_boundary)

res_predperf_no_cei_flo_01 <- run_predperf_cohort_01(
  data_in = d_predperf_no_cei_flo_patients_01,
  cohort_label = "No ceiling/floor cohort - patients without any boundary value",
  fit_mode = "abc"
)

# **********************************************************************
## All cohort ----
# **********************************************************************
# Full cohort (all available rows).
d_predperf_all_01 <- d_predperf_base_01

res_predperf_all_01 <- run_predperf_cohort_01(
  data_in = d_predperf_all_01,
  cohort_label = "All cohort",
  fit_mode = "abc"
)


# **********************************************************************
## Ceiling/floor cohort: ceiling-only  ----
# **********************************************************************
# Keep only patients whose observed values are boundary-only across visits.
d_predperf_ceiling_only_01 <- d_predperf_base_01 |>
  group_by(outcome, patient_id) |>
  mutate(
    n_obs_patient = sum(!is.na(y)),
    n_boundary_patient = sum(boundary_flag == 1L, na.rm = TRUE),
    patient_all_boundary = n_obs_patient > 0 & n_obs_patient == n_boundary_patient
  ) |>
  ungroup() |>
  filter(patient_all_boundary) |>
  select(-n_obs_patient, -n_boundary_patient, -patient_all_boundary)

res_predperf_ceiling_only_01 <- run_predperf_cohort_01(
  data_in = d_predperf_ceiling_only_01,
  cohort_label = "Ceiling/floor cohort: boundary-only patients",
  fit_mode = "c_only"
)

# **********************************************************************
## Ceiling/floor occurence in cohort  ----
# **********************************************************************
# Keep patients who reached outcome-specific boundary at least once.
d_predperf_boundary_occurrence_01 <- d_predperf_base_01 |>
  group_by(outcome, patient_id) |>
  mutate(patient_has_boundary = any(boundary_flag == 1L, na.rm = TRUE)) |>
  ungroup() |>
  filter(patient_has_boundary) |>
  select(-patient_has_boundary)

res_predperf_boundary_occurrence_01 <- run_predperf_cohort_01(
  data_in = d_predperf_boundary_occurrence_01,
  cohort_label = "Ceiling/floor occurrence cohort (at least one boundary visit)",
  fit_mode = "abc"
)

# **********************************************************************
### Figures ----
# **********************************************************************
# Build two multi-figures per cohort:
# 1) observed vs predicted comparison across models
# 2) modelr-like evaluation on an explanatory ig_gradient grid
res_predperf_accuracy_multi_no_cei_flo_01 <- combine_plot_list_01(
  res_predperf_no_cei_flo_01$accuracy_plot_list,
  title_txt = "Prediction accuracy: observed vs predicted",
  subtitle_txt = "No ceiling/floor cohort - patients without any boundary value"
)
res_predperf_grid_multi_no_cei_flo_01 <- combine_plot_list_01(
  res_predperf_no_cei_flo_01$grid_plot_list,
  title_txt = paste0("Prediction evaluation on ", x_lab_01, " grid"),
  subtitle_txt = "No ceiling/floor cohort - patients without any boundary value"
)

res_predperf_accuracy_multi_all_01 <- combine_plot_list_01(
  res_predperf_all_01$accuracy_plot_list,
  title_txt = "Prediction accuracy: observed vs predicted",
  subtitle_txt = "All cohort"
)
res_predperf_grid_multi_all_01 <- combine_plot_list_01(
  res_predperf_all_01$grid_plot_list,
  title_txt = paste0("Prediction evaluation on ", x_lab_01, " grid"),
  subtitle_txt = "All cohort"
)

res_predperf_accuracy_multi_ceiling_only_01 <- combine_plot_list_01(
  res_predperf_ceiling_only_01$accuracy_plot_list,
  title_txt = "Prediction accuracy: observed vs predicted",
  subtitle_txt = "Ceiling/floor cohort: boundary-only patients (Model 3 only)"
)
res_predperf_grid_multi_ceiling_only_01 <- combine_plot_list_01(
  res_predperf_ceiling_only_01$grid_plot_list,
  title_txt = paste0("Prediction evaluation on ", x_lab_01, " grid"),
  subtitle_txt = "Ceiling/floor cohort: boundary-only patients (Model 3 only)"
)

res_predperf_accuracy_multi_boundary_occ_01 <- combine_plot_list_01(
  res_predperf_boundary_occurrence_01$accuracy_plot_list,
  title_txt = "Prediction accuracy: observed vs predicted",
  subtitle_txt = "Ceiling/floor occurrence cohort (at least one boundary visit)"
)
res_predperf_grid_multi_boundary_occ_01 <- combine_plot_list_01(
  res_predperf_boundary_occurrence_01$grid_plot_list,
  title_txt = paste0("Prediction evaluation on ", x_lab_01, " grid"),
  subtitle_txt = "Ceiling/floor occurrence cohort (at least one boundary visit)"
)

# 3) Metric comparison with uncertainty (bootstrap 95% CI for RMSE and MAE).
res_predperf_uncertainty_multi_no_cei_flo_01 <- combine_plot_list_01(
  res_predperf_no_cei_flo_01$uncertainty_plot_list,
  title_txt = "Metric comparison with uncertainty (bootstrap 95% CI)",
  subtitle_txt = "No ceiling/floor cohort - patients without any boundary value"
)
res_predperf_uncertainty_multi_all_01 <- combine_plot_list_01(
  res_predperf_all_01$uncertainty_plot_list,
  title_txt = "Metric comparison with uncertainty (bootstrap 95% CI)",
  subtitle_txt = "All cohort"
)
res_predperf_uncertainty_multi_ceiling_only_01 <- combine_plot_list_01(
  res_predperf_ceiling_only_01$uncertainty_plot_list,
  title_txt = "Metric comparison with uncertainty (bootstrap 95% CI)",
  subtitle_txt = "Ceiling/floor cohort: boundary-only patients (Model 3 only)"
)
res_predperf_uncertainty_multi_boundary_occ_01 <- combine_plot_list_01(
  res_predperf_boundary_occurrence_01$uncertainty_plot_list,
  title_txt = "Metric comparison with uncertainty (bootstrap 95% CI)",
  subtitle_txt = "Ceiling/floor occurrence cohort (at least one boundary visit)"
)

save_predperf_plot_01 <- function(plot_obj, file_name, width = 16, height = 11) {
  fp <- file.path(output_fig_dir_ceil_01, file_name)
  ok <- if (!is.null(plot_obj)) {
    save_plot_safe_01(
      plot_obj = plot_obj,
      file_path = fp,
      width = width,
      height = height,
      dpi = 300
    )
  } else {
    FALSE
  }
  tibble(file_path = fp, saved = ok)
}

run_date_predperf_01 <- format(Sys.Date(), "%y%m%d")

res_predperf_fig_manifest_01 <- bind_rows(
  save_predperf_plot_01(res_predperf_accuracy_multi_no_cei_flo_01, paste0(run_date_predperf_01, "_predperf_accuracy_no_cei_flo_no_boundary_01.png")),
  save_predperf_plot_01(res_predperf_grid_multi_no_cei_flo_01, paste0(run_date_predperf_01, "_predperf_grid_no_cei_flo_no_boundary_01.png")),
  save_predperf_plot_01(res_predperf_uncertainty_multi_no_cei_flo_01, paste0(run_date_predperf_01, "_predperf_uncertainty_no_cei_flo_no_boundary_01.png")),
  save_predperf_plot_01(res_predperf_accuracy_multi_all_01, paste0(run_date_predperf_01, "_predperf_accuracy_all_cohort_01.png")),
  save_predperf_plot_01(res_predperf_grid_multi_all_01, paste0(run_date_predperf_01, "_predperf_grid_all_cohort_01.png")),
  save_predperf_plot_01(res_predperf_uncertainty_multi_all_01, paste0(run_date_predperf_01, "_predperf_uncertainty_all_cohort_01.png")),
  save_predperf_plot_01(res_predperf_accuracy_multi_ceiling_only_01, paste0(run_date_predperf_01, "_predperf_accuracy_ceiling_only_01.png")),
  save_predperf_plot_01(res_predperf_grid_multi_ceiling_only_01, paste0(run_date_predperf_01, "_predperf_grid_ceiling_only_01.png")),
  save_predperf_plot_01(res_predperf_uncertainty_multi_ceiling_only_01, paste0(run_date_predperf_01, "_predperf_uncertainty_ceiling_only_01.png")),
  save_predperf_plot_01(res_predperf_accuracy_multi_boundary_occ_01, paste0(run_date_predperf_01, "_predperf_accuracy_boundary_occurrence_01.png")),
  save_predperf_plot_01(res_predperf_grid_multi_boundary_occ_01, paste0(run_date_predperf_01, "_predperf_grid_boundary_occurrence_01.png")),
  save_predperf_plot_01(res_predperf_uncertainty_multi_boundary_occ_01, paste0(run_date_predperf_01, "_predperf_uncertainty_boundary_occurrence_01.png"))
)

# **********************************************************************
### Tables & export ----
# **********************************************************************
res_predperf_tab_export_01 <- list(
  cohort1_status = janitor::clean_names(res_predperf_no_cei_flo_01$status),
  cohort1_pred_metrics = janitor::clean_names(res_predperf_no_cei_flo_01$pred_metrics),
  cohort1_metric_ci_boot = janitor::clean_names(res_predperf_no_cei_flo_01$metric_ci_boot),
  cohort1_model_compare = janitor::clean_names(res_predperf_no_cei_flo_01$compare),

  cohort2_status = janitor::clean_names(res_predperf_all_01$status),
  cohort2_pred_metrics = janitor::clean_names(res_predperf_all_01$pred_metrics),
  cohort2_metric_ci_boot = janitor::clean_names(res_predperf_all_01$metric_ci_boot),
  cohort2_model_compare = janitor::clean_names(res_predperf_all_01$compare),

  cohort3_status = janitor::clean_names(res_predperf_ceiling_only_01$status),
  cohort3_pred_metrics = janitor::clean_names(res_predperf_ceiling_only_01$pred_metrics),
  cohort3_metric_ci_boot = janitor::clean_names(res_predperf_ceiling_only_01$metric_ci_boot),

  cohort4_status = janitor::clean_names(res_predperf_boundary_occurrence_01$status),
  cohort4_pred_metrics = janitor::clean_names(res_predperf_boundary_occurrence_01$pred_metrics),
  cohort4_metric_ci_boot = janitor::clean_names(res_predperf_boundary_occurrence_01$metric_ci_boot),
  cohort4_model_compare = janitor::clean_names(res_predperf_boundary_occurrence_01$compare),

  figure_manifest = janitor::clean_names(res_predperf_fig_manifest_01)
)

rio::export(
  res_predperf_tab_export_01,
  file.path(
    output_tab_dir_ceil_01,
    paste0(run_date_predperf_01, "_prediction_performance_model_comparison_01.xlsx")
  )
)




