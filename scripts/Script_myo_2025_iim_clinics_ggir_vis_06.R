# ************************************************************
# TIS vs GGIR explanatory value for SF-36 PCS
# For each GGIR parameter compare:
# 1) TIS only
# 2) GGIR only
# 3) TIS + GGIR
# Outcome: log(sf36_pcs)
# Covariates: exam + age_m0 + sex + sf36_mcs
# ************************************************************

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  dplyr,
  purrr,
  tibble,
  stringr,
  tidyr,
  ggplot2,
  ggbeeswarm,
  patchwork,
  future,
  furrr,
  easystats,
  rio,
  broom.mixed,
  lme4
)

safe_future_workers_vis06 <- as.integer(getOption(
  "myo_ncores",
  max(1L, future::availableCores() - 1L)
))

if (!is.finite(safe_future_workers_vis06) || is.na(safe_future_workers_vis06)) {
  safe_future_workers_vis06 <- max(1L, future::availableCores() - 1L)
}

safe_future_workers_vis06 <- max(1L, safe_future_workers_vis06)

fixed_model_colors_vis06 <- c(
  tis_only = "#F8766D",
  ggir_only = "#00BA38",
  tis_ggir = "#619CFF"
)

scale01_vis06 <- function(x, higher_better = TRUE) {
  x <- as.numeric(x)

  if (!higher_better) {
    x <- -x
  }

  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)

  if (!any(ok)) {
    return(out)
  }

  rng <- range(x[ok], na.rm = TRUE)

  if (diff(rng) == 0) {
    out[ok] <- 0.5
  } else {
    out[ok] <- (x[ok] - rng[1]) / diff(rng)
  }

  out
}

resolve_first_existing_vis06 <- function(data, candidates, label) {
  hit <- intersect(candidates, names(data))

  if (length(hit) == 0) {
    stop(
      "V datech chybi pozadovany sloupec pro `", label, "`. ",
      "Zkouseno: ", paste(candidates, collapse = ", "), "."
    )
  }

  hit[1]
}

sanitize_name_vis06 <- function(x) {
  stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
}

merge_tis_into_ggir_vis06 <- function(ggir_data, tis_data) {
  tis_join <- tis_data |>
    select(any_of(c("immet_id", "exam", "tis_total"))) |>
    mutate(
      immet_id = as.character(immet_id),
      exam = as.character(exam),
      tis_total = as.numeric(tis_total)
    )

  ggir_data |>
    mutate(
      immet_id = as.character(.data[["immet_id"]]),
      exam = as.character(.data[["exam"]])
    ) |>
    left_join(tis_join, by = c("immet_id", "exam"))
}

prepare_tis_ggir_eff_data_vis06 <- function(ggir_data, tis_data) {
  merged_data <- merge_tis_into_ggir_vis06(ggir_data, tis_data)

  immet_col <- resolve_first_existing_vis06(merged_data, c("immet_id"), "immet_id")
  exam_col <- resolve_first_existing_vis06(merged_data, c("exam_order", "exam"), "exam")
  sex_col <- resolve_first_existing_vis06(merged_data, c("sex", "pohlavi", "gender"), "sex")
  age_col <- resolve_first_existing_vis06(merged_data, c("age_m0", "age"), "age_m0")
  sf36_pcs_col <- resolve_first_existing_vis06(merged_data, c("sf36_pcs"), "sf36_pcs")
  sf36_mcs_col <- resolve_first_existing_vis06(merged_data, c("sf36_mcs"), "sf36_mcs")
  tis_col <- resolve_first_existing_vis06(merged_data, c("tis_total"), "tis_total")

  ggir_candidates <- c(
    "ig_gradient_enmo",
    "ig_intercept_enmo",
    "acc_day_spt_wei",
    "m5_enmo",
    "mvpa_t100_time_min",
    "p9931_enmo"
  )
  ggir_cols <- intersect(ggir_candidates, names(merged_data))

  merged_data |>
    transmute(
      immet_id = as.character(.data[[immet_col]]),
      exam = as.character(.data[[exam_col]]),
      sex = as.factor(.data[[sex_col]]),
      age_m0 = as.numeric(.data[[age_col]]),
      sf36_pcs = as.numeric(.data[[sf36_pcs_col]]),
      sf36_mcs = as.numeric(.data[[sf36_mcs_col]]),
      tis_total = as.numeric(.data[[tis_col]]),
      across(all_of(ggir_cols), as.numeric)
    ) |>
    mutate(
      exam = factor(
        exam,
        levels = unique(exam)[order(as.numeric(sub("^M", "", unique(exam))))]
      )
    )
}

fit_lmer_safe_vis06 <- function(formula_obj, data) {
  tryCatch(
    lme4::lmer(formula_obj, data = data, REML = FALSE),
    error = function(e) NULL
  )
}

relabel_compare_performance_vis06 <- function(perf_obj, model_labels) {
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
    perf_obj[[cc]] <- dplyr::recode(perf_obj[[cc]], !!!stats::setNames(model_labels, paste0("Model ", seq_along(model_labels))))
  }

  perf_obj
}

build_models_for_parameter_vis06 <- function(data, ggir_parameter, log_outcome = FALSE) {
  dat <- data |>
    select(immet_id, exam, age_m0, sex, sf36_pcs, sf36_mcs, tis_total, all_of(ggir_parameter)) |>
    na.omit()

  if (isTRUE(log_outcome)) {
    dat <- dat |>
      filter(sf36_pcs > 0)
  }

  outcome_expr <- if (isTRUE(log_outcome)) "log(sf36_pcs)" else "sf36_pcs"

  formulas <- list(
    tis_only = stats::as.formula(
      paste0(
        outcome_expr,
        " ~ tis_total + exam + age_m0 + sex + sf36_mcs + (1 | immet_id)"
      )
    ),
    ggir_only = stats::as.formula(
      paste0(
        outcome_expr,
        " ~ ", ggir_parameter, " + exam + age_m0 + sex + sf36_mcs + (1 | immet_id)"
      )
    ),
    tis_ggir = stats::as.formula(
      paste0(
        outcome_expr,
        " ~ tis_total + ", ggir_parameter, " + exam + age_m0 + sex + sf36_mcs + (1 | immet_id)"
      )
    )
  )

  models <- purrr::map(formulas, fit_lmer_safe_vis06, data = dat)

  list(
    data = dat,
    formulas = formulas,
    models = models
  )
}

compare_models_vis06 <- function(models_named, ggir_parameter) {
  valid_models <- models_named[!purrr::map_lgl(models_named, is.null)]

  if (length(valid_models) < 2) {
    return(NULL)
  }

  cmp_obj <- tryCatch(
    do.call(
      performance::compare_performance,
      c(valid_models, list(rank = TRUE, verbose = FALSE))
    ),
    error = function(e) NULL
  )

  if (is.null(cmp_obj)) {
    return(NULL)
  }

  relabel_compare_performance_vis06(cmp_obj, model_labels = names(valid_models)) |>
    as_tibble() |>
    mutate(ggir_parameter = ggir_parameter, .before = 1)
}

tidy_models_vis06 <- function(models_named, ggir_parameter) {
  purrr::imap_dfr(models_named, function(mod, model_name) {
    if (is.null(mod)) {
      return(tibble())
    }

    broom.mixed::tidy(
      mod,
      effects = "fixed",
      conf.int = TRUE,
      conf.method = "Wald"
    ) |>
      mutate(
        ggir_parameter = ggir_parameter,
        model = model_name,
        .before = 1
      )
  })
}

build_compare_plot_vis06 <- function(models_named, ggir_parameter, log_outcome = TRUE) {
  valid_models <- models_named[!purrr::map_lgl(models_named, is.null)]

  if (length(valid_models) < 2) {
    return(NULL)
  }

  cmp_obj <- tryCatch(
    do.call(
      performance::compare_performance,
      c(valid_models, list(rank = TRUE, verbose = FALSE))
    ),
    error = function(e) NULL
  )

  if (is.null(cmp_obj)) {
    return(NULL)
  }

  cmp_obj <- relabel_compare_performance_vis06(cmp_obj, model_labels = names(valid_models))

  p_obj <- tryCatch(
    plot(cmp_obj),
    error = function(e) NULL
  )

  if (is.null(p_obj)) {
    return(NULL)
  }

  if (inherits(p_obj, "ggplot")) {
    outcome_label <- if (isTRUE(log_outcome)) "log(sf36_pcs)" else "sf36_pcs"

    p_obj <- p_obj +
      ggplot2::labs(
        title = paste0("compare_performance: ", ggir_parameter),
        subtitle = paste0(
          outcome_label, " ~ tis_total + ", ggir_parameter,
          " + exam + age_m0 + sex + sf36_mcs + (1 | immet_id)"
        )
      ) +
      ggplot2::theme_minimal(base_size = 28) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 40, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 35, face = "italic", colour = "black"),
        axis.title = ggplot2::element_text(size = 28, face = "bold"),
        axis.text = ggplot2::element_text(size = 28, face = "bold"),
        legend.title = ggplot2::element_text(size = 28, face = "bold"),
        legend.text = ggplot2::element_text(size = 28, face = "bold")
      )
  }

  p_obj
}

build_prediction_df_vis06 <- function(models_named, data, log_outcome = TRUE) {
  obs_vec <- if (isTRUE(log_outcome)) log(data$sf36_pcs) else data$sf36_pcs

  purrr::imap_dfr(models_named, function(mod, model_name) {
    if (is.null(mod)) {
      return(tibble())
    }

    fitted_vals <- stats::predict(mod)
    resid_vals <- stats::residuals(mod)

    tibble(
      model = factor(model_name, levels = names(models_named)),
      observed = as.numeric(obs_vec),
      predicted = as.numeric(fitted_vals),
      residual = as.numeric(resid_vals),
      fitted = as.numeric(fitted_vals)
    )
  })
}

compute_metrics_vis06 <- function(observed, predicted) {
  tibble(
    rmse = sqrt(mean((observed - predicted)^2, na.rm = TRUE)),
    mae = mean(abs(observed - predicted), na.rm = TRUE),
    r2 = suppressWarnings(stats::cor(observed, predicted, use = "complete.obs")^2)
  )
}

build_cluster_bootstrap_data_vis06 <- function(data, id_col = "immet_id") {
  ids <- unique(as.character(data[[id_col]]))
  sampled_ids <- sample(ids, size = length(ids), replace = TRUE)

  purrr::map2_dfr(sampled_ids, seq_along(sampled_ids), function(id_i, draw_i) {
    data |>
      filter(.data[[id_col]] == id_i) |>
      mutate(immet_id = paste0(id_i, "_boot_", draw_i))
  })
}

bootstrap_model_metrics_vis06 <- function(
    formulas,
    data,
    log_outcome = TRUE,
    n_boot = 500,
    parallel = TRUE,
    workers = safe_future_workers_vis06,
    seed = 123
) {
  workers <- as.integer(workers)

  if (!is.finite(workers) || is.na(workers)) {
    workers <- safe_future_workers_vis06
  }

  workers <- max(1L, workers)

  if (isTRUE(parallel) && workers > 1L) {
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan("sequential"), add = TRUE)

    furrr::future_map_dfr(
      seq_len(n_boot),
      function(iter_i) {
        boot_data <- build_cluster_bootstrap_data_vis06(data)
        obs_vec <- if (isTRUE(log_outcome)) log(boot_data$sf36_pcs) else boot_data$sf36_pcs

        purrr::imap_dfr(formulas, function(formula_obj, model_name) {
          mod <- fit_lmer_safe_vis06(formula_obj, boot_data)

          if (is.null(mod)) {
            return(tibble())
          }

          pred_vec <- stats::predict(mod)
          metric_tbl <- compute_metrics_vis06(obs_vec, pred_vec)

          metric_tbl |>
            mutate(
              bootstrap = iter_i,
              model = model_name,
              .before = 1
            )
        })
      },
      .options = furrr::furrr_options(seed = seed)
    )
  } else {
    on.exit(future::plan("sequential"), add = TRUE)

    purrr::map_dfr(seq_len(n_boot), function(iter_i) {
      boot_data <- build_cluster_bootstrap_data_vis06(data)
      obs_vec <- if (isTRUE(log_outcome)) log(boot_data$sf36_pcs) else boot_data$sf36_pcs

      purrr::imap_dfr(formulas, function(formula_obj, model_name) {
        mod <- fit_lmer_safe_vis06(formula_obj, boot_data)

        if (is.null(mod)) {
          return(tibble())
        }

        pred_vec <- stats::predict(mod)
        metric_tbl <- compute_metrics_vis06(obs_vec, pred_vec)

        metric_tbl |>
          mutate(
            bootstrap = iter_i,
            model = model_name,
            .before = 1
          )
      })
    })
  }
}

summarize_bootstrap_metrics_vis06 <- function(boot_tbl) {
  boot_tbl |>
    group_by(model) |>
    summarise(
      rmse = mean(rmse, na.rm = TRUE),
      rmse_low = stats::quantile(rmse, probs = 0.025, na.rm = TRUE),
      rmse_high = stats::quantile(rmse, probs = 0.975, na.rm = TRUE),
      mae = mean(mae, na.rm = TRUE),
      mae_low = stats::quantile(mae, probs = 0.025, na.rm = TRUE),
      mae_high = stats::quantile(mae, probs = 0.975, na.rm = TRUE),
      r2 = mean(r2, na.rm = TRUE),
      r2_low = stats::quantile(r2, probs = 0.025, na.rm = TRUE),
      r2_high = stats::quantile(r2, probs = 0.975, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      rmse_display = sprintf("%.3f [%.3f; %.3f]", rmse, rmse_low, rmse_high),
      mae_display = sprintf("%.3f [%.3f; %.3f]", mae, mae_low, mae_high),
      r2_display = sprintf("%.3f [%.3f; %.3f]", r2, r2_low, r2_high)
    )
}

build_actual_vs_pred_plot_vis06 <- function(pred_df, ggir_parameter, log_outcome = TRUE) {
  outcome_label <- if (isTRUE(log_outcome)) "Observed log(sf36_pcs)" else "Observed sf36_pcs"
  pred_label <- if (isTRUE(log_outcome)) "Predicted log(sf36_pcs)" else "Predicted sf36_pcs"

  ggplot(pred_df, aes(x = observed, y = predicted, color = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.8, colour = "grey50") +
    geom_point(alpha = 0.3, size = 2.2) +
    stat_smooth(method = "lm", se = FALSE, linewidth = 1.0) +
    scale_color_manual(values = fixed_model_colors_vis06, drop = FALSE) +
    labs(
      title = "Actual vs predicted",
      subtitle = paste0("log(sf36_pcs) models for ", ggir_parameter),
      x = outcome_label,
      y = pred_label,
      color = "Model"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      plot.subtitle = element_text(size = 16, face = "italic"),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14, face = "bold")
    )
}

build_residual_plot_vis06 <- function(pred_df, ggir_parameter) {
  ggplot(pred_df, aes(x = model, y = residual, color = model)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.8, colour = "grey50") +
    ggbeeswarm::geom_quasirandom(alpha = 0.7, size = 2.0, width = 0.22) +
    scale_color_manual(values = fixed_model_colors_vis06, drop = FALSE) +
    labs(
      title = "Residual diagnostics",
      subtitle = paste0("Residual beeswarm for ", ggir_parameter),
      x = "Model",
      y = "Residuals",
      color = "Model"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      plot.subtitle = element_text(size = 16, face = "italic"),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
}

build_bootstrap_boxplot_vis06 <- function(bootstrap_tbl, ggir_parameter) {
  metric_long <- bootstrap_tbl |>
    select(model, bootstrap, rmse, mae, r2) |>
    pivot_longer(
      cols = c(rmse, mae, r2),
      names_to = "metric",
      values_to = "value"
    ) |>
    group_by(metric) |>
    mutate(
      value_norm = dplyr::case_when(
        metric %in% c("rmse", "mae") ~ scale01_vis06(value, higher_better = FALSE),
        metric == "r2" ~ scale01_vis06(value, higher_better = TRUE),
        TRUE ~ NA_real_
      )
    ) |>
    ungroup() |>
    mutate(
      metric = factor(metric, levels = c("rmse", "mae", "r2"), labels = c("RMSE", "MAE", "R2")),
      model = factor(model, levels = c("tis_only", "ggir_only", "tis_ggir"))
    )

  ggplot(metric_long, aes(x = metric, y = value_norm, fill = model, group = interaction(metric, model))) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.75, outlier.alpha = 0.25) +
    scale_fill_manual(values = fixed_model_colors_vis06, drop = FALSE) +
    labs(
      title = "Bootstrap performance",
      subtitle = paste0("Normalized RMSE, MAE and R2 for ", ggir_parameter),
      tag = "Higher normalized value = better model for RMSE, MAE and R2",
      x = "Metric",
      y = "Normalized value",
      fill = "Model"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      plot.subtitle = element_text(size = 16, face = "italic"),
      plot.tag = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14, face = "bold")
    )
}

build_multi_panel_plot_vis06 <- function(pred_df, bootstrap_tbl, ggir_parameter, log_outcome = TRUE) {
  p_actual <- build_actual_vs_pred_plot_vis06(pred_df, ggir_parameter, log_outcome = log_outcome)
  p_resid <- build_residual_plot_vis06(pred_df, ggir_parameter)
  p_boot <- build_bootstrap_boxplot_vis06(bootstrap_tbl, ggir_parameter)

  p_actual + p_resid + p_boot +
    patchwork::plot_layout(ncol = 3, widths = c(1.15, 1.35, 1.1))
}

export_parameter_results_vis06 <- function(
    compare_tbl,
    tidy_tbl,
    plot_obj,
    bootstrap_tbl,
    metrics_summary,
    diagnostic_plot,
    ggir_parameter,
    suffix = "01",
    width = 11,
    height = 8,
    dpi = 320
) {
  dir.create(file.path("output", "tis_ggir_eff", "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("output", "tis_ggir_eff", "figures"), recursive = TRUE, showWarnings = FALSE)

  prefix <- paste0(format(Sys.Date(), "%y%m%d"), "_", sanitize_name_vis06(ggir_parameter), "_tis_ggir_eff")
  xlsx_file <- file.path("output", "tis_ggir_eff", "tables", paste0(prefix, "_", suffix, ".xlsx"))
  plot_file <- file.path("output", "tis_ggir_eff", "figures", paste0(prefix, "_", suffix, ".png"))
  diagnostic_file <- file.path("output", "tis_ggir_eff", "figures", paste0(prefix, "_diagnostic_", suffix, ".png"))

  if (file.exists(xlsx_file)) {
    unlink(xlsx_file, force = TRUE)
  }

  if (file.exists(plot_file)) {
    unlink(plot_file, force = TRUE)
  }

  if (file.exists(diagnostic_file)) {
    unlink(diagnostic_file, force = TRUE)
  }

  rio::export(
    list(
      compare_performance = compare_tbl |>
        mutate(across(where(is.numeric), ~ round(.x, 4))),
      tidy_fixed = tidy_tbl |>
        mutate(across(where(is.numeric), ~ round(.x, 4))),
      bootstrap_metrics = metrics_summary |>
        mutate(across(where(is.numeric), ~ round(.x, 4)))
    ),
    xlsx_file
  )

  if (!is.null(plot_obj)) {
    grDevices::png(
      filename = plot_file,
      width = width * dpi,
      height = height * dpi,
      res = dpi,
      bg = "white",
      type = "cairo"
    )
    print(plot_obj)
    grDevices::dev.off()
  }

  if (!is.null(diagnostic_plot)) {
    grDevices::png(
      filename = diagnostic_file,
      width = 18 * dpi,
      height = 6.5 * dpi,
      res = dpi,
      bg = "white",
      type = "cairo"
    )
    print(diagnostic_plot)
    grDevices::dev.off()
  }

  list(
    xlsx_file = normalizePath(xlsx_file, winslash = "/", mustWork = FALSE),
    plot_file = normalizePath(plot_file, winslash = "/", mustWork = FALSE),
    diagnostic_file = normalizePath(diagnostic_file, winslash = "/", mustWork = FALSE)
  )
}

run_tis_ggir_eff_workflow_vis06 <- function(
    ggir_data = d18_ggir_data,
    tis_data = d18_tis_comp_01,
    ggir_parameters = c(
      "ig_gradient_enmo",
      "ig_intercept_enmo",
      "acc_day_spt_wei",
      "m5_enmo",
      "mvpa_t100_time_min",
      "p9931_enmo"
    ),
    suffix = "01",
    export_results = TRUE,
    bootstrap_parallel = TRUE,
    bootstrap_workers = safe_future_workers_vis06
) {
  prepared_data <- prepare_tis_ggir_eff_data_vis06(ggir_data, tis_data)
  ggir_available <- intersect(ggir_parameters, names(prepared_data))

  parameter_results <- purrr::map(ggir_available, function(ggir_parameter) {
    raw_bundle <- build_models_for_parameter_vis06(
      data = prepared_data,
      ggir_parameter = ggir_parameter,
      log_outcome = FALSE
    )

    log_bundle <- build_models_for_parameter_vis06(
      data = prepared_data,
      ggir_parameter = ggir_parameter,
      log_outcome = TRUE
    )

    compare_tbl <- compare_models_vis06(log_bundle$models, ggir_parameter)
    tidy_tbl <- tidy_models_vis06(log_bundle$models, ggir_parameter)
    plot_obj <- build_compare_plot_vis06(
      models_named = log_bundle$models,
      ggir_parameter = ggir_parameter,
      log_outcome = TRUE
    )
    pred_df <- build_prediction_df_vis06(
      models_named = log_bundle$models,
      data = log_bundle$data,
      log_outcome = TRUE
    )
    bootstrap_tbl <- bootstrap_model_metrics_vis06(
      formulas = log_bundle$formulas,
      data = log_bundle$data,
      log_outcome = TRUE,
      n_boot = 500,
      parallel = bootstrap_parallel,
      workers = bootstrap_workers
    )
    metrics_summary <- summarize_bootstrap_metrics_vis06(bootstrap_tbl)
    diagnostic_plot <- build_multi_panel_plot_vis06(
      pred_df = pred_df,
      bootstrap_tbl = bootstrap_tbl,
      ggir_parameter = ggir_parameter,
      log_outcome = TRUE
    )

    export_paths <- NULL
    if (isTRUE(export_results) && !is.null(compare_tbl)) {
      export_paths <- export_parameter_results_vis06(
        compare_tbl = compare_tbl,
        tidy_tbl = tidy_tbl,
        plot_obj = plot_obj,
        bootstrap_tbl = bootstrap_tbl,
        metrics_summary = metrics_summary,
        diagnostic_plot = diagnostic_plot,
        ggir_parameter = ggir_parameter,
        suffix = suffix
      )
    }

    list(
      compare_performance = compare_tbl,
      tidy = tidy_tbl,
      compare_plot = plot_obj,
      prediction_df = pred_df,
      bootstrap_metrics_raw = bootstrap_tbl,
      bootstrap_metrics_summary = metrics_summary,
      diagnostic_plot = diagnostic_plot,
      raw_models = raw_bundle$models,
      raw_formulas = raw_bundle$formulas,
      log_models = log_bundle$models,
      log_formulas = log_bundle$formulas,
      raw_data = raw_bundle$data,
      log_data = log_bundle$data,
      export_paths = export_paths
    )
  })

  names(parameter_results) <- ggir_available

  complete_models <- purrr::imap(parameter_results, function(res, nm) {
    list(
      ggir_parameter = nm,
      raw = list(
        formulas = res$raw_formulas,
        models = res$raw_models,
        data = res$raw_data
      ),
      log_sf36_pcs = list(
        formulas = res$log_formulas,
        models = res$log_models,
        data = res$log_data
      )
    )
  })

  list(
    prepared_data = prepared_data,
    merged_tis_data = merge_tis_into_ggir_vis06(ggir_data, tis_data),
    parameter_results = parameter_results,
    complete_models = complete_models
  )
}

if (exists("d18_ggir_data") && exists("d18_tis_comp_01")) {
  ggir_vis06_results <- run_tis_ggir_eff_workflow_vis06(
    ggir_data = d18_ggir_data,
    tis_data = d18_tis_comp_01,
    export_results = isTRUE(getOption("myo.ggir.vis06.export", TRUE))
  )

  future::plan("sequential")

  ggir_vis06_complete_models <- ggir_vis06_results$complete_models

  message(
    "vis06 exported to: ",
    normalizePath(file.path("output", "tis_ggir_eff"), winslash = "/", mustWork = FALSE)
  )
}
