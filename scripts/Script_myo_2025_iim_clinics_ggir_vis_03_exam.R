# ************************************************************
# SF-36 radar comparisons for GGIR / fi2 / mmt8 models
# Adds exam as a fixed effect to the original vis_03 workflow.
# ************************************************************

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  dplyr,
  purrr,
  tibble,
  stringr,
  fmsb,
  paletteer,
  easystats,
  rio,
  lme4
)

scale01_vis03 <- function(x, higher_better = TRUE) {
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

weighted_mean_na_vis03 <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)

  if (!any(ok)) {
    return(NA_real_)
  }

  sum(x[ok] * w[ok]) / sum(w[ok])
}

fixed_parameter_colors_vis03 <- stats::setNames(
  as.character(paletteer::paletteer_d("ggsci::default_igv", n = 8)),
  c(
    "fi2",
    "mmt8",
    "ig_gradient_enmo",
    "acc_day_spt_wei",
    "p9931_enmo",
    "m5_enmo",
    "mvpa_t100_time_min",
    "ig_intercept_enmo"
  )
)

get_fixed_colors_vis03 <- function(labels) {
  cols <- fixed_parameter_colors_vis03[labels]
  missing_idx <- is.na(cols)

  if (any(missing_idx)) {
    cols[missing_idx] <- grDevices::hcl.colors(sum(missing_idx), palette = "Dynamic")
  }

  unname(cols)
}

resolve_first_existing_vis03 <- function(data, candidates, label) {
  hit <- intersect(candidates, names(data))

  if (length(hit) == 0) {
    stop(
      "V datech chybi pozadovany sloupec pro `", label, "`. ",
      "Zkouseno: ", paste(candidates, collapse = ", "), "."
    )
  }

  hit[1]
}

prepare_ggir_sf36_data_vis03 <- function(data) {
  immet_col <- resolve_first_existing_vis03(data, c("immet_id"), "immet_id")
  exam_col <- resolve_first_existing_vis03(data, c("exam_order", "exam"), "exam")
  sex_col <- resolve_first_existing_vis03(data, c("sex", "pohlavi", "gender"), "sex")
  age_col <- resolve_first_existing_vis03(data, c("age_m0", "age"), "age_m0")
  fi2_col <- resolve_first_existing_vis03(data, c("fi2", "fi_2"), "fi2")
  mmt8_col <- resolve_first_existing_vis03(data, c("mmt8", "mmt8_total"), "mmt8")

  sf36_candidates <- c(
    "sf36_pcs", "sf36_mcs", "sf36pf", "sf36rp", "sf36bp",
    "sf36gh", "sf36vt", "sf36sf", "sf36re", "sf36mf"
  )
  sf36_cols <- intersect(sf36_candidates, names(data))

  ggir_candidates <- c(
    "ig_gradient_enmo",
    "ig_intercept_enmo",
    "acc_day_spt_wei",
    "m5_enmo",
    "mvpa_t100_time_min",
    "p9931_enmo"
  )
  ggir_cols <- intersect(ggir_candidates, names(data))

  data |>
    transmute(
      immet_id = as.character(.data[[immet_col]]),
      exam = as.character(.data[[exam_col]]),
      sex = as.factor(.data[[sex_col]]),
      age_m0 = as.numeric(.data[[age_col]]),
      fi2 = as.numeric(.data[[fi2_col]]),
      mmt8 = as.numeric(.data[[mmt8_col]]),
      across(all_of(sf36_cols), as.numeric),
      across(all_of(ggir_cols), as.numeric)
    ) |>
    mutate(
      exam = factor(
        exam,
        levels = unique(exam)[order(as.numeric(sub("^M", "", unique(exam))))]
      )
    )
}

compute_semipartial_r2_vis03 <- function(mod, focal_term) {
  reduced_mod <- stats::update(
    mod,
    formula = stats::as.formula(paste(". ~ . -", focal_term))
  )

  full_r2 <- performance::r2_nakagawa(mod, verbose = FALSE)
  reduced_r2 <- performance::r2_nakagawa(reduced_mod, verbose = FALSE)

  pmax(
    0,
    unname(full_r2$R2_marginal) - unname(reduced_r2$R2_marginal)
  )
}

extract_single_model_metrics_vis03 <- function(mod, response, focal_term, ci = 0.95) {
  std_tab <- parameters::standardize_parameters(
    mod,
    method = "refit",
    ci = ci,
    verbose = FALSE
  ) |>
    as.data.frame()

  std_col <- grep("^Std_", names(std_tab), value = TRUE)[1]

  if (is.na(std_col)) {
    stop("Nepodarilo se najit sloupec se standardizovanym koeficientem.")
  }

  term_std <- std_tab |>
    filter(.data$Parameter == focal_term)

  if (nrow(term_std) != 1) {
    stop(
      "Term `", focal_term, "` nebyl nalezen jednoznacne. ",
      "Skript podporuje jednoduche 1-df hlavni efekty."
    )
  }

  beta_std <- term_std[[std_col]][1]
  ci_low <- term_std$CI_low[1]
  ci_high <- term_std$CI_high[1]
  ci_width <- abs(ci_high - ci_low)

  support_prop <- dplyr::case_when(
    beta_std > 0 & ci_high > 0 ~ pmax(ci_high, 0) / ci_width,
    beta_std < 0 & ci_low < 0 ~ abs(pmin(ci_low, 0)) / ci_width,
    TRUE ~ 0.5
  )

  support_prop <- pmin(pmax(support_prop, 0), 1)
  r2_tbl <- performance::r2_nakagawa(mod, verbose = FALSE)

  tibble(
    response = response,
    parameter = focal_term,
    beta_std = beta_std,
    beta_abs = abs(beta_std),
    ci_low = ci_low,
    ci_high = ci_high,
    ci_width = ci_width,
    ci_support_prop = support_prop,
    ci_precision_raw = (1 / ci_width) * support_prop,
    semi_partial_r2 = compute_semipartial_r2_vis03(mod, focal_term),
    r2_marginal = unname(r2_tbl$R2_marginal),
    r2_conditional = unname(r2_tbl$R2_conditional),
    aic = stats::AIC(mod)
  )
}

score_model_metrics_vis03 <- function(raw_tbl) {
  raw_tbl |>
    mutate(
      score_beta = scale01_vis03(beta_abs, higher_better = TRUE),
      score_ci = scale01_vis03(ci_precision_raw, higher_better = TRUE),
      score_semi_partial_r2 = scale01_vis03(semi_partial_r2, higher_better = TRUE),
      score_r2_marginal = scale01_vis03(r2_marginal, higher_better = TRUE),
      score_r2_conditional = scale01_vis03(r2_conditional, higher_better = TRUE),
      score_aic = scale01_vis03(aic, higher_better = FALSE)
    ) |>
    rowwise() |>
    mutate(
      composite_score = weighted_mean_na_vis03(
        c(
          score_beta,
          score_ci,
          score_semi_partial_r2,
          score_r2_marginal,
          score_r2_conditional,
          score_aic
        ),
        c(1, 1, 1, 1, 1, 1)
      )
    ) |>
    ungroup() |>
    arrange(desc(composite_score)) |>
    mutate(rank = row_number())
}

format_model_table_vis03 <- function(score_tbl) {
  score_tbl |>
    select(
      rank,
      parameter,
      response,
      beta_std,
      ci_low,
      ci_high,
      ci_width,
      ci_support_prop,
      semi_partial_r2,
      r2_marginal,
      r2_conditional,
      aic,
      score_beta,
      score_ci,
      score_semi_partial_r2,
      score_r2_marginal,
      score_r2_conditional,
      score_aic,
      composite_score
    ) |>
    rename(
      Parameter = parameter,
      Outcome = response,
      `Std. beta` = beta_std,
      `CI low` = ci_low,
      `CI high` = ci_high,
      `CI width` = ci_width,
      `CI support proportion` = ci_support_prop,
      `R2 - semi-partial` = semi_partial_r2,
      `R2 marginal` = r2_marginal,
      `R2 conditional` = r2_conditional,
      AIC = aic,
      `Score |beta|` = score_beta,
      `Score CI precision` = score_ci,
      `Score R2 - semi-partial` = score_semi_partial_r2,
      `Score R2 marginal` = score_r2_marginal,
      `Score R2 conditional` = score_r2_conditional,
      `Score AIC` = score_aic,
      `Composite score` = composite_score
    ) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
}

build_fmsb_radar_vis03 <- function(score_tbl, title, subtitle) {
  radar_metrics <- c(
    "score_beta",
    "score_ci",
    "score_semi_partial_r2",
    "score_r2_marginal",
    "score_r2_conditional",
    "score_aic"
  )

  pretty_labels <- c(
    score_beta = "|beta| std.",
    score_ci = "CI precision",
    score_semi_partial_r2 = "R2 - semi-partial",
    score_r2_marginal = "R2 marginal",
    score_r2_conditional = "R2 conditional",
    score_aic = "AIC (reverse)"
  )

  radar_df <- score_tbl |>
    select(parameter, all_of(radar_metrics)) |>
    as.data.frame()

  rownames(radar_df) <- radar_df$parameter
  radar_df$parameter <- NULL
  radar_df[] <- lapply(radar_df, as.numeric)
  names(radar_df) <- pretty_labels[radar_metrics]

  radar_df <- rbind(
    rep(1, ncol(radar_df)),
    rep(0, ncol(radar_df)),
    radar_df
  ) |>
    as.data.frame()

  rownames(radar_df)[1:2] <- c("max", "min")

  structure(
    list(
      data = radar_df,
      axis_labels = pretty_labels[radar_metrics],
      legend_labels = score_tbl$parameter,
      title = title,
      subtitle = subtitle
    ),
    class = "fmsb_radar_vis03"
  )
}

print.fmsb_radar_vis03 <- function(x, ...) {
  title_cex <- 60 / 12
  subtitle_cex <- 50 / 12
  axis_cex <- 32 / 12
  y_axis_cex <- 36 / 12
  legend_cex <- 28 / 12
  point_cex <- 36 / 12

  cols <- get_fixed_colors_vis03(x$legend_labels)
  fill_cols <- grDevices::adjustcolor(cols, alpha.f = 0.18)

  op <- graphics::par(
    mar = c(5, 7, 8, 7),
    plt = c(0.08, 0.50, 0.12, 0.80),
    xpd = TRUE,
    font = 2
  )
  on.exit(graphics::par(op), add = TRUE)

  fmsb::radarchart(
    x$data,
    vlabels = x$axis_labels,
    axistype = 1,
    seg = 4,
    pcol = cols,
    pfcol = fill_cols,
    plwd = 2,
    plty = 1,
    cglcol = "grey85",
    cglty = 1,
    cglwd = 1,
    axislabcol = "grey50",
    vlcex = axis_cex,
    calcex = y_axis_cex
  )

  graphics::title(main = x$title, cex.main = title_cex, font.main = 2, line = 3)
  graphics::mtext(x$subtitle, side = 3, line = 1.2, cex = subtitle_cex, font = 3)
  graphics::legend(
    "topright",
    inset = c(-0.34, 0.02),
    legend = x$legend_labels,
    bty = "n",
    pch = 16,
    col = cols,
    pt.cex = point_cex,
    cex = legend_cex,
    text.font = 2
  )

  invisible(x)
}

compare_adjusted_models_vis03 <- function(data, response, base_predictor, tested_parameters, ci = 0.95) {
  models <- map(tested_parameters, function(parameter_name) {
    model_data <- data |>
      select(immet_id, exam, age_m0, sex, all_of(response), all_of(base_predictor), all_of(parameter_name)) |>
      na.omit()

    fml <- stats::as.formula(
      paste0(
        response, " ~ ", base_predictor, " + ", parameter_name,
        " + exam + age_m0 + sex + (1 | immet_id)"
      )
    )

    lme4::lmer(fml, data = model_data, REML = FALSE)
  })

  names(models) <- tested_parameters

  raw_tbl <- imap_dfr(models, function(mod, parameter_name) {
    extract_single_model_metrics_vis03(
      mod = mod,
      response = response,
      focal_term = parameter_name,
      ci = ci
    )
  })

  score_tbl <- score_model_metrics_vis03(raw_tbl)

  list(
    models = models,
    raw_metrics = raw_tbl,
    composite_table = score_tbl,
    composite_table_export = format_model_table_vis03(score_tbl),
    radar_plot = build_fmsb_radar_vis03(
      score_tbl = score_tbl,
      title = paste0("Radar comparison of focal predictors: ", response),
      subtitle = paste0(
        response, " ~ ", base_predictor,
        " + tested variable + exam + age_m0 + sex + (1 | immet_id)"
      )
    )
  )
}

compare_single_predictor_models_vis03 <- function(data, response, tested_parameters, ci = 0.95) {
  models <- map(tested_parameters, function(parameter_name) {
    model_data <- data |>
      select(immet_id, exam, age_m0, sex, all_of(response), all_of(parameter_name)) |>
      na.omit()

    fml <- stats::as.formula(
      paste0(response, " ~ ", parameter_name, " + exam + age_m0 + sex + (1 | immet_id)")
    )

    lme4::lmer(fml, data = model_data, REML = FALSE)
  })

  names(models) <- tested_parameters

  raw_tbl <- imap_dfr(models, function(mod, parameter_name) {
    extract_single_model_metrics_vis03(
      mod = mod,
      response = response,
      focal_term = parameter_name,
      ci = ci
    )
  })

  score_tbl <- score_model_metrics_vis03(raw_tbl)

  list(
    models = models,
    raw_metrics = raw_tbl,
    composite_table = score_tbl,
    composite_table_export = format_model_table_vis03(score_tbl),
    radar_plot = build_fmsb_radar_vis03(
      score_tbl = score_tbl,
      title = paste0("Radar comparison of focal predictors: ", response),
      subtitle = paste0(
        response, " ~ tested variable + exam + age_m0 + sex + (1 | immet_id)"
      )
    )
  )
}

export_result_vis03 <- function(result, table_file, plot_file, width = 14, height = 8, dpi = 320) {
  dir.create(dirname(table_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(plot_file), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(table_file)) {
    unlink(table_file, force = TRUE)
  }

  if (file.exists(plot_file)) {
    unlink(plot_file, force = TRUE)
  }

  rio::export(
    list(
      composite_table = result$composite_table_export,
      raw_metrics = result$raw_metrics |>
        mutate(across(where(is.numeric), ~ round(.x, 4)))
    ),
    table_file
  )

  grDevices::png(
    filename = plot_file,
    width = width * dpi,
    height = height * dpi,
    res = dpi,
    bg = "white",
    type = "cairo"
  )
  print(result$radar_plot)
  grDevices::dev.off()

  invisible(
    list(
      table_file = normalizePath(table_file, winslash = "/", mustWork = FALSE),
      plot_file = normalizePath(plot_file, winslash = "/", mustWork = FALSE)
    )
  )
}

sanitize_name_vis03 <- function(x) {
  stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
}

run_sf36_workflow_vis03_exam <- function(
    data = d18_ggir_data,
    sf36_outcomes = c(
      "sf36_pcs", "sf36_mcs", "sf36pf", "sf36rp", "sf36bp",
      "sf36gh", "sf36vt", "sf36sf", "sf36re", "sf36mf"
    ),
    ggir_parameters = c(
      "ig_gradient_enmo",
      "ig_intercept_enmo",
      "acc_day_spt_wei",
      "m5_enmo",
      "mvpa_t100_time_min",
      "p9931_enmo"
    ),
    single_parameters = c(
      "fi2", "mmt8",
      "ig_gradient_enmo",
      "ig_intercept_enmo",
      "acc_day_spt_wei",
      "m5_enmo",
      "mvpa_t100_time_min",
      "p9931_enmo"
    ),
    export_results = TRUE,
    suffix = "01"
) {
  prepared_data <- prepare_ggir_sf36_data_vis03(data)
  sf36_available <- intersect(sf36_outcomes, names(prepared_data))
  ggir_available <- intersect(ggir_parameters, names(prepared_data))
  single_available <- intersect(single_parameters, names(prepared_data))

  dir.create(file.path("output", "exam_incl", "sf36", "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("output", "exam_incl", "sf36", "figures"), recursive = TRUE, showWarnings = FALSE)

  results <- map(sf36_available, function(outcome_name) {
    res_fi2 <- compare_adjusted_models_vis03(
      data = prepared_data,
      response = outcome_name,
      base_predictor = "fi2",
      tested_parameters = ggir_available
    )

    res_mmt8 <- compare_adjusted_models_vis03(
      data = prepared_data,
      response = outcome_name,
      base_predictor = "mmt8",
      tested_parameters = ggir_available
    )

    res_single <- compare_single_predictor_models_vis03(
      data = prepared_data,
      response = outcome_name,
      tested_parameters = single_available
    )

    export_paths <- NULL

    if (isTRUE(export_results)) {
      prefix <- paste0(format(Sys.Date(), "%y%m%d"), "_", sanitize_name_vis03(outcome_name), "_sf36")

      export_paths <- list(
        fi2 = export_result_vis03(
          result = res_fi2,
          table_file = file.path("output", "exam_incl", "sf36", "tables", paste0(prefix, "_fi2adj_exam_", suffix, ".xlsx")),
          plot_file = file.path("output", "exam_incl", "sf36", "figures", paste0(prefix, "_fi2adj_exam_", suffix, ".png"))
        ),
        mmt8 = export_result_vis03(
          result = res_mmt8,
          table_file = file.path("output", "exam_incl", "sf36", "tables", paste0(prefix, "_mmt8adj_exam_", suffix, ".xlsx")),
          plot_file = file.path("output", "exam_incl", "sf36", "figures", paste0(prefix, "_mmt8adj_exam_", suffix, ".png"))
        ),
        single = export_result_vis03(
          result = res_single,
          table_file = file.path("output", "exam_incl", "sf36", "tables", paste0(prefix, "_single_exam_", suffix, ".xlsx")),
          plot_file = file.path("output", "exam_incl", "sf36", "figures", paste0(prefix, "_single_exam_", suffix, ".png"))
        )
      )
    }

    list(
      fi2_adjusted = res_fi2,
      mmt8_adjusted = res_mmt8,
      single_predictor = res_single,
      export_paths = export_paths
    )
  })

  names(results) <- sf36_available
  results
}

if (exists("d18_ggir_data")) {
  ggir_vis03_exam_results <- run_sf36_workflow_vis03_exam(
    data = d18_ggir_data,
    export_results = isTRUE(getOption("myo.ggir.vis03.exam.export", TRUE))
  )

  message("vis03_exam exported to: ", normalizePath(file.path("output", "exam_incl", "sf36"), winslash = "/", mustWork = FALSE))
}
