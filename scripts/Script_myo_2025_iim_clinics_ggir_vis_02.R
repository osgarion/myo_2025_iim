# ************************************************************
# Triplet comparison of GGIR models for PhVAS
# For each tested parameter compare:
# 1) PhVAS ~ base_predictor + age_m0 + sex + (1 | immet_id)
# 2) PhVAS ~ tested_parameter + age_m0 + sex + (1 | immet_id)
# 3) PhVAS ~ base_predictor + tested_parameter + age_m0 + sex + (1 | immet_id)
# Supports: lmerMod
# ************************************************************

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  dplyr,
  tidyr,
  purrr,
  tibble,
  stringr,
  fmsb,
  paletteer,
  easystats,
  rio,
  lme4
)

scale01_vis02 <- function(x, higher_better = TRUE) {
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

weighted_mean_na_vis02 <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)

  if (!any(ok)) {
    return(NA_real_)
  }

  sum(x[ok] * w[ok]) / sum(w[ok])
}

fixed_parameter_colors_vis02 <- stats::setNames(
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

get_fixed_colors_vis02 <- function(labels) {
  cols <- fixed_parameter_colors_vis02[labels]
  missing_idx <- is.na(cols)

  if (any(missing_idx)) {
    cols[missing_idx] <- grDevices::hcl.colors(sum(missing_idx), palette = "Dynamic")
  }

  unname(cols)
}

resolve_first_existing_vis02 <- function(data, candidates, label) {
  hit <- intersect(candidates, names(data))

  if (length(hit) == 0) {
    stop(
      "V datech chybi pozadovany sloupec pro `", label, "`. ",
      "Zkoušeno: ", paste(candidates, collapse = ", "), "."
    )
  }

  hit[1]
}

prepare_ggir_triplet_data_vis02 <- function(data) {
  immet_col <- resolve_first_existing_vis02(data, c("immet_id"), "immet_id")
  exam_col <- resolve_first_existing_vis02(data, c("exam_order", "exam"), "exam")
  sex_col <- resolve_first_existing_vis02(data, c("sex", "pohlavi", "gender"), "sex")
  age_col <- resolve_first_existing_vis02(data, c("age_m0", "age"), "age_m0")
  phvas_col <- resolve_first_existing_vis02(data, c("ph_vas", "physician_vas", "PhVAS"), "PhVAS")
  mmt8_col <- resolve_first_existing_vis02(data, c("mmt8", "mmt8_total"), "mmt8")
  fi2_col <- resolve_first_existing_vis02(data, c("fi2", "fi_2"), "fi2")

  parameter_candidates <- c(
    "ig_gradient_enmo",
    "ig_intercept_enmo",
    "acc_day_spt_wei",
    "m5_enmo",
    "mvpa_t100_time_min",
    "p9931_enmo"
  )

  parameter_cols <- intersect(parameter_candidates, names(data))

  if (length(parameter_cols) != length(parameter_candidates)) {
    stop(
      "V `d18_ggir_data` chybi nektere pozadovane GGIR parametry. ",
      "Nalezeno: ", paste(parameter_cols, collapse = ", "), "."
    )
  }

  data |>
    transmute(
      immet_id = as.character(.data[[immet_col]]),
      exam_order = as.character(.data[[exam_col]]),
      sex = as.factor(.data[[sex_col]]),
      age_m0 = as.numeric(.data[[age_col]]),
      PhVAS = as.numeric(.data[[phvas_col]]),
      mmt8 = as.numeric(.data[[mmt8_col]]),
      fi2 = as.numeric(.data[[fi2_col]]),
      across(all_of(parameter_cols), as.numeric)
    ) |>
    mutate(
      exam_order = factor(
        exam_order,
        levels = unique(exam_order)[order(as.numeric(sub("^M", "", unique(exam_order))))]
      )
    )
}

extract_single_term_metrics_vis02 <- function(mod, term_name, ci = 0.95) {
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
    filter(.data$Parameter == term_name)

  if (nrow(term_std) != 1) {
    stop(
      "Term `", term_name, "` nebyl nalezen jednoznacne. ",
      "Skript predpoklada jednoduche 1-df hlavni efekty."
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

  tibble(
    beta_std = beta_std,
    beta_abs = abs(beta_std),
    ci_low = ci_low,
    ci_high = ci_high,
    ci_width = ci_width,
    ci_support_prop = support_prop,
    ci_precision_raw = (1 / ci_width) * support_prop
  )
}

fit_triplet_models_vis02 <- function(data, base_predictor, tested_parameter) {
  formulas <- list(
    base_only = stats::as.formula(
      paste0("PhVAS ~ ", base_predictor, " + age_m0 + sex + (1 | immet_id)")
    ),
    tested_only = stats::as.formula(
      paste0("PhVAS ~ ", tested_parameter, " + age_m0 + sex + (1 | immet_id)")
    ),
    combined = stats::as.formula(
      paste0(
        "PhVAS ~ ", base_predictor, " + ", tested_parameter,
        " + age_m0 + sex + (1 | immet_id)"
      )
    )
  )

  target_terms <- c(
    base_only = base_predictor,
    tested_only = tested_parameter,
    combined = tested_parameter
  )

  model_labels <- c(
    base_only = paste(base_predictor, "only"),
    tested_only = paste(tested_parameter, "only"),
    combined = paste(base_predictor, "+", tested_parameter)
  )

  models <- imap(formulas, function(fml, nm) {
    needed_terms <- unique(c("PhVAS", "age_m0", "sex", "immet_id", all.vars(fml)))
    model_data <- data |>
      select(any_of(needed_terms)) |>
      na.omit()

    lme4::lmer(fml, data = model_data, REML = FALSE)
  })

  list(
    models = models,
    target_terms = target_terms,
    formulas = formulas,
    model_labels = model_labels
  )
}

extract_triplet_metrics_vis02 <- function(model_bundle, base_predictor, tested_parameter, ci = 0.95) {
  imap_dfr(model_bundle$models, function(mod, model_id) {
    term_name <- unname(model_bundle$target_terms[[model_id]])
    term_metrics <- extract_single_term_metrics_vis02(mod, term_name = term_name, ci = ci)
    r2_tbl <- performance::r2_nakagawa(mod, verbose = FALSE)

    tibble(
      comparison_set = paste(base_predictor, tested_parameter, sep = "__"),
      base_predictor = base_predictor,
      tested_parameter = tested_parameter,
      model_id = model_id,
      model_label = unname(model_bundle$model_labels[[model_id]]),
      focus_term = term_name,
      formula_label = paste(deparse(stats::formula(mod)), collapse = " "),
      aic = stats::AIC(mod),
      r2_marginal = unname(r2_tbl$R2_marginal),
      r2_conditional = unname(r2_tbl$R2_conditional)
    ) |>
      bind_cols(term_metrics)
  })
}

score_triplet_metrics_vis02 <- function(raw_tbl) {
  raw_tbl |>
    mutate(
      score_beta = scale01_vis02(beta_abs, higher_better = TRUE),
      score_ci = scale01_vis02(ci_precision_raw, higher_better = TRUE),
      score_r2_marginal = scale01_vis02(r2_marginal, higher_better = TRUE),
      score_r2_conditional = scale01_vis02(r2_conditional, higher_better = TRUE),
      score_aic = scale01_vis02(aic, higher_better = FALSE)
    ) |>
    rowwise() |>
    mutate(
      composite_score = weighted_mean_na_vis02(
        c(score_beta, score_ci, score_r2_marginal, score_r2_conditional, score_aic),
        c(1, 1, 1, 1, 1)
      )
    ) |>
    ungroup() |>
    arrange(desc(composite_score)) |>
    mutate(rank = row_number())
}

format_triplet_table_vis02 <- function(score_tbl) {
  score_tbl |>
    select(
      rank,
      base_predictor,
      tested_parameter,
      model_id,
      model_label,
      focus_term,
      beta_std,
      ci_low,
      ci_high,
      ci_width,
      ci_support_prop,
      r2_marginal,
      r2_conditional,
      aic,
      score_beta,
      score_ci,
      score_r2_marginal,
      score_r2_conditional,
      score_aic,
      composite_score
    ) |>
    rename(
      `Base predictor` = base_predictor,
      `Tested parameter` = tested_parameter,
      `Model ID` = model_id,
      `Model label` = model_label,
      `Focus term` = focus_term,
      `Std. beta` = beta_std,
      `CI low` = ci_low,
      `CI high` = ci_high,
      `CI width` = ci_width,
      `CI support proportion` = ci_support_prop,
      `R2 marginal` = r2_marginal,
      `R2 conditional` = r2_conditional,
      AIC = aic,
      `Score |beta|` = score_beta,
      `Score CI precision` = score_ci,
      `Score R2 marginal` = score_r2_marginal,
      `Score R2 conditional` = score_r2_conditional,
      `Score AIC` = score_aic,
      `Composite score` = composite_score
    ) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
}

build_fmsb_triplet_radar_vis02 <- function(score_tbl) {
  radar_metrics <- c(
    "score_beta",
    "score_ci",
    "score_r2_marginal",
    "score_r2_conditional",
    "score_aic"
  )

  pretty_labels <- c(
    score_beta = "|beta| std.",
    score_ci = "CI precision",
    score_r2_marginal = "R2 marginal",
    score_r2_conditional = "R2 conditional",
    score_aic = "AIC (reverse)"
  )

  radar_df <- score_tbl |>
    select(model_label, all_of(radar_metrics)) |>
    as.data.frame()

  rownames(radar_df) <- radar_df$model_label
  radar_df$model_label <- NULL
  radar_df[] <- lapply(radar_df, as.numeric)
  names(radar_df) <- pretty_labels[radar_metrics]

  radar_df <- rbind(
    rep(1, ncol(radar_df)),
    rep(0, ncol(radar_df)),
    radar_df
  ) |>
    as.data.frame()

  rownames(radar_df)[1:2] <- c("max", "min")

  base_predictor <- unique(score_tbl$base_predictor)
  tested_parameter <- unique(score_tbl$tested_parameter)

  structure(
    list(
      data = radar_df,
      legend_labels = score_tbl$model_label,
      title = paste0("Triplet model comparison: ", base_predictor, " vs ", tested_parameter),
      subtitle = paste(
        "1) PhVAS ~", base_predictor, "+ age_m0 + sex + (1 | immet_id)",
        "\n2) PhVAS ~ tested variable + age_m0 + sex + (1 | immet_id)",
        "\n3) PhVAS ~", base_predictor, "+ tested variable + age_m0 + sex + (1 | immet_id)"
      )
    ),
    class = "fmsb_radar_vis02"
  )
}

print.fmsb_radar_vis02 <- function(x, ...) {
  title_cex <- 60 / 12
  subtitle_cex <- 50 / 12
  axis_cex <- 32 / 12
  y_axis_cex <- 36 / 12
  legend_cex <- 28 / 12
  point_cex <- 36 / 12

  cols <- get_fixed_colors_vis02(x$legend_labels)
  fill_cols <- grDevices::adjustcolor(cols, alpha.f = 0.18)

  op <- graphics::par(
    mar = c(3, 3, 10, 7),
    plt = c(0.08, 0.50, 0.10, 0.78),
    xpd = TRUE,
    font = 2
  )
  on.exit(graphics::par(op), add = TRUE)

  fmsb::radarchart(
    x$data,
    axistype = 1,
    seg = 4,
    pcol = cols,
    pfcol = fill_cols,
    plwd = 3,
    plty = 1,
    cglcol = "grey85",
    cglty = 1,
    cglwd = 1,
    axislabcol = "grey50",
    vlcex = axis_cex,
    calcex = y_axis_cex
  )

  graphics::title(main = x$title, cex.main = title_cex, font.main = 2, line = 5)
  graphics::mtext(x$subtitle, side = 3, line = 1.5, cex = subtitle_cex, font = 3)
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

compare_triplet_models_vis02 <- function(data, base_predictor, tested_parameter, ci = 0.95) {
  model_bundle <- fit_triplet_models_vis02(
    data = data,
    base_predictor = base_predictor,
    tested_parameter = tested_parameter
  )

  raw_tbl <- extract_triplet_metrics_vis02(
    model_bundle = model_bundle,
    base_predictor = base_predictor,
    tested_parameter = tested_parameter,
    ci = ci
  )

  score_tbl <- score_triplet_metrics_vis02(raw_tbl)

  list(
    radar_plot = build_fmsb_triplet_radar_vis02(score_tbl),
    composite_table = score_tbl,
    composite_table_export = format_triplet_table_vis02(score_tbl),
    raw_metrics = raw_tbl,
    formulas = model_bundle$formulas
  )
}

export_triplet_result_vis02 <- function(result, xlsx_file, plot_file, width = 9, height = 9, dpi = 320) {
  dir.create(dirname(xlsx_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(plot_file), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(xlsx_file)) {
    unlink(xlsx_file, force = TRUE)
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
    xlsx_file
  )

  grDevices::png(
    filename = plot_file,
    width = width * dpi,
    height = height * dpi,
    res = dpi,
    bg = "white"
  )
  print(result$radar_plot)
  grDevices::dev.off()

  invisible(
    list(
      xlsx_file = normalizePath(xlsx_file, winslash = "/", mustWork = FALSE),
      plot_file = normalizePath(plot_file, winslash = "/", mustWork = FALSE)
    )
  )
}

sanitize_name_vis02 <- function(x) {
  stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
}

run_triplet_set_vis02 <- function(data, base_predictor, parameters, export_results = TRUE, suffix = "01") {
  results <- map(parameters, function(tested_parameter) {
    compare_triplet_models_vis02(
      data = data,
      base_predictor = base_predictor,
      tested_parameter = tested_parameter
    )
  })
  names(results) <- parameters

  export_index <- imap_dfr(results, function(result, tested_parameter) {
    prefix <- paste0(
      format(Sys.Date(), "%y%m%d"),
      "_ggir_triplet_",
      sanitize_name_vis02(base_predictor),
      "_",
      sanitize_name_vis02(tested_parameter),
      "_",
      suffix
    )

    paths <- NULL

    if (isTRUE(export_results)) {
      paths <- export_triplet_result_vis02(
        result = result,
        xlsx_file = file.path("output", "tables", paste0(prefix, ".xlsx")),
        plot_file = file.path("output", "figures", paste0(prefix, ".png"))
      )
    }

    print(result$radar_plot)
    print(result$composite_table_export)

    tibble(
      base_predictor = base_predictor,
      tested_parameter = tested_parameter,
      table_file = if (is.null(paths)) NA_character_ else paths$xlsx_file,
      plot_file = if (is.null(paths)) NA_character_ else paths$plot_file
    )
  })

  summary_tbl <- bind_rows(map(results, "composite_table_export"), .id = "tested_parameter_key")

  summary_file <- file.path(
    "output",
    "tables",
    paste0(format(Sys.Date(), "%y%m%d"), "_ggir_triplet_", sanitize_name_vis02(base_predictor), "_summary_", suffix, ".xlsx")
  )

  if (isTRUE(export_results)) {
    if (file.exists(summary_file)) {
      unlink(summary_file, force = TRUE)
    }

    rio::export(
      list(
        export_index = export_index,
        all_triplets = summary_tbl
      ),
      summary_file
    )
  }

  list(
    base_predictor = base_predictor,
    results = results,
    export_index = export_index,
    summary_table = summary_tbl,
    summary_file = if (isTRUE(export_results)) normalizePath(summary_file, winslash = "/", mustWork = FALSE) else NA_character_
  )
}

prepare_ggir_single_data_vis02 <- function(data) {
  immet_col <- resolve_first_existing_vis02(data, c("immet_id"), "immet_id")
  exam_col <- resolve_first_existing_vis02(data, c("exam_order", "exam"), "exam")
  sex_col <- resolve_first_existing_vis02(data, c("sex", "pohlavi", "gender"), "sex")
  age_col <- resolve_first_existing_vis02(data, c("age_m0", "age"), "age_m0")
  phvas_col <- resolve_first_existing_vis02(data, c("ph_vas", "physician_vas", "PhVAS"), "PhVAS")
  mmt8_col <- resolve_first_existing_vis02(data, c("mmt8", "mmt8_total"), "mmt8")
  fi2_col <- resolve_first_existing_vis02(data, c("fi2", "fi_2"), "fi2")

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
      exam_order = as.character(.data[[exam_col]]),
      sex = as.factor(.data[[sex_col]]),
      age_m0 = as.numeric(.data[[age_col]]),
      PhVAS = as.numeric(.data[[phvas_col]]),
      fi2 = as.numeric(.data[[fi2_col]]),
      mmt8 = as.numeric(.data[[mmt8_col]]),
      across(all_of(ggir_cols), as.numeric)
    ) |>
    mutate(
      exam_order = factor(
        exam_order,
        levels = unique(exam_order)[order(as.numeric(sub("^M", "", unique(exam_order))))]
      )
    )
}

compute_semipartial_r2_vis02 <- function(mod, focal_term) {
  reduced_mod <- stats::update(
    mod,
    formula = stats::as.formula(paste(". ~ . -", focal_term))
  )

  full_r2 <- performance::r2_nakagawa(mod, verbose = FALSE)
  reduced_r2 <- performance::r2_nakagawa(reduced_mod, verbose = FALSE)

  pmax(0, unname(full_r2$R2_marginal) - unname(reduced_r2$R2_marginal))
}

extract_single_model_metrics_vis02 <- function(mod, focal_term, ci = 0.95) {
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
    parameter = focal_term,
    response = "PhVAS",
    beta_std = beta_std,
    beta_abs = abs(beta_std),
    ci_low = ci_low,
    ci_high = ci_high,
    ci_width = ci_width,
    ci_support_prop = support_prop,
    ci_precision_raw = (1 / ci_width) * support_prop,
    semi_partial_r2 = compute_semipartial_r2_vis02(mod, focal_term),
    r2_marginal = unname(r2_tbl$R2_marginal),
    r2_conditional = unname(r2_tbl$R2_conditional),
    aic = stats::AIC(mod)
  )
}

fit_single_predictor_models_vis02 <- function(data, parameters) {
  models <- map(parameters, function(parameter_name) {
    model_data <- data |>
      select(immet_id, age_m0, sex, PhVAS, all_of(parameter_name)) |>
      na.omit()

    fml <- stats::as.formula(
      paste0("PhVAS ~ ", parameter_name, " + age_m0 + sex + (1 | immet_id)")
    )

    lme4::lmer(fml, data = model_data, REML = FALSE)
  })

  names(models) <- parameters
  models
}

score_single_predictor_metrics_vis02 <- function(raw_tbl) {
  raw_tbl |>
    mutate(
      score_beta = scale01_vis02(beta_abs, higher_better = TRUE),
      score_ci = scale01_vis02(ci_precision_raw, higher_better = TRUE),
      score_semi_partial_r2 = scale01_vis02(semi_partial_r2, higher_better = TRUE),
      score_r2_marginal = scale01_vis02(r2_marginal, higher_better = TRUE),
      score_r2_conditional = scale01_vis02(r2_conditional, higher_better = TRUE),
      score_aic = scale01_vis02(aic, higher_better = FALSE)
    ) |>
    rowwise() |>
    mutate(
      composite_score = weighted_mean_na_vis02(
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

format_single_predictor_table_vis02 <- function(score_tbl) {
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

build_fmsb_radar_vis02 <- function(score_tbl) {
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
      title = "Radar comparison of focal predictors (lmer)",
      subtitle = "PhVAS ~ tested variable + age_m0 + sex + (1 | immet_id)"
    ),
    class = "fmsb_radar_vis02"
  )
}

print.fmsb_radar_vis02 <- function(x, ...) {
  title_cex <- 60 / 12
  subtitle_cex <- 50 / 12
  axis_cex <- 32 / 12
  y_axis_cex <- 36 / 12
  legend_cex <- 28 / 12
  point_cex <- 36 / 12

  cols <- get_fixed_colors_vis02(x$legend_labels)
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

compare_single_predictor_models_vis02 <- function(
    data,
    parameters = c(
      "fi2",
      "mmt8",
      "ig_gradient_enmo",
      "ig_intercept_enmo",
      "acc_day_spt_wei",
      "m5_enmo",
      "mvpa_t100_time_min",
      "p9931_enmo"
    ),
    ci = 0.95
) {
  prepared_data <- prepare_ggir_single_data_vis02(data)
  available_parameters <- intersect(parameters, names(prepared_data))
  models <- fit_single_predictor_models_vis02(prepared_data, available_parameters)

  raw_tbl <- imap_dfr(models, function(mod, parameter_name) {
    extract_single_model_metrics_vis02(
      mod = mod,
      focal_term = parameter_name,
      ci = ci
    )
  })

  score_tbl <- score_single_predictor_metrics_vis02(raw_tbl)

  list(
    prepared_data = prepared_data,
    models = models,
    raw_metrics = raw_tbl,
    composite_table = score_tbl,
    composite_table_export = format_single_predictor_table_vis02(score_tbl),
    radar_plot = build_fmsb_radar_vis02(score_tbl)
  )
}

export_single_predictor_vis02 <- function(
    result,
    prefix = paste0(format(Sys.Date(), "%y%m%d"), "_ggir_single_predictor_02"),
    suffix = "01",
    width = 14,
    height = 9,
    dpi = 320
) {
  xlsx_file <- file.path("output", "tables", paste0(prefix, "_", suffix, ".xlsx"))
  plot_file <- file.path("output", "figures", paste0(prefix, "_", suffix, ".png"))

  dir.create(dirname(xlsx_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(plot_file), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(xlsx_file)) {
    unlink(xlsx_file, force = TRUE)
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
    xlsx_file
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
      xlsx_file = normalizePath(xlsx_file, winslash = "/", mustWork = FALSE),
      plot_file = normalizePath(plot_file, winslash = "/", mustWork = FALSE)
    )
  )
}

run_ggir_vis02_workflow <- function(
    data = d18_ggir_data,
    parameters = c(
      "fi2",
      "mmt8",
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
  result <- compare_single_predictor_models_vis02(
    data = data,
    parameters = parameters
  )

  export_paths <- NULL

  if (isTRUE(export_results)) {
    export_paths <- export_single_predictor_vis02(
      result = result,
      suffix = suffix
    )
  }

  list(
    result = result,
    export_paths = export_paths
  )
}

if (exists("d18_ggir_data")) {
  ggir_vis02_results <- run_ggir_vis02_workflow(
    data = d18_ggir_data,
    export_results = isTRUE(getOption("myo.ggir.vis02.export", TRUE))
  )

  ggir_vis02_model_results <- ggir_vis02_results$result
  res_vis02 <- ggir_vis02_model_results

  print(ggir_vis02_model_results$radar_plot)
  print(ggir_vis02_model_results$composite_table_export)

  if (!is.null(ggir_vis02_results$export_paths)) {
    message("vis02 table: ", ggir_vis02_results$export_paths$xlsx_file)
    message("vis02 plot: ", ggir_vis02_results$export_paths$plot_file)
  }
}
