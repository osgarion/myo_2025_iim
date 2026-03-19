# ************************************************************
# Radar comparison of focal predictors across comparable models
# Supports: lm, lmerMod
# Author: ChatGPT
# ************************************************************

# instalace pacman pouze pokud chybí
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

# načtení balíčků (instalují se jen pokud chybí)
pacman::p_load(
  dplyr,
  tidyr,
  purrr,
  tibble,
  ggplot2,
  fmsb,
  paletteer,
  easystats,  # zahrnuje parameters, performance, effectsize, insight
  rio,
  lme4
)

# ************************************************************
# Helpers ----
# ************************************************************

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

scale01 <- function(x, higher_better = TRUE) {
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

weighted_mean_na <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  
  if (!any(ok)) {
    return(NA_real_)
  }
  
  sum(x[ok] * w[ok]) / sum(w[ok])
}

fixed_parameter_colors_vis <- stats::setNames(
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

get_fixed_colors_vis <- function(labels) {
  cols <- fixed_parameter_colors_vis[labels]
  missing_idx <- is.na(cols)

  if (any(missing_idx)) {
    cols[missing_idx] <- grDevices::hcl.colors(sum(missing_idx), palette = "Dynamic")
  }

  unname(cols)
}

radar_polygon_area <- function(x) {
  x <- as.numeric(x)
  
  if (anyNA(x)) {
    return(NA_real_)
  }
  
  n <- length(x)
  theta <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)]
  
  xy <- cbind(
    x * cos(theta),
    x * sin(theta)
  )
  
  xy <- rbind(xy, xy[1, , drop = FALSE])
  
  0.5 * abs(sum(
    xy[-1, 1] * xy[-nrow(xy), 2] -
      xy[-nrow(xy), 1] * xy[-1, 2]
  ))
}

detect_model_family <- function(mod) {
  if (inherits(mod, "lmerMod")) {
    return("lmer")
  }
  
  if (inherits(mod, "lm") && !inherits(mod, "glm")) {
    return("lm")
  }
  
  stop("Podporované jsou zatím jen modely tříd `lm` a `lmerMod`.")
}

get_fixed_terms <- function(mod) {
  tt <- insight::find_terms(mod)
  terms <- sort(unique(tt$conditional %||% character(0)))
  setdiff(terms, "(Intercept)")
}

get_random_terms <- function(mod) {
  tt <- insight::find_terms(mod)
  sort(unique(tt$random %||% character(0)))
}

compute_semipartial_r2_lmer <- function(mod, focal_term) {
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

# ************************************************************
# 1) Validation of model set ----
# ************************************************************

validate_model_set <- function(models, allow_one_alt_response = TRUE) {
  if (!is.list(models) || length(models) < 2) {
    stop("`models` musí být list alespoň dvou modelů.")
  }
  
  if (is.null(names(models)) || any(names(models) == "")) {
    stop("`models` musí být pojmenovaný list.")
  }
  
  families <- map_chr(models, detect_model_family)
  
  if (dplyr::n_distinct(families) != 1) {
    stop("Všechny modely musí být stejného typu (`lm` nebo `lmerMod`).")
  }
  
  family <- families[1]
  
  responses <- map_chr(models, ~ insight::find_response(.x, combine = TRUE)[1])
  response_tab <- sort(table(responses), decreasing = TRUE)
  
  # max one different response
  if (length(response_tab) > 2) {
    stop("Modely obsahují více než dvě různé závislé proměnné.")
  }
  
  if (length(response_tab) == 2 && min(response_tab) > 1) {
    stop("Odlišná závislá proměnná se vyskytuje ve více než jednom modelu.")
  }
  
  response_note <- "Všechny modely mají stejnou závislou proměnnou."
  
  if (length(response_tab) == 2) {
    response_note <- paste(
      "Jeden model má odlišnou závislou proměnnou.",
      "Výstup je pak vhodný jen deskriptivně."
    )
    
    if (!allow_one_alt_response) {
      stop(response_note)
    }
    
    warning(response_note, call. = FALSE)
  }
  
  nobs_vec <- map_dbl(models, stats::nobs)
  
  if (dplyr::n_distinct(nobs_vec) != 1) {
    warning(
      "Modely nemají stejný počet pozorování. To zhoršuje porovnatelnost AIC a R2.",
      call. = FALSE
    )
  }
  
  fixed_terms <- map(models, get_fixed_terms)
  common_fixed_terms <- reduce(fixed_terms, intersect)
  
  focal_terms <- imap_chr(fixed_terms, function(x, nm) {
    diff_terms <- setdiff(x, common_fixed_terms)
    
    if (length(diff_terms) != 1) {
      stop(
        "Model `", nm, "` se musí lišit právě jedním fixním termem. ",
        "Nalezeno: ", length(diff_terms), "."
      )
    }
    
    diff_terms
  })
  
  if (anyDuplicated(focal_terms)) {
    stop("Rozdílné termy musí být unikátní.")
  }
  
  if (family == "lmer") {
    random_terms <- map(models, get_random_terms)
    random_ref <- random_terms[[1]]
    
    same_random <- map_lgl(random_terms, ~ setequal(.x, random_ref))
    
    if (!all(same_random)) {
      stop("U `lmer` musí být random struktura ve všech modelech stejná.")
    }
    
    reml_flags <- map_lgl(models, lme4::isREML)
    
    if (any(reml_flags)) {
      stop(
        "U `lmer` modelů pro porovnání AIC mezi modely s různými fixními efekty ",
        "refituj modely s `REML = FALSE`."
      )
    }
  }
  
  tibble(
    model = names(models),
    family = family,
    response = responses,
    parameter = focal_terms,
    response_note = response_note
  )
}

# ************************************************************
# 2) Metric extraction ----
# ************************************************************

extract_focal_metrics <- function(
    models,
    validation_tbl,
    ci = 0.95,
    standardize_method = "refit"
) {
  map_dfr(validation_tbl$model, function(model_name) {
    mod <- models[[model_name]]
    focal_term <- validation_tbl$parameter[validation_tbl$model == model_name]
    family <- validation_tbl$family[validation_tbl$model == model_name]
    response <- validation_tbl$response[validation_tbl$model == model_name]
    response_note <- validation_tbl$response_note[validation_tbl$model == model_name]
    
    std_tab <- parameters::standardize_parameters(
      mod,
      method = standardize_method,
      ci = ci,
      verbose = FALSE
    ) |>
      as.data.frame()
    
    std_col <- grep("^Std_", names(std_tab), value = TRUE)[1]
    
    if (is.na(std_col)) {
      stop("Nepodařilo se najít sloupec se standardizovaným koeficientem.")
    }
    
    term_std <- std_tab |>
      filter(.data$Parameter == focal_term)
    
    if (nrow(term_std) != 1) {
      stop(
        "Rozdílný term `", focal_term, "` nebyl nalezen jednoznačně. ",
        "Funkce je určená hlavně pro jednoduché 1-df hlavní efekty ",
        "(ne faktory s více úrovněmi, interakce nebo spline termy)."
      )
    }
    
    beta_std <- term_std[[std_col]][1]
    ci_low <- term_std$CI_low[1]
    ci_high <- term_std$CI_high[1]
    ci_width <- abs(ci_high - ci_low)
    ci_crosses_zero <- ci_low <= 0 && ci_high >= 0
    
    # smoother CI penalty:
    # base = narrower CI is better
    # if CI crosses zero, reduce score according to how much of CI lies
    # on the "supporting" side of the estimated effect
    support_prop <- dplyr::case_when(
      beta_std > 0 & ci_high > 0 ~ pmax(ci_high, 0) / ci_width,
      beta_std < 0 & ci_low < 0 ~ abs(pmin(ci_low, 0)) / ci_width,
      TRUE ~ 0.5
    )
    
    support_prop <- pmin(pmax(support_prop, 0), 1)
    ci_precision_raw <- (1 / ci_width) * support_prop
    
    aic <- stats::AIC(mod)
    
    if (family == "lm") {
      sr2_tab <- effectsize::r2_semipartial(
        mod,
        type = "terms",
        ci = ci,
        alternative = "two.sided"
      ) |>
        as.data.frame()
      
      sr2_term <- sr2_tab |>
        filter(.data$Parameter == focal_term)
      
      if (nrow(sr2_term) != 1) {
        stop(
          "Nepodařilo se jednoznačně vytáhnout partial R2 pro term `", focal_term, "`."
        )
      }
      
      sr2_col <- setdiff(
        names(sr2_term)[vapply(sr2_term, is.numeric, logical(1))],
        c("CI", "CI_low", "CI_high")
      )[1]
      
      predictor_r2 <- sr2_term[[sr2_col]][1]
      
      tibble(
        model = model_name,
        family = family,
        response = response,
        parameter = focal_term,
        beta_std = beta_std,
        beta_abs = abs(beta_std),
        ci_low = ci_low,
        ci_high = ci_high,
        ci_width = ci_width,
        ci_crosses_zero = ci_crosses_zero,
        ci_support_prop = support_prop,
        ci_precision_raw = ci_precision_raw,
        predictor_r2 = predictor_r2,
        semi_partial_r2 = predictor_r2,
        r2_marginal = NA_real_,
        r2_conditional = NA_real_,
        aic = aic,
        response_note = response_note
      )
    } else {
      r2_tab <- performance::r2_nakagawa(mod, verbose = FALSE)
      semi_partial_r2 <- compute_semipartial_r2_lmer(mod, focal_term)
      
      tibble(
        model = model_name,
        family = family,
        response = response,
        parameter = focal_term,
        beta_std = beta_std,
        beta_abs = abs(beta_std),
        ci_low = ci_low,
        ci_high = ci_high,
        ci_width = ci_width,
        ci_crosses_zero = ci_crosses_zero,
        ci_support_prop = support_prop,
        ci_precision_raw = ci_precision_raw,
        predictor_r2 = NA_real_,
        semi_partial_r2 = semi_partial_r2,
        r2_marginal = unname(r2_tab$R2_marginal),
        r2_conditional = unname(r2_tab$R2_conditional),
        aic = aic,
        response_note = response_note
      )
    }
  })
}

# ************************************************************
# 3) Scoring + composite table ----
# ************************************************************

score_focal_metrics <- function(raw_tbl, weights = NULL) {
  family <- unique(raw_tbl$family)
  
  if (length(family) != 1) {
    stop("Všechny modely musí být stejné rodiny.")
  }
  
  family <- family[1]
  
  score_tbl <- raw_tbl |>
    mutate(
      score_beta = scale01(beta_abs, higher_better = TRUE),
      score_ci = scale01(ci_precision_raw, higher_better = TRUE),
      score_aic = scale01(aic, higher_better = FALSE),
      score_predictor_r2 = scale01(predictor_r2, higher_better = TRUE),
      score_semi_partial_r2 = scale01(semi_partial_r2, higher_better = TRUE),
      score_r2_marginal = scale01(r2_marginal, higher_better = TRUE),
      score_r2_conditional = scale01(r2_conditional, higher_better = TRUE)
    )
  
  if (family == "lm") {
    radar_metrics <- c(
      "score_beta",
      "score_ci",
      "score_semi_partial_r2",
      "score_predictor_r2",
      "score_aic"
    )
  } else {
    radar_metrics <- c(
      "score_beta",
      "score_ci",
      "score_semi_partial_r2",
      "score_r2_marginal",
      "score_r2_conditional",
      "score_aic"
    )
  }
  
  if (is.null(weights)) {
    weights <- rep(1, length(radar_metrics))
    names(weights) <- radar_metrics
  } else {
    if (!all(radar_metrics %in% names(weights))) {
      stop("`weights` musí být pojmenovaný vektor obsahující všechny radar metriky.")
    }
    
    weights <- weights[radar_metrics]
  }
  
  score_tbl <- score_tbl |>
    rowwise() |>
    mutate(
      composite_score = weighted_mean_na(
        c_across(all_of(radar_metrics)),
        weights[radar_metrics]
      ),
      radar_area = radar_polygon_area(
        c_across(all_of(radar_metrics))
      )
    ) |>
    ungroup() |>
    arrange(desc(composite_score)) |>
    mutate(rank = row_number())
  
  attr(score_tbl, "radar_metrics") <- radar_metrics
  score_tbl
}

# ************************************************************
# 4) Radar plot ----
# ************************************************************

plot_radar_focal_metrics <- function(score_tbl) {
  radar_metrics <- attr(score_tbl, "radar_metrics")
  family <- unique(score_tbl$family)
  
  pretty_labels <- c(
    score_beta = "|β| std.",
    score_ci = "CI precision",
    score_predictor_r2 = "partial R²",
    score_r2_marginal = "R² marginal",
    score_r2_conditional = "R² conditional",
    score_aic = "AIC (reverse)"
  )
  
  plot_tbl <- score_tbl |>
    select(parameter, all_of(radar_metrics)) |>
    pivot_longer(
      cols = all_of(radar_metrics),
      names_to = "metric",
      values_to = "score"
    ) |>
    mutate(
      metric = factor(
        metric,
        levels = radar_metrics,
        labels = pretty_labels[radar_metrics]
      )
    )
  
  plot_tbl_closed <- plot_tbl |>
    group_by(parameter) |>
    arrange(metric, .by_group = TRUE) |>
    group_modify(~ bind_rows(.x, .x[1, ])) |>
    ungroup()
  
  ggplot(
    plot_tbl_closed,
    aes(
      x = metric,
      y = score,
      group = parameter,
      colour = parameter,
      fill = parameter
    )
  ) +
    geom_polygon(alpha = 0.12, linewidth = 0.9) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    coord_polar() +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0.25, 0.50, 0.75, 1.00)
    ) +
    labs(
      title = ifelse(
        family == "lm",
        "Radar comparison of focal predictors (lm)",
        "Radar comparison of focal predictors (lmer)"
      ),
      subtitle = "Legenda obsahuje pouze rozdílné parametry",
      x = NULL,
      y = NULL,
      colour = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "right"
    )
}

# ************************************************************
# 5) Nice export table ----
# ************************************************************

format_composite_table_for_export <- function(composite_tbl) {
  family <- unique(composite_tbl$family)
  
  if (length(family) != 1) {
    stop("Tabulka musí obsahovat jen jednu rodinu modelů.")
  }
  
  family <- family[1]
  
  if (family == "lm") {
    out <- composite_tbl |>
      select(
        rank,
        parameter,
        response,
        beta_std,
        ci_low,
        ci_high,
        ci_width,
        ci_support_prop,
        predictor_r2,
        semi_partial_r2,
        aic,
        score_beta,
        score_ci,
        score_semi_partial_r2,
        score_predictor_r2,
        score_aic,
        composite_score,
        radar_area
      ) |>
      rename(
        Parameter = parameter,
        Outcome = response,
        `Std. beta` = beta_std,
        `CI low` = ci_low,
        `CI high` = ci_high,
        `CI width` = ci_width,
        `CI support proportion` = ci_support_prop,
        `Partial R2` = predictor_r2,
        `Semi-partial R2` = semi_partial_r2,
        AIC = aic,
        `Score |beta|` = score_beta,
        `Score CI precision` = score_ci,
        `Score Semi-partial R2` = score_semi_partial_r2,
        `Score Partial R2` = score_predictor_r2,
        `Score AIC` = score_aic,
        `Composite score` = composite_score,
        `Radar area` = radar_area
      )
  } else {
    out <- composite_tbl |>
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
        composite_score,
        radar_area
      ) |>
      rename(
        Parameter = parameter,
        Outcome = response,
        `Std. beta` = beta_std,
        `CI low` = ci_low,
        `CI high` = ci_high,
        `CI width` = ci_width,
        `CI support proportion` = ci_support_prop,
        `Semi-partial R2` = semi_partial_r2,
        `R2 marginal` = r2_marginal,
        `R2 conditional` = r2_conditional,
        AIC = aic,
        `Score |beta|` = score_beta,
        `Score CI precision` = score_ci,
        `Score Semi-partial R2` = score_semi_partial_r2,
        `Score R2 marginal` = score_r2_marginal,
        `Score R2 conditional` = score_r2_conditional,
        `Score AIC` = score_aic,
        `Composite score` = composite_score,
        `Radar area` = radar_area
      )
  }
  
  out |>
    mutate(
      across(where(is.numeric), ~ round(.x, 3))
    )
}

# ************************************************************
# 6) Main wrapper ----
# ************************************************************

compare_focal_predictors <- function(
    models,
    ci = 0.95,
    standardize_method = "refit",
    allow_one_alt_response = TRUE,
    weights = NULL
) {
  validation_tbl <- validate_model_set(
    models = models,
    allow_one_alt_response = allow_one_alt_response
  )
  
  raw_tbl <- extract_focal_metrics(
    models = models,
    validation_tbl = validation_tbl,
    ci = ci,
    standardize_method = standardize_method
  )
  
  composite_tbl <- score_focal_metrics(
    raw_tbl = raw_tbl,
    weights = weights
  )

  model_subtitle <- paste(deparse(stats::formula(models[[1]])), collapse = " ")
  model_subtitle <- gsub("\\s+", " ", model_subtitle)
  model_subtitle <- sub(
    validation_tbl$parameter[[1]],
    "tested variable",
    model_subtitle,
    fixed = TRUE
  )
  attr(composite_tbl, "model_subtitle") <- model_subtitle
  
  radar_plot <- plot_radar_focal_metrics(composite_tbl)
  
  list(
    radar_plot = radar_plot,
    composite_table = composite_tbl,
    composite_table_export = format_composite_table_for_export(composite_tbl),
    validation_table = validation_tbl,
    raw_metrics = raw_tbl
  )
}

# ************************************************************
# 7) Export ----
# ************************************************************

export_focal_comparison <- function(
    result,
    xlsx_file = file.path(
      "output",
      "tables",
      paste0(format(Sys.Date(), "%y%m%d"), "_radar_comparison_01.xlsx")
    ),
    plot_file = file.path(
      "output",
      "figures",
      paste0(format(Sys.Date(), "%y%m%d"), "_radar_comparison_01.png")
    ),
    width = 8,
    height = 8,
    dpi = 320
) {
  dir.create(dirname(xlsx_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(plot_file), recursive = TRUE, showWarnings = FALSE)

  export_list <- list(
    composite_table = result$composite_table_export,
    raw_metrics = result$raw_metrics |>
      mutate(across(where(is.numeric), ~ round(.x, 4))),
    validation = result$validation_table
  )
  
  rio::export(export_list, xlsx_file)
  
  ggplot2::ggsave(
    filename = plot_file,
    plot = result$radar_plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  invisible(
    list(
      xlsx_file = normalizePath(xlsx_file, winslash = "/", mustWork = FALSE),
      plot_file = normalizePath(plot_file, winslash = "/", mustWork = FALSE)
    )
  )
}

# ************************************************************
# 8) Example usage ----
# ************************************************************

# **********************************************************************
## LM example ----
# **********************************************************************
# models_lm <- list(
#   age = model_age,
#   crp = model_crp,
#   bmi = model_bmi,
#   ck  = model_ck
# )
#
# res_lm <- compare_focal_predictors(models_lm)
# res_lm$radar_plot
# res_lm$composite_table_export
#
# export_focal_comparison(
#   result = res_lm,
#   xlsx_file = "radar_comparison_lm.xlsx",
#   plot_file = "radar_comparison_lm.png"
# )

# **********************************************************************
## LMER example ----
# **********************************************************************
# IMPORTANT: use REML = FALSE
# models_lmer <- list(
#   age = model_age_lmer,
#   crp = model_crp_lmer,
#   bmi = model_bmi_lmer,
#   ck  = model_ck_lmer
# )
#
# res_lmer <- compare_focal_predictors(models_lmer)
# res_lmer$radar_plot
# res_lmer$composite_table_export
#
# export_focal_comparison(
#   result = res_lmer,
#   xlsx_file = "radar_comparison_lmer.xlsx",
#   plot_file = "radar_comparison_lmer.png"
# )

# ************************************************************
# 9) Optional custom weights ----
# ************************************************************

# Example for lm:
# weights_lm <- c(
#   score_beta = 0.35,
#   score_ci = 0.25,
#   score_predictor_r2 = 0.25,
#   score_aic = 0.15
# )
#
# res_lm_w <- compare_focal_predictors(
#   models = models_lm,
#   weights = weights_lm
# )

# Example for lmer:
# weights_lmer <- c(
#   score_beta = 0.30,
#   score_ci = 0.25,
#   score_r2_marginal = 0.20,
#   score_r2_conditional = 0.15,
#   score_aic = 0.10
# )
#
# res_lmer_w <- compare_focal_predictors(
#   models = models_lmer,
#   weights = weights_lmer
# )

# ************************************************************
# 10) Application to d18_ggir_data ----
# ************************************************************

resolve_first_existing <- function(data, candidates, label) {
  hit <- intersect(candidates, names(data))

  if (length(hit) == 0) {
    stop(
      "V datech chybí požadovaný sloupec pro `", label, "`. ",
      "Zkoušeno: ", paste(candidates, collapse = ", "), "."
    )
  }

  hit[1]
}

prepare_ggir_model_data <- function(data) {
  immet_col <- resolve_first_existing(data, c("immet_id"), "immet_id")
  exam_col <- resolve_first_existing(data, c("exam_order", "exam"), "exam")
  sex_col <- resolve_first_existing(data, c("sex", "pohlavi", "gender"), "sex")
  age_col <- resolve_first_existing(data, c("age_m0", "age"), "age_m0")
  phvas_col <- resolve_first_existing(
    data,
    c("ph_vas", "physician_vas", "PhVAS"),
    "PhVAS"
  )
  mmt8_col <- resolve_first_existing(data, c("mmt8", "mmt8_total"), "mmt8")
  fi2_col <- resolve_first_existing(data, c("fi2", "fi_2"), "fi2")

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
      "V `d18_ggir_data` chybí některé požadované GGIR parametry. ",
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

build_focal_lmer_models <- function(data, base_predictor, parameters) {
  models <- map(parameters, function(parameter_name) {
    model_data <- data |>
      select(immet_id, exam_order, sex, age_m0, PhVAS, all_of(base_predictor), all_of(parameter_name)) |>
      na.omit()

    stats::formula(
      paste0(
        "PhVAS ~ ", base_predictor, " + ", parameter_name,
        " + age_m0 + sex + (1 | immet_id)"
      )
    ) |>
      lme4::lmer(data = model_data, REML = FALSE)
  })

  names(models) <- parameters
  models
}

collect_best_model_row <- function(result, base_predictor) {
  result$composite_table_export |>
    slice(1) |>
    mutate(base_predictor = base_predictor, .before = 1)
}

export_ggir_focal_results <- function(
    res_mmt8,
    res_fi2,
    best_models,
    prefix = paste0(format(Sys.Date(), "%y%m%d"), "_ggir_phvas_focal"),
    suffix = "01"
) {
  dir.create(file.path("output", "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("output", "figures"), recursive = TRUE, showWarnings = FALSE)

  mmt8_paths <- export_focal_comparison(
    result = res_mmt8,
    xlsx_file = file.path("output", "tables", paste0(prefix, "_mmt8_", suffix, ".xlsx")),
    plot_file = file.path("output", "figures", paste0(prefix, "_mmt8_", suffix, ".png"))
  )

  fi2_paths <- export_focal_comparison(
    result = res_fi2,
    xlsx_file = file.path("output", "tables", paste0(prefix, "_fi2_", suffix, ".xlsx")),
    plot_file = file.path("output", "figures", paste0(prefix, "_fi2_", suffix, ".png"))
  )

  best_models_file <- file.path(
    "output",
    "tables",
    paste0(prefix, "_best_models_", suffix, ".xlsx")
  )

  rio::export(
    list(
      best_models = best_models,
      mmt8_composite = res_mmt8$composite_table_export,
      fi2_composite = res_fi2$composite_table_export
    ),
    best_models_file
  )

  invisible(
    list(
      mmt8 = mmt8_paths,
      fi2 = fi2_paths,
      best_models_file = normalizePath(best_models_file, winslash = "/", mustWork = FALSE)
    )
  )
}

# Example workflow for current project data:
# source("scripts/functions/OBJ_01.R")
# ggir_model_data <- prepare_ggir_model_data(d18_ggir_data)
# ggir_parameters <- c(
#   "ig_gradient_enmo",
#   "ig_intercept_enmo",
#   "dur_day_spt_wei",
#   "m5_enmo",
#   "mvpa_t100_time_min",
#   "p9931_enmo"
# )
# models_mmt8 <- build_focal_lmer_models(
#   data = ggir_model_data,
#   base_predictor = "mmt8",
#   parameters = ggir_parameters
# )
# models_fi2 <- build_focal_lmer_models(
#   data = ggir_model_data,
#   base_predictor = "fi2",
#   parameters = ggir_parameters
# )
# res_mmt8 <- compare_focal_predictors(models_mmt8)
# res_fi2 <- compare_focal_predictors(models_fi2)
# best_models_export <- bind_rows(
#   collect_best_model_row(res_mmt8, "mmt8"),
#   collect_best_model_row(res_fi2, "fi2")
# )
# res_mmt8$radar_plot
# res_mmt8$composite_table_export
# res_fi2$radar_plot
# res_fi2$composite_table_export
# best_models_export
# export_ggir_focal_results(
#   res_mmt8 = res_mmt8,
#   res_fi2 = res_fi2,
#   best_models = best_models_export
# )

run_ggir_focal_workflow <- function(
    data = d18_ggir_data,
    parameters = c(
      "ig_gradient_enmo",
      "ig_intercept_enmo",
      "acc_day_spt_wei",
      "m5_enmo",
      "mvpa_t100_time_min",
      "p9931_enmo"
    ),
    export_results = FALSE,
    prefix = paste0(format(Sys.Date(), "%y%m%d"), "_ggir_phvas_focal"),
    suffix = "01"
) {
  ggir_model_data <- prepare_ggir_model_data(data)

  models_mmt8 <- build_focal_lmer_models(
    data = ggir_model_data,
    base_predictor = "mmt8",
    parameters = parameters
  )

  models_fi2 <- build_focal_lmer_models(
    data = ggir_model_data,
    base_predictor = "fi2",
    parameters = parameters
  )

  res_mmt8 <- compare_focal_predictors(models_mmt8)
  res_fi2 <- compare_focal_predictors(models_fi2)

  best_models_export <- bind_rows(
    collect_best_model_row(res_mmt8, "mmt8"),
    collect_best_model_row(res_fi2, "fi2")
  )

  export_paths <- NULL

  if (isTRUE(export_results)) {
    export_paths <- export_ggir_focal_results(
      res_mmt8 = res_mmt8,
      res_fi2 = res_fi2,
      best_models = best_models_export,
      prefix = prefix,
      suffix = suffix
    )
  }

  list(
    ggir_model_data = ggir_model_data,
    models_mmt8 = models_mmt8,
    models_fi2 = models_fi2,
    res_mmt8 = res_mmt8,
    res_fi2 = res_fi2,
    best_models_export = best_models_export,
    export_paths = export_paths
  )
}

show_ggir_focal_results <- function(results) {
  print(results$res_mmt8$radar_plot)
  print(results$res_fi2$radar_plot)
  print(results$res_mmt8$composite_table_export)
  print(results$res_fi2$composite_table_export)
  print(results$best_models_export)

  if (interactive()) {
    utils::View(results$res_mmt8$composite_table_export, title = "GGIR mmt8 composite")
    utils::View(results$res_fi2$composite_table_export, title = "GGIR fi2 composite")
    utils::View(results$best_models_export, title = "GGIR best models")
  }

  invisible(results)
}

# Override the first plotting helpers so the radar closes correctly in coord_polar()
# and previews stay inside the standard RStudio console/Plots workflow.
plot_radar_focal_metrics <- function(score_tbl) {
  radar_metrics <- attr(score_tbl, "radar_metrics")
  family <- unique(score_tbl$family)
  n_metrics <- length(radar_metrics)

  pretty_labels <- c(
    score_beta = "|beta| std.",
    score_ci = "CI precision",
    score_semi_partial_r2 = "semi-partial R2",
    score_predictor_r2 = "partial R2",
    score_r2_marginal = "R2 marginal",
    score_r2_conditional = "R2 conditional",
    score_aic = "AIC (reverse)"
  )

  plot_tbl <- score_tbl |>
    select(parameter, all_of(radar_metrics)) |>
    pivot_longer(
      cols = all_of(radar_metrics),
      names_to = "metric",
      values_to = "score"
    ) |>
    mutate(metric_id = match(metric, radar_metrics))

  plot_tbl_closed <- plot_tbl |>
    group_by(parameter) |>
    arrange(metric_id, .by_group = TRUE) |>
    group_modify(~ bind_rows(
      .x,
      mutate(.x[1, , drop = FALSE], metric_id = n_metrics + 1)
    )) |>
    ungroup()

  ggplot(
    plot_tbl_closed,
    aes(
      x = metric_id,
      y = score,
      group = parameter,
      colour = parameter,
      fill = parameter
    )
  ) +
    geom_polygon(alpha = 0.12, linewidth = 0.9) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    coord_polar(theta = "x", clip = "off") +
    scale_x_continuous(
      breaks = seq_len(n_metrics),
      labels = pretty_labels[radar_metrics],
      limits = c(1, n_metrics + 1)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0.25, 0.50, 0.75, 1.00)
    ) +
    labs(
      title = ifelse(
        family == "lm",
        "Radar comparison of focal predictors (lm)",
        "Radar comparison of focal predictors (lmer)"
      ),
      subtitle = "Legenda obsahuje pouze rozdilne parametry",
      x = NULL,
      y = NULL,
      colour = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "right"
    )
}

show_ggir_focal_results <- function(results) {
  print(results$res_mmt8$radar_plot)
  print(results$res_fi2$radar_plot)
  print(results$res_mmt8$composite_table_export)
  print(results$res_fi2$composite_table_export)
  print(results$best_models_export)

  invisible(results)
}

build_fmsb_radar <- function(score_tbl) {
  radar_metrics <- attr(score_tbl, "radar_metrics")
  family <- unique(score_tbl$family)
  model_subtitle <- attr(score_tbl, "model_subtitle")

  pretty_labels <- c(
    score_beta = "|beta| std.",
    score_ci = "CI precision",
    score_semi_partial_r2 = "R2 - semi-partial",
    score_predictor_r2 = "partial R2",
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
      family = family,
      parameters = score_tbl$parameter,
      axis_labels = pretty_labels[radar_metrics],
      title = ifelse(
        family == "lm",
        "Radar comparison of focal predictors (lm)",
        "Radar comparison of focal predictors (lmer)"
      ),
      subtitle = if (is.null(model_subtitle)) "tested variable model" else model_subtitle
    ),
    class = "fmsb_radar"
  )
}

print.fmsb_radar <- function(x, ...) {
  title_cex <- 60 / 12
  subtitle_cex <- 50 / 12
  axis_cex <- 32 / 12
  y_axis_cex <- 36 / 12
  legend_cex <- 28 / 12
  point_cex <- 36 / 12

  cols <- get_fixed_colors_vis(x$parameters)
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
    legend = x$parameters,
    bty = "n",
    pch = 16,
    col = cols,
    pt.cex = point_cex,
    cex = legend_cex,
    text.font = 2
  )

  invisible(x)
}

plot_radar_focal_metrics <- function(score_tbl) {
  build_fmsb_radar(score_tbl)
}

export_focal_comparison <- function(
    result,
    xlsx_file = file.path(
      "output",
      "tables",
      paste0(format(Sys.Date(), "%y%m%d"), "_radar_comparison_01.xlsx")
    ),
    plot_file = file.path(
      "output",
      "figures",
      paste0(format(Sys.Date(), "%y%m%d"), "_radar_comparison_01.png")
    ),
    width = 14,
    height = 8,
    dpi = 320
) {
  dir.create(dirname(xlsx_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(plot_file), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(xlsx_file)) {
    unlink(xlsx_file, force = TRUE)
  }

  if (file.exists(plot_file)) {
    unlink(plot_file, force = TRUE)
  }

  export_list <- list(
    composite_table = result$composite_table_export,
    raw_metrics = result$raw_metrics |>
      mutate(across(where(is.numeric), ~ round(.x, 4))),
    validation = result$validation_table
  )

  rio::export(export_list, xlsx_file)

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

# ************************************************************
# 11) Run block for this project ----
# ************************************************************

# Expected usage:
# source("scripts/functions/OBJ_01.R")
# source("scripts/Script_myo_2025_iim_clinics_ggir_vis_01.R")
#
# Objects available for preview after sourcing this script:
# ggir_focal_results$res_mmt8$radar_plot
# ggir_focal_results$res_mmt8$composite_table_export
# ggir_focal_results$res_fi2$radar_plot
# ggir_focal_results$res_fi2$composite_table_export
# ggir_focal_results$best_models_export
#
# To disable automatic export during sourcing, define beforehand:
# options(myo.ggir.focal.export = FALSE)

if (exists("d18_ggir_data")) {
  ggir_focal_results <- run_ggir_focal_workflow(
    data = d18_ggir_data,
    export_results = isTRUE(getOption("myo.ggir.focal.export", TRUE))
  )

  res_mmt8 <- ggir_focal_results$res_mmt8
  res_fi2 <- ggir_focal_results$res_fi2
  best_models_export <- ggir_focal_results$best_models_export

  show_ggir_focal_results(ggir_focal_results)

  if (!is.null(ggir_focal_results$export_paths)) {
    message("Exported mmt8 table: ", ggir_focal_results$export_paths$mmt8$xlsx_file)
    message("Exported mmt8 plot: ", ggir_focal_results$export_paths$mmt8$plot_file)
    message("Exported fi2 table: ", ggir_focal_results$export_paths$fi2$xlsx_file)
    message("Exported fi2 plot: ", ggir_focal_results$export_paths$fi2$plot_file)
    message("Exported best-models table: ", ggir_focal_results$export_paths$best_models_file)
  }
}
