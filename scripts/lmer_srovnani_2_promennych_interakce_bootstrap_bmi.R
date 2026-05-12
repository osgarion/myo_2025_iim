# ============================================================
# GGIR bootstrap lmer
# Model:
# ig_gradient_enmo ~ exam + age + sex + bcm_between + bmi_between
#
# Bootstrap:
# - n = 5000
# - v každém cyklu náhodný výběr 50 pacientů
# - defaultně s návratem = klasický bootstrap
# - celý preprocessing se opakuje v každém cyklu
# - paralelizace přes future/furrr
# ============================================================

# ************************************************************
# 0) Balíčky ----
# ************************************************************
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(furrr)
  library(future)
  library(progressr)
  library(lme4)
  library(datawizard)
  library(parameters)
  library(performance)
  library(rio)
  library(ggplot2)
  library(stringr)
})

# ************************************************************
# 0A) CONFIG ----
# ************************************************************
cfg <- list(
  
  # ------------------------------
  # cesty
  # ------------------------------
  data_path = "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/GENEActivII/rCCA_rozdelene_sady_parametru/ggir_vs_myokiny_bez_imputace/20263003.xlsx",
  
  output_dir = "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/GENEActivII/lmer_custom_multivariable_phase_angle/residualizace/bootstrap",
  
  import_which = "sheet1",
  
  # ------------------------------
  # proměnné
  # ------------------------------
  id_var       = "filename_2",
  response_var = "ig_gradient_enmo",
  time_var_raw = "m",
  time_var     = "exam",
  sex_var      = "sex",
  
  continuous_vars = c(
    "age",
    "bcm",
    "bmi"
  ),
  
  # ------------------------------
  # bootstrap nastavení
  # ------------------------------
  n_boot = 5000,
  n_patients_per_boot = 50,
  
  # TRUE = klasický bootstrap pacientů s návratem
  # FALSE = náhodný výběr 50 unikátních pacientů bez návratu
  sample_replace = TRUE,
  
  seed = 1234,
  
  # ------------------------------
  # paralelizace
  # ------------------------------
  n_workers = max(1, parallel::detectCores() - 1),
  
  # ------------------------------
  # model
  # ------------------------------
  model_rhs = "exam + age + sex + bcm_between + bmi_between",
  
  # volitelně fitovat i standardizovaný model
  make_standardized_model = TRUE,
  z_standardize_response = TRUE,
  
  # ------------------------------
  # výstupy
  # ------------------------------
  dpi = 300,
  plot_width = 10,
  plot_height = 7
)

# ************************************************************
# 0B) Kontroly ----
# ************************************************************
if (!file.exists(cfg$data_path)) {
  stop("Vstupní soubor neexistuje: ", cfg$data_path)
}

if (!dir.exists(cfg$output_dir)) {
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(cfg$output_dir)) {
  stop("Výstupní složku se nepodařilo vytvořit: ", cfg$output_dir)
}

message("Vstupní data: ", normalizePath(cfg$data_path, winslash = "/"))
message("Výstupy budou uloženy do: ", normalizePath(cfg$output_dir, winslash = "/"))
message("Počet bootstrap cyklů: ", cfg$n_boot)
message("Počet pacientů v jednom cyklu: ", cfg$n_patients_per_boot)
message("Vzorkování s návratem: ", cfg$sample_replace)
message("Počet paralelních workerů: ", cfg$n_workers)

# ************************************************************
# 0C) Pomocné funkce ----
# ************************************************************

is_01_col <- function(x) {
  ux <- unique(stats::na.omit(x))
  all(ux %in% c(0, 1))
}

coerce_model_vars <- function(data, exam_var = "exam") {
  
  if (exam_var %in% names(data)) {
    if (is.numeric(data[[exam_var]]) || is.integer(data[[exam_var]])) {
      data[[exam_var]] <- factor(
        data[[exam_var]],
        levels = c(0, 3, 6, 18, 30, 48),
        labels = c("M0", "M3", "M6", "M18", "M30", "M48")
      )
    } else {
      data[[exam_var]] <- factor(
        data[[exam_var]],
        levels = c("M0", "M3", "M6", "M18", "M30", "M48")
      )
    }
  }
  
  candidate_vars <- names(data)[vapply(data, is.numeric, logical(1))]
  binary_vars <- candidate_vars[vapply(data[candidate_vars], is_01_col, logical(1))]
  binary_vars <- setdiff(binary_vars, exam_var)
  
  for (v in binary_vars) {
    if (v == "sex") {
      data[[v]] <- factor(
        data[[v]],
        levels = c(0, 1),
        labels = c("Female", "Male")
      )
    } else if (v == "ild") {
      data[[v]] <- factor(
        data[[v]],
        levels = c(0, 1),
        labels = c("No", "Yes")
      )
    } else {
      data[[v]] <- factor(
        data[[v]],
        levels = c(0, 1),
        labels = c("0", "1")
      )
    }
  }
  
  data
}

safe_scale <- function(x) {
  sx <- stats::sd(x, na.rm = TRUE)
  
  if (is.na(sx) || sx == 0) {
    return(rep(NA_real_, length(x)))
  }
  
  as.numeric(scale(x))
}

extract_fixed_params <- function(model, boot_id, model_type = "raw") {
  
  parameters::model_parameters(
    model,
    effects = "fixed",
    ci_random = FALSE
  ) |>
    as.data.frame() |>
    tibble::as_tibble() |>
    dplyr::mutate(
      boot_id = boot_id,
      model_type = model_type,
      converged = TRUE,
      singular = lme4::isSingular(model, tol = 1e-4),
      n_obs = stats::nobs(model),
      n_groups = length(lme4::getME(model, "flist")[[1]]),
      error_message = NA_character_
    )
}

make_boot_sample <- function(dat, eligible_ids, cfg) {
  
  sampled_ids <- sample(
    eligible_ids,
    size = cfg$n_patients_per_boot,
    replace = cfg$sample_replace
  )
  
  boot_key <- tibble::tibble(
    original_id = sampled_ids,
    boot_cluster_number = seq_along(sampled_ids),
    boot_id = paste0(sampled_ids, "__bootcluster_", boot_cluster_number)
  )
  
  dat |>
    dplyr::inner_join(
      boot_key,
      by = setNames("original_id", cfg$id_var),
      relationship = "many-to-many"
    )
}

prepare_boot_data <- function(dat_sample, cfg) {
  
  # Důležité:
  # within-between složky se počítají znovu v každém bootstrap cyklu.
  # Používáme boot_id, protože při vzorkování s návratem může být stejný pacient
  # ve vzorku víckrát. Každá kopie pacienta je pak samostatný bootstrap cluster.
  
  dat_wb <- dat_sample |>
    datawizard::demean(
      select = c("bcm", "bmi"),
      by = "boot_id"
    ) |>
    coerce_model_vars(exam_var = cfg$time_var)
  
  dat_wb |>
    dplyr::filter(
      !is.na(.data[[cfg$response_var]]),
      !is.na(.data[[cfg$time_var]]),
      !is.na(age),
      !is.na(.data[[cfg$sex_var]]),
      !is.na(bcm_between),
      !is.na(bmi_between),
      !is.na(boot_id)
    )
}

fit_one_boot <- function(i, dat, eligible_ids, cfg) {
  
  tryCatch({
    
    dat_sample <- make_boot_sample(
      dat = dat,
      eligible_ids = eligible_ids,
      cfg = cfg
    )
    
    dat_boot <- prepare_boot_data(
      dat_sample = dat_sample,
      cfg = cfg
    )
    
    if (dplyr::n_distinct(dat_boot$boot_id) < 5) {
      stop("Příliš málo bootstrap clusterů po odstranění chybějících hodnot.")
    }
    
    if (nrow(dat_boot) < 10) {
      stop("Příliš málo řádků po odstranění chybějících hodnot.")
    }
    
    formula_raw <- stats::as.formula(
      paste0(
        cfg$response_var,
        " ~ ",
        cfg$model_rhs,
        " + (1 | boot_id)"
      )
    )
    
    mod_raw <- lme4::lmer(
      formula = formula_raw,
      data = dat_boot,
      REML = FALSE,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5),
        calc.derivs = FALSE
      )
    )
    
    out_raw <- extract_fixed_params(
      model = mod_raw,
      boot_id = i,
      model_type = "raw"
    )
    
    if (isTRUE(cfg$make_standardized_model)) {
      
      dat_boot_z <- dat_boot |>
        dplyr::mutate(
          age_z = safe_scale(age),
          bcm_between_z = safe_scale(bcm_between),
          bmi_between_z = safe_scale(bmi_between)
        )
      
      response_z <- cfg$response_var
      
      if (isTRUE(cfg$z_standardize_response)) {
        response_z <- paste0(cfg$response_var, "_z")
        dat_boot_z[[response_z]] <- safe_scale(dat_boot_z[[cfg$response_var]])
      }
      
      dat_boot_z <- dat_boot_z |>
        dplyr::filter(
          !is.na(.data[[response_z]]),
          !is.na(age_z),
          !is.na(bcm_between_z),
          !is.na(bmi_between_z)
        )
      
      formula_z <- stats::as.formula(
        paste0(
          response_z,
          " ~ exam + age_z + sex + bcm_between_z + bmi_between_z + (1 | boot_id)"
        )
      )
      
      mod_z <- lme4::lmer(
        formula = formula_z,
        data = dat_boot_z,
        REML = FALSE,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 2e5),
          calc.derivs = FALSE
        )
      )
      
      out_z <- extract_fixed_params(
        model = mod_z,
        boot_id = i,
        model_type = "standardized"
      )
      
      dplyr::bind_rows(out_raw, out_z)
      
    } else {
      
      out_raw
    }
    
  }, error = function(e) {
    
    tibble::tibble(
      boot_id = i,
      model_type = c("raw", "standardized"),
      Parameter = NA_character_,
      Coefficient = NA_real_,
      SE = NA_real_,
      CI = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      t = NA_real_,
      df_error = NA_real_,
      p = NA_real_,
      Effects = NA_character_,
      converged = FALSE,
      singular = NA,
      n_obs = NA_integer_,
      n_groups = NA_integer_,
      error_message = conditionMessage(e)
    )
  })
}

summarise_boot <- function(boot_params) {
  
  boot_params |>
    dplyr::filter(
      converged,
      !is.na(Parameter),
      !is.na(Coefficient)
    ) |>
    dplyr::group_by(model_type, Parameter) |>
    dplyr::summarise(
      n_success = dplyr::n(),
      estimate_mean = mean(Coefficient, na.rm = TRUE),
      estimate_median = median(Coefficient, na.rm = TRUE),
      boot_se = stats::sd(Coefficient, na.rm = TRUE),
      ci_low_percentile = stats::quantile(Coefficient, probs = 0.025, na.rm = TRUE),
      ci_high_percentile = stats::quantile(Coefficient, probs = 0.975, na.rm = TRUE),
      p_boot_two_sided = 2 * pmin(
        mean(Coefficient <= 0, na.rm = TRUE),
        mean(Coefficient >= 0, na.rm = TRUE)
      ),
      prop_positive = mean(Coefficient > 0, na.rm = TRUE),
      prop_negative = mean(Coefficient < 0, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(model_type, Parameter)
}

make_density_plot <- function(boot_params, model_type_selected = "standardized") {
  
  boot_params |>
    dplyr::filter(
      converged,
      model_type == model_type_selected,
      !is.na(Parameter),
      !is.na(Coefficient),
      Parameter != "(Intercept)"
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(x = Coefficient)
    ) +
    ggplot2::geom_density(linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~ Parameter, scales = "free") +
    ggplot2::labs(
      title = paste0("Bootstrap distributions of fixed effects: ", model_type_selected),
      x = "Bootstrap coefficient",
      y = "Density"
    ) +
    ggplot2::theme_minimal()
}

make_forest_plot <- function(boot_summary, model_type_selected = "standardized") {
  
  boot_summary |>
    dplyr::filter(
      model_type == model_type_selected,
      Parameter != "(Intercept)"
    ) |>
    dplyr::mutate(
      Parameter = factor(Parameter, levels = rev(unique(Parameter)))
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = estimate_median,
        y = Parameter
      )
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = ci_low_percentile,
        xmax = ci_high_percentile
      ),
      height = 0.2
    ) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::labs(
      title = paste0("Bootstrap percentile CI: ", model_type_selected),
      x = "Median bootstrap coefficient with 95% percentile CI",
      y = NULL
    ) +
    ggplot2::theme_minimal()
}

# ************************************************************
# 1) Načtení dat ----
# ************************************************************
if (is.null(cfg$import_which)) {
  dat_raw <- rio::import(cfg$data_path)
} else {
  dat_raw <- rio::import(cfg$data_path, which = cfg$import_which)
}

dat_raw <- tibble::as_tibble(dat_raw)

# přejmenování času: m -> exam
if (cfg$time_var_raw %in% names(dat_raw) && !(cfg$time_var %in% names(dat_raw))) {
  dat_raw <- dat_raw |>
    dplyr::rename(!!cfg$time_var := !!rlang::sym(cfg$time_var_raw))
}

required_vars <- unique(c(
  cfg$id_var,
  cfg$response_var,
  cfg$time_var,
  cfg$sex_var,
  cfg$continuous_vars
))

missing_required <- setdiff(required_vars, names(dat_raw))

if (length(missing_required) > 0) {
  stop(
    "V datech chybí tyto požadované proměnné: ",
    paste(missing_required, collapse = ", ")
  )
}

# základní kontrola typů
dat_raw <- dat_raw |>
  dplyr::filter(!is.na(.data[[cfg$id_var]]))

# ************************************************************
# 2) Určení způsobilých pacientů ----
# ************************************************************
# Způsobilý pacient = po vytvoření bcm_between a bmi_between má alespoň jeden
# použitelný řádek pro model.

dat_check <- dat_raw |>
  datawizard::demean(
    select = c("bcm", "bmi"),
    by = cfg$id_var
  ) |>
  coerce_model_vars(exam_var = cfg$time_var) |>
  dplyr::filter(
    !is.na(.data[[cfg$response_var]]),
    !is.na(.data[[cfg$time_var]]),
    !is.na(age),
    !is.na(.data[[cfg$sex_var]]),
    !is.na(bcm_between),
    !is.na(bmi_between)
  )

eligible_ids <- dat_check |>
  dplyr::distinct(!!rlang::sym(cfg$id_var)) |>
  dplyr::pull(!!rlang::sym(cfg$id_var))

if (length(eligible_ids) < cfg$n_patients_per_boot && !isTRUE(cfg$sample_replace)) {
  stop(
    "Pro výběr bez návratu není dost způsobilých pacientů. Počet způsobilých pacientů: ",
    length(eligible_ids),
    "; požadováno: ",
    cfg$n_patients_per_boot
  )
}

message("Počet způsobilých pacientů: ", length(eligible_ids))

# ************************************************************
# 3) Full-sample model pro orientaci ----
# ************************************************************
dat_full <- dat_raw |>
  dplyr::mutate(boot_id = .data[[cfg$id_var]]) |>
  prepare_boot_data(cfg = cfg)

full_formula <- stats::as.formula(
  paste0(
    cfg$response_var,
    " ~ ",
    cfg$model_rhs,
    " + (1 | boot_id)"
  )
)

full_model <- lme4::lmer(
  formula = full_formula,
  data = dat_full,
  REML = FALSE,
  control = lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

full_model_params <- parameters::model_parameters(
  full_model,
  effects = "fixed",
  ci_random = FALSE
) |>
  as.data.frame() |>
  tibble::as_tibble()

full_collinearity <- performance::check_collinearity(full_model) |>
  as.data.frame() |>
  tibble::as_tibble()

saveRDS(
  full_model,
  file.path(cfg$output_dir, "full_sample_model_raw.rds")
)

# ************************************************************
# 4) Paralelní bootstrap ----
# ************************************************************
set.seed(cfg$seed)

future::plan(
  future::multisession,
  workers = cfg$n_workers
)

options(
  future.globals.maxSize = 4 * 1024^3
)

progressr::handlers(global = TRUE)

message("Spouštím bootstrap...")

boot_params <- progressr::with_progress({
  
  p <- progressr::progressor(along = seq_len(cfg$n_boot))
  
  furrr::future_map_dfr(
    .x = seq_len(cfg$n_boot),
    .f = function(i) {
      p()
      fit_one_boot(
        i = i,
        dat = dat_raw,
        eligible_ids = eligible_ids,
        cfg = cfg
      )
    },
    .options = furrr::furrr_options(
      seed = cfg$seed,
      packages = c(
        "dplyr",
        "tidyr",
        "tibble",
        "lme4",
        "datawizard",
        "parameters",
        "performance",
        "rlang"
      )
    )
  )
})

future::plan(future::sequential)

message("Bootstrap dokončen.")

# ************************************************************
# 5) Souhrny bootstrapu ----
# ************************************************************
boot_summary <- summarise_boot(boot_params)

boot_status <- boot_params |>
  dplyr::distinct(
    boot_id,
    model_type,
    converged,
    singular,
    n_obs,
    n_groups,
    error_message
  ) |>
  dplyr::group_by(model_type) |>
  dplyr::summarise(
    n_total = dplyr::n(),
    n_converged = sum(converged, na.rm = TRUE),
    n_failed = sum(!converged, na.rm = TRUE),
    n_singular = sum(singular %in% TRUE, na.rm = TRUE),
    median_n_obs = median(n_obs, na.rm = TRUE),
    median_n_groups = median(n_groups, na.rm = TRUE),
    .groups = "drop"
  )

failed_boots <- boot_params |>
  dplyr::filter(!converged) |>
  dplyr::distinct(
    boot_id,
    model_type,
    error_message
  )

# ************************************************************
# 6) Grafy ----
# ************************************************************
p_density_raw <- make_density_plot(
  boot_params = boot_params,
  model_type_selected = "raw"
)

ggplot2::ggsave(
  filename = file.path(cfg$output_dir, "bootstrap_density_raw.png"),
  plot = p_density_raw,
  width = cfg$plot_width,
  height = cfg$plot_height,
  dpi = cfg$dpi
)

p_forest_raw <- make_forest_plot(
  boot_summary = boot_summary,
  model_type_selected = "raw"
)

ggplot2::ggsave(
  filename = file.path(cfg$output_dir, "bootstrap_forest_raw.png"),
  plot = p_forest_raw,
  width = cfg$plot_width,
  height = cfg$plot_height,
  dpi = cfg$dpi
)

if (isTRUE(cfg$make_standardized_model)) {
  
  p_density_z <- make_density_plot(
    boot_params = boot_params,
    model_type_selected = "standardized"
  )
  
  ggplot2::ggsave(
    filename = file.path(cfg$output_dir, "bootstrap_density_standardized.png"),
    plot = p_density_z,
    width = cfg$plot_width,
    height = cfg$plot_height,
    dpi = cfg$dpi
  )
  
  p_forest_z <- make_forest_plot(
    boot_summary = boot_summary,
    model_type_selected = "standardized"
  )
  
  ggplot2::ggsave(
    filename = file.path(cfg$output_dir, "bootstrap_forest_standardized.png"),
    plot = p_forest_z,
    width = cfg$plot_width,
    height = cfg$plot_height,
    dpi = cfg$dpi
  )
}

# ************************************************************
# 7) Export tabulek ----
# ************************************************************
export_list <- list(
  config = tibble::enframe(cfg),
  full_sample_parameters_raw = full_model_params,
  full_sample_collinearity_raw = full_collinearity,
  bootstrap_summary = boot_summary,
  bootstrap_status = boot_status,
  bootstrap_all_parameters = boot_params,
  failed_boots = failed_boots
)

rio::export(
  export_list,
  file.path(cfg$output_dir, "bootstrap_lmer_report.xlsx")
)

readr::write_csv(
  boot_summary,
  file.path(cfg$output_dir, "bootstrap_summary.csv")
)

readr::write_csv(
  boot_params,
  file.path(cfg$output_dir, "bootstrap_all_parameters.csv")
)

readr::write_csv(
  boot_status,
  file.path(cfg$output_dir, "bootstrap_status.csv")
)

if (nrow(failed_boots) > 0) {
  readr::write_csv(
    failed_boots,
    file.path(cfg$output_dir, "failed_boots.csv")
  )
}

# ************************************************************
# 8) Uložení RDS ----
# ************************************************************
saveRDS(
  boot_params,
  file.path(cfg$output_dir, "bootstrap_all_parameters.rds")
)

saveRDS(
  boot_summary,
  file.path(cfg$output_dir, "bootstrap_summary.rds")
)

saveRDS(
  boot_status,
  file.path(cfg$output_dir, "bootstrap_status.rds")
)

# ************************************************************
# 9) Závěrečná zpráva ----
# ************************************************************
message("Hotovo.")
message("Výstupy uloženy do: ", normalizePath(cfg$output_dir, winslash = "/"))
message("Hlavní Excel report: ", file.path(cfg$output_dir, "bootstrap_lmer_report.xlsx"))
