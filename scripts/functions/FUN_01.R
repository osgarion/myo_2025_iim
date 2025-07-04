# Libraries & functions ----
## libraries ----
pacman::p_load(update = T,  
               equatiomatic, beepr, tictoc, 
               tidyverse, purrr, furrr, easystats, rio, janitor, ggthemes, car,
               gtsummary, skimr, sjPlot, flextable, ggpubr, rstatix, tidymodels,
               kableExtra, skimr, GGally, testthat, factoextra, gplots, uwot,
               lmerTest, dlookr
)

# Missing values, multivariate analyses
pacman::p_load(naniar, MVN) 

## function specification ----
conflicted::conflicts_prefer(
  janitor::remove_empty,
  dplyr::filter,
  dplyr::mutate,
  dplyr::rename,
  dplyr::summarize,
  dplyr::summarise,
  dplyr::select,
  purrr::map,
  tidyr::extract,
  janitor::clean_names,
  dplyr::relocate,
  lmerTest::lmer,
  GLMMadaptive::predict
)

## Conflicted functions ----
ls_conflicted_01 <- conflicted::conflict_scout()

# Options ----
## Furrr
furrr_options(seed = TRUE,
              scheduling = Inf)
options(knitr.kable.NA = '') # empty space in cells with NAs in kable

## Markdown ----
options(knitr.kable.NA = '',     # empty space in cells with NAs in kable
        scipen = 999)            # non-academic format of numbers




## Number of cores ----
ncores <- availableCores() -1

## Okabe & Ito palette - colorblind palette ----
pal_okabe_ito <- colorblind_pal()(8)[2:8] # ggthemes

# Figures ----
## Group description
# group_colors <- c(C = "#1b9e77", SE = "#d95f02")
# group_fills  <- c(C = "#a6dba0", SE = "#fc8d62")  
# group_shapes <- c(C = 21, SE = 22)               
# group_names  <- c(C = "Control", SE = "Status Epilepticus")



# any ggplot
median_cl_boot <- function(x, conf = 0.95) {
  lconf <- (1 - conf)/2
  uconf <- 1 - lconf
  require(boot)
  bmedian <- function(x, ind) median(x[ind])
  bt <- boot(x, bmedian, 50000)
  bb <- boot.ci(bt, type = "perc")
  data.frame(y = median(x), ymin = quantile(bt$t, lconf), ymax = quantile(bt$t, 
                                                                          uconf))
}

# GGally
my_fn <- function(data, mapping, method="loess", ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(alpha=0.3) + 
    geom_smooth(method=method, ...)
  p
}


# UMAP
kabsch <- function(pm, qm) {
  pm_dims <- dim(pm)
  if (!all(dim(qm) == pm_dims)) {
    stop(call. = TRUE, "Point sets must have the same dimensions")
  }
  # The rotation matrix will have (ncol - 1) leading ones in the diagonal
  diag_ones <- rep(1, pm_dims[2] - 1)
  
  # center the points
  pm <- scale(pm, center = TRUE, scale = FALSE)
  qm <- scale(qm, center = TRUE, scale = FALSE)
  
  am <- crossprod(pm, qm)
  
  svd_res <- svd(am)
  # use the sign of the determinant to ensure a right-hand coordinate system
  d <- determinant(tcrossprod(svd_res$v, svd_res$u))$sign
  dm <- diag(c(diag_ones, d))
  
  # rotation matrix
  um <- svd_res$v %*% tcrossprod(dm, svd_res$u)
  
  # Rotate and then translate to the original centroid location of qm
  sweep(t(tcrossprod(um, pm)), 2, -attr(qm, "scaled:center"))
}

plot_umap <- function(coords, col, pca, main = NULL) {
  plot(kabsch(coords, pca), col = col, xlab = "", ylab = "", 
       main = main)
  legend("topright", legend = unique(col))
  }


fig_mixMODmmt <- function(data, mod, var_indep_name, var_dep_name, upper) {
  if (is.null(mod)) return(NULL)  # 
  
  data_set <- data |> 
    mutate(
      var_indep_value_scl = log(var_indep_value +1),
      y_capped = pmin(var_dep_value, upper),
      ind = as.integer(var_dep_value >= upper)
    ) |> 
    drop_na(var_dep_value, var_indep_value) |>
    group_by(projekt_id) |>
    filter(n() >= 3) |>
    ungroup()
  
  pred_df <- GLMMadaptive::effectPlotData(mod, newdata = na.omit(data_set))
  
  fig <- ggplot(pred_df, aes(x = var_indep_value, y = var_dep_value)) +
    geom_point(aes(color = poradie_vysetrenia), alpha = 0.6, show.legend = FALSE) +
    geom_line(aes(y = pred, color = poradie_vysetrenia), size = 1, show.legend = FALSE) +
    geom_ribbon(aes(ymax = upp, ymin = low, fill = poradie_vysetrenia), 
                alpha = 0.3, linetype = 0) +
    facet_wrap(~ poradie_vysetrenia, scales = "free") +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    labs(
      x = var_indep_name,
      y = var_dep_name,
      fill = "Time point"
    ) +
    theme_sjplot2() +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "gray90", color = NA)
    )
  
  return(fig)
}


# tables ----
## column with 0, 1, and NA ----
is_01_col <- function(x) {
  all(unique(x) %in% c(0, 1, NA))
}

## skimr ----
my_skim <- skimr::skim_with(numeric = sfl(median = ~ median(., na.rm = TRUE),
                                          mad = ~ mad(., na.rm = TRUE)), 
                            append = T)

## import ----
file_to_load <- function(folder_path, phrase) {
  # Získání všech souborů (včetně cesty)
  files <- list.files(folder_path, recursive = TRUE, full.names = TRUE)
  
  # Filtrování souborů podle fráze v názvu
  matching_files <- files[grepl(phrase, basename(files), ignore.case = TRUE)]
  
  return(matching_files)
}

## datafor mix_model ----
prepare_for_censored_model <- function(data, dep_var, indep_var, group_var, upper = NULL) {
  # Základní kontrola vstupů
  if (!all(c(dep_var, indep_var, group_var) %in% names(data))) {
    stop("Zadané názvy proměnných neodpovídají sloupcům v datech.")
  }
  
  # Horní mez pro y_capped (implicitně maximum závislé proměnné)
  if (is.null(upper)) {
    upper <- max(data[[dep_var]], na.rm = TRUE)
  }
  
  # Příprava dat
  data_prepared <- data |>
    mutate(
      var_dep_value        = .data[[dep_var]],
      var_indep_value      = .data[[indep_var]],
      projekt_id           = .data[[group_var]],
      var_indep_value_log  = log(var_indep_value + 1),
      var_indep_value_scl  = as.numeric(scale(var_indep_value_log)),
      y_capped             = pmin(var_dep_value, upper),
      ind                  = as.integer(var_dep_value >= upper)
    ) |>
    drop_na(var_dep_value, var_indep_value, var_indep_value_scl, y_capped, ind) |>
    group_by(projekt_id) |>
    filter(n() >= 3) |>
    mutate(
      sd_y = sd(y_capped),
      sd_x = sd(var_indep_value_scl)
    ) |>
    filter(sd_y > 0, sd_x > 0) |>
    ungroup() |>
    select(-sd_y, -sd_x)
  
  return(data_prepared)
}

mod_mixMODmmt <- function(data, upper) {
  
  
  data_set <- data |>
    mutate(
      var_indep_value_scl = log(var_indep_value + 1),
      y_capped = pmin(var_dep_value, upper),
      ind = as.integer(var_dep_value >= upper)
    ) |>
    drop_na(var_dep_value, var_indep_value) |>
    group_by(projekt_id) |>
    filter(n() >= 3) |>
    mutate(
      sd_y = sd(y_capped),
      sd_x = sd(var_indep_value_scl)
    ) |>
    filter(sd_y > 0, sd_x > 0) |>  
    ungroup() |>
    select(-sd_y, -sd_x)  
  
  fm <- mixed_model(
    fixed  = cbind(y_capped, ind) ~ poradie_vysetrenia + var_indep_value_scl,
    random = ~ 1 | projekt_id,
    family = censored.normal(),
    data   = data_set
  )
  return(fm)
}

## Possibly ----
possLMER <- possibly(.f=lmer, otherwise = NA_real_)
possmixMODmmt <- possibly(mod_mixMODmmt, otherwise = NULL)
