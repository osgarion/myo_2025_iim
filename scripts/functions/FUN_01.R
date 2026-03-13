# Libraries & functions ----
## libraries ----
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(recipes::update)
}
pacman_update_01 <- getOption("myo_pacman_update", TRUE)
pacman_install_01 <- getOption("myo_pacman_install", TRUE)

pacman::p_load(update = pacman_update_01, install = pacman_install_01,
               equatiomatic, beepr, tictoc, 
               tidyverse, purrr, furrr, easystats, rio, janitor, ggthemes, car,
               gtsummary, skimr, sjPlot, flextable, ggpubr, rstatix, tidymodels,
               kableExtra, skimr, GGally, testthat, factoextra, gplots, uwot,
               lmerTest, dlookr, multilevelmod,furrr,ggforce, lazyWeave, paletteer,
               emmeans, openxlsx, GLMMadaptive, FactoMineR, cluster, mixOmics,
               mice, MOFA2, patchwork
)

# Missing values, multivariate analyses
pacman::p_load(
  naniar, MVN, visdat,
  update = pacman_update_01,
  install = pacman_install_01
)

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
  dplyr::recode,
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
ncores_default_01 <- max(1L, min(4L, future::availableCores() - 1L))
ncores <- as.integer(getOption("myo_ncores", ncores_default_01))
if (!is.finite(ncores) || is.na(ncores)) ncores <- ncores_default_01
ncores <- max(1L, ncores)

## Okabe & Ito palette - colorblind palette ----
pal_okabe_ito <- colorblind_pal()(8)[2:8] # ggthemes

# Figures ----
## Group description
# group_colors <- c(C = "#1b9e77", SE = "#d95f02")
# group_fills  <- c(C = "#a6dba0", SE = "#fc8d62")  
# group_shapes <- c(C = 21, SE = 22)               
# group_names  <- c(C = "Control", SE = "Status Epilepticus")



# any ggplot
median_cl_boot <- function(x, conf = 0.95, B = 50000) {
  lconf <- (1 - conf)/2
  uconf <- 1 - lconf
  require(boot)
  bmedian <- function(x, ind) median(x[ind])
  bt <- boot(x, bmedian, B)
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
## pca + umap + dbscan
plot_dbscan_pca_umap <- function(data_numeric, groups, eps = 0.5, minPts = 5) {
  # Checks
  stopifnot(nrow(data_numeric) == length(groups))
  
  # Packages
  require(ggplot2)
  require(dbscan)
  require(umap)
  require(gridExtra)
  
  # Run DBSCAN
  db <- dbscan::dbscan(data_numeric, eps = eps, minPts = minPts)
  cluster <- as.factor(db$cluster)  # 0 means noise
  
  group_factor <- as.factor(groups)
  
  # PCA
  pca <- prcomp(data_numeric, scale. = TRUE)
  pca_df <- as.data.frame(pca$x[, 1:2])
  colnames(pca_df) <- c("PC1", "PC2")
  pca_df$Cluster <- cluster
  pca_df$Group <- group_factor
  
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Cluster)) +
    geom_point(size = 3, alpha = 0.3) +
    theme_minimal() +
    labs(title = "PCA + DBSCAN", color = "Group", shape = "DBSCAN Cluster")
  
  # UMAP
  umap_result <- umap::umap(data_numeric)
  umap_df <- as.data.frame(umap_result$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$Cluster <- cluster
  umap_df$Group <- group_factor
  
  p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group, shape = Cluster)) +
    geom_point(size = 3, alpha = 0.3) +
    theme_minimal() +
    labs(title = "UMAP + DBSCAN", color = "Group", shape = "DBSCAN Cluster")
  
  # Combine
  gridExtra::grid.arrange(p_pca, p_umap, nrow = 1)
}


# plot mixed effect
plot_mixed_terms_01 <- function(model, 
                                var_pattern = "var_indep_",
                                xlab = NULL, 
                                ylab = NULL,
                                p_label_size = 6,
                                ...) {
  library(broom.mixed)
  library(sjPlot)
  library(ggplot2)
  library(sjmisc)    # for theme_sjplot2()
  
  # 1) grab fixed‐effect terms + p‐values
  fe <- broom.mixed::tidy(as_lmerModLmerTest(model$fit), effects = "fixed")
  
  # 2) select those matching your pattern
  sel_terms <- grep(var_pattern, fe$term, value = TRUE)
  if (length(sel_terms) == 0) {
    stop("No fixed-effect terms matched pattern: ", var_pattern)
  }
  if (length(sel_terms) > 1) {
    warning("More than one term matched; plotting all. p-value will use the first.")
  }
  
  # 3) extract the p-value of the first matched term
  pval <- fe %>% 
    filter(term == sel_terms[1]) %>% 
    pull(p.value)
  
  # 4) build the emmeans plot
  p <- sjPlot::plot_model(
    model,
    type      = "emm",
    terms     = sel_terms,
    show.data = TRUE
  )
  
  # 5) apply axis labels if given
  if (!is.null(xlab) || !is.null(ylab)) {
    p <- p + labs(
      x = xlab %||% waiver(),
      y = ylab %||% waiver()
    )
  }
  
  # 6) annotate p-value in top-right
  p <- p +
    annotate(
      "text",
      x     = Inf,
      y     = Inf,
      label = paste0("p = ", signif(pval, 2)),
      hjust = 1.1,
      vjust = 1.5,
      size  = p_label_size
    ) +
    
    # 7) custom theme
    theme_sjplot2() +
    theme(
      plot.title       = element_blank(),
      plot.subtitle    = element_text(size = 20, color = "grey40", hjust = 0),
      axis.title       = element_text(face = "bold", size = 20),
      axis.text        = element_text(size = 18, color = "grey20"),
      panel.grid.major = element_line(linetype = "dotted", color = "grey80"),
      panel.grid.minor = element_blank(),
      legend.position  = "none"
    )
  
  return(p)
}

plot_mixed_terms_02 <- function(model, 
                                var_pattern = "var_indep_", 
                                xlab        = NULL, 
                                ylab        = NULL) {
  # model        : a parsnip model_fit (with $fit = lmerMod) or a bare lmerMod
  # var_pattern  : regex/substring to match raw vs scaled predictor
  # xlab, ylab   : axis labels; NULL will default to the matched term / "Predicted response"
  
  library(ggeffects)
  library(ggplot2)
  library(broom.mixed)
  library(lmerTest)
  
  # 1) pull out the lmerMod
  engine <- if (inherits(model, "model_fit")) model$fit else model
  
  # 2) find the exact fixed‐effect term(s)
  fe      <- broom.mixed::tidy(as_lmerModLmerTest(engine), effects = "fixed")
  matches <- grep(var_pattern, fe$term, value = TRUE)
  if (length(matches) == 0) {
    stop("No fixed-effect term matches: ", var_pattern)
  }
  # prefer the scaled variant if present:
  predictor <- if ("var_indep_value_scl" %in% matches) {
    "var_indep_value_scl"
  } else matches[1]
  
  # 3) extract its p-value
  pval <- fe %>% 
    filter(term == predictor) %>% 
    pull(p.value) %>% 
    signif(2)
  
  # 4) get ggpredict data
  preds <- ggpredict(
    engine,
    terms = c(predictor, "podtyp_nemoci_zjednoduseny")
  )
  
  # 5) build the base plot
  p <- plot(preds, show_data = TRUE) +
    facet_wrap(~ group, scales = "free_y") +
    labs(
      title  = paste0("p-value = ", pval),
      x      = xlab %||% predictor,
      y      = ylab %||% "Predicted response",
      colour = "Subtype"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text       = element_text(face = "bold", size = 12),
      axis.title       = element_text(face = "bold", size = 14),
      axis.text        = element_text(size = 12),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# podtyp nemoci + poradi vysetreni
plot_mixed_terms_03 <- function(model, 
                                var_pattern = "var_indep_",
                                var_indep2 = "podtyp_nemoci_zjednoduseny",
                                xlab        = NULL, 
                                ylab        = NULL) {
  # model        : a parsnip model_fit (with $fit = lmerMod) or a bare lmerMod
  # var_pattern  : regex/substring to match raw vs scaled predictor
  # xlab, ylab   : axis labels; NULL will default to the matched term / "Predicted response"
  
  library(ggeffects)
  library(ggplot2)
  library(broom.mixed)
  library(lmerTest)
  
  # 1) pull out the lmerMod
  engine <- if (inherits(model, "model_fit")) model$fit else model
  
  # 2) find the exact fixed‐effect term(s)
  fe      <- broom.mixed::tidy(as_lmerModLmerTest(engine), effects = "fixed")
  matches <- grep(var_pattern, fe$term, value = TRUE)
  if (length(matches) == 0) {
    stop("No fixed-effect term matches: ", var_pattern)
  }
  # prefer the scaled variant if present:
  predictor <- if ("var_indep_value_scl" %in% matches) {
    "var_indep_value_scl"
  } else matches[1]
  
  # 3) extract its p-value
  pval <- fe %>% 
    filter(term == predictor) %>% 
    pull(p.value) %>% 
    signif(2)
  
  
  # 4) get ggpredict data
  preds <- ggpredict(
    engine,
    terms = c(predictor, var_indep2)
  )
  
  # 5) build the base plot
  fml <- as.formula(paste0("~", var_indep2))
  
  ct <- emtrends(engine,
                 fml,
                 var =predictor) |>
    summary(infer = c(TRUE, TRUE)) |> 
    data.frame()
  
  preds2 <- ggpredict(
    engine,
    terms = c(predictor, var_indep2, "poradie_vysetrenia")
  )
  
  raw_data <- attr(preds2, "rawdata")
  
  # 1) Build a lookup table of p‐values per facet
  pval_df <- ct %>%
    mutate(
      group   = as.character(.data[[var_indep2]]),         # vezmeme hodnoty ve sloupci "pohlavi"
      p.label = paste0("p = ", signif(p.value, 3))         # naformátujeme p‑hodnotu
    ) %>%
    select(group, p.label)
  
  
  colors_20 <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000",
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
    "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
    "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
  )
  
  # 2) Your base plot
  p <- plot(preds) +
    geom_point(
      data        = raw_data,
      aes(x = x, y = response, colour = facet, fill = facet, shape = facet),
      size = 3,
      alpha       = 0.3,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(values = colors_20) +
    scale_fill_manual(values = colors_20) +
    scale_shape_manual(values = c(21, 22, 23, 25,21, 22, 23, 25), name = "Examination") +
    # paletteer::scale_colour_paletteer_d("tvthemes::Steven") +
    # paletteer::scale_fill_paletteer_d("tvthemes::Steven") +
    facet_wrap(~ group, scales = "free_y") +
    labs(
      title  = paste0("p-value = ", pval),
      x      = xlab %||% predictor,
      y      = ylab %||% "Predicted response",
      colour = ""
    ) +
    theme_minimal(base_size = 18) +
    theme(
      strip.text       = element_text(face = "bold", size = 22),
      axis.title       = element_text(face = "bold", size = 22),
      axis.text        = element_text(size = 21),
      plot.title       = element_text(face = "bold", size = 25, hjust = 0.5),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
  
  # 3) Add a geom_text layer to annotate each facet
  p <- p + 
    geom_text(
      data = pval_df,
      aes(x = Inf, y = Inf, label = p.label),
      hjust = 1.1,    # nudge left from right edge
      vjust = 1.1,    # nudge down from top
      size  = 6, 
      inherit.aes = FALSE
    )
  
  return(p)
  
}

# podtyp nemoci + odpoved po 6 mesicich
plot_mixed_terms_04 <- function(model, 
                                var_pattern = "var_indep_", 
                                var_indep2 = "podtyp_nemoci_zjednoduseny",
                                xlab        = NULL, 
                                ylab        = NULL,
                                data        = data) {
  # model        : a parsnip model_fit (with $fit = lmerMod) or a bare lmerMod
  # var_pattern  : regex/substring to match raw vs scaled predictor
  # xlab, ylab   : axis labels; NULL will default to the matched term / "Predicted response"
  
  library(ggeffects)
  library(ggplot2)
  library(broom.mixed)
  library(lmerTest)
  
  colors_20 <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000",
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
    "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
    "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
  )
  
  
  # 1) pull out the lmerMod
  engine <- if (inherits(model, "model_fit")) model$fit else model
  
  # 2) find the exact fixed‐effect term(s)
  fe      <- broom.mixed::tidy(as_lmerModLmerTest(engine), effects = "fixed")
  matches <- grep(var_pattern, fe$term, value = TRUE)
  if (length(matches) == 0) {
    stop("No fixed-effect term matches: ", var_pattern)
  }
  # prefer the scaled variant if present:
  predictor <- if ("var_indep_value_scl" %in% matches) {
    "var_indep_value_scl"
  } else matches[1]
  
  # 3) extract its p-value
  pval <- fe %>% 
    filter(term == predictor) %>% 
    pull(p.value) %>% 
    signif(2)
  
  
  
  # 4) get ggpredict data
  preds <- ggpredict(
    engine,
    terms = c(predictor, "podtyp_nemoci_zjednoduseny")
  )
  
  # 5) build the base plot
  # ct <- emmeans(engine, "podtyp_nemoci_zjednoduseny") |> 
  #   contrast(adjust="none") |> data.frame()
  
  
  fml <- as.formula(paste0("~", var_indep2))
  
  ct <- emtrends(engine,
                 fml,
                 var =predictor) |>
    summary(infer = c(TRUE, TRUE)) |> 
    data.frame()
  
  preds2 <- ggpredict(
    engine,
    terms = c(predictor, "podtyp_nemoci_zjednoduseny", "poradie_vysetrenia")
  )
  
  raw_data <- attr(preds2, "rawdata")
  
  
  # 1) Build a lookup table of p‐values per facet
  pval_df <- ct  |> 
    mutate(
      group   = as.character(.data[[var_indep2]]),         
      p.label = paste0("p = ", signif(p.value, 3))         
    )  |> 
    select(group, p.label)
  
  
  data_new <- data
  dep_var <- deparse(extract_fit_engine(model)@call$formula)[1] |> 
    str_extract("^[^ ]+")
  
  # 2) Your base plot
  p <- plot(preds) +
    geom_point(
      data = data_new |> 
        filter(!is.na(podtyp_nemoci_zjednoduseny)) |> 
        mutate(group = podtyp_nemoci_zjednoduseny),
      aes(x = !!sym(predictor), 
          y = !!sym(dep_var), 
          colour = odpoved_na_terapii_m0_vs_m6, 
          fill = odpoved_na_terapii_m0_vs_m6, 
          shape = odpoved_na_terapii_m0_vs_m6),
      size = 3,
      alpha       = 0.3,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(values = colors_20) +
    scale_fill_manual(values = colors_20) +
    scale_shape_manual(values = c(21, 22, 23, 24, 25,21, 22, 23, 24, 25), name = "Examination") +
    # paletteer::scale_colour_paletteer_d("tvthemes::Steven") +
    # paletteer::scale_fill_paletteer_d("tvthemes::Steven") +
    facet_wrap(~ group, scales = "free_y") +
    labs(
      title  = paste0("p-value = ", pval),
      x      = xlab %||% predictor,
      y      = ylab %||% "Predicted response",
      colour = ""
    ) +
    theme_minimal(base_size = 18) +
    theme(
      strip.text       = element_text(face = "bold", size = 22),
      axis.title       = element_text(face = "bold", size = 22),
      axis.text        = element_text(size = 21),
      plot.title       = element_text(face = "bold", size = 25, hjust = 0.5),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
  
  # 3) Add a geom_text layer to annotate each facet
  p <- p + 
    geom_text(
      data = pval_df,
      aes(x = Inf, y = Inf, label = p.label),
      hjust = 1.1,    # nudge left from right edge
      vjust = 1.1,    # nudge down from top
      size  = 6, 
      inherit.aes = FALSE
    )
  
  return(p)
  
}

plot_mixed_terms_05 <- function(model, 
                                var_pattern = "var_indep_",
                                xlab = NULL, 
                                ylab = NULL,
                                p_label_size = 6,
                                ...) {
  library(broom.mixed)
  library(sjPlot)
  library(ggplot2)
  library(sjmisc)    # for theme_sjplot2()
  
  # 1) grab fixed‐effect terms + p‐values
  fe <- broom.mixed::tidy(as_lmerModLmerTest(model$fit), effects = "fixed")
  
  # 2) select those matching your pattern
  sel_terms <- grep(var_pattern, fe$term, value = TRUE)
  if (length(sel_terms) == 0) {
    stop("No fixed-effect terms matched pattern: ", var_pattern)
  }
  if (length(sel_terms) > 1) {
    warning("More than one term matched; plotting all. p-value will use the first.")
  }
  
  # 3) extract the p-value of the first matched term
  pval <- fe %>% 
    filter(term == sel_terms[1])  |>  
    pull(p.value)
  
  
  
  
  fml <- Reduce(paste, 
                deparse(parsnip::extract_fit_engine(model) |> formula())) |> 
    str_replace("\\+\\s*poradie_vysetrenia",
                "* poradie_vysetrenia") |> 
    as.formula()
  
  
  model_int <- parsnip::fit(model$spec, fml, data = model$fit@frame)
  
  sjPlot::plot_model(
    model_int,
    type      = "int",
    show.data = TRUE
  )
  
  
  # 4) build the emmeans plot
  p <- sjPlot::plot_model(
    model_int,
    type      = "int",
    show.data = TRUE
  )
  
  # 5) apply axis labels if given
  if (!is.null(xlab) || !is.null(ylab)) {
    p <- p + labs(
      x = xlab %||% waiver(),
      y = ylab %||% waiver()
    )
  }
  
  # 6) annotate p-value in top-right
  p <- p +
    annotate(
      "text",
      x     = Inf,
      y     = Inf,
      label = paste0("p = ", signif(pval, 2)),
      hjust = 1.1,
      vjust = 1.5,
      size  = p_label_size
    ) +
    
    # 7) custom theme
    theme_sjplot2() +
    theme(
      plot.title       = element_blank(),
      plot.subtitle    = element_text(size = 20, color = "grey40", hjust = 0),
      axis.title       = element_text(face = "bold", size = 20),
      axis.text        = element_text(size = 18, color = "grey20"),
      panel.grid.major = element_line(linetype = "dotted", color = "grey80"),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}


# eda - podtyp_zjednodus across periods
plot_period_type_01 <- function(data,
                                ylab = NULL) {
  pjd <- position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)
  pd01 <- position_dodge(width = 0.6)
  
  fig <- data |>
    ggplot(aes(
      x = poradie_vysetrenia,
      y = var_dep_value,
      group = podtyp_nemoci_zjednoduseny,
      colour = podtyp_nemoci_zjednoduseny,
      na.rm = TRUE
    )) +
    theme_sjplot2() +
    geom_jitter(alpha = 0.2, size = 2, position = pjd, show.legend = FALSE) +
    stat_summary(
      fun.data = median_cl_boot, geom = "errorbar", B = 1000,   # B refers to boostraps iters 
      position = pd01, na.rm = TRUE, width = 0.3, size = 1
    ) +
    stat_summary(
      fun = median, geom = "crossbar", 
      position = pd01, na.rm = TRUE, width = 0.6, size = 0.3
    ) +
    stat_summary(
      fun = median, geom = "point",
      position = pd01, size = 2, show.legend = FALSE
    )
  
  
  if (!is.null(xlab) || !is.null(ylab)) {
    fig <- fig + labs(
      x = "Poradie vyšetrenia",
      y = ylab %||% waiver()
    )
  }
  
  fig <- fig +
    theme(axis.title = element_text(size = 22, face = "bold"),
          axis.text = element_text(size = 20)) +
    scale_colour_discrete(drop = TRUE, na.translate = FALSE) +
    scale_fill_discrete(drop = TRUE, na.translate = FALSE)
  
  return(fig)
}


plot_lmer_grid_basic_01 <- function(data,
                                    xlab = NULL, 
                                    ylab = NULL,
                                    ...) {
  
  library(dplyr)
  library(lmerTest)  # lmer with p-values
  library(emmeans)   # emtrends for simple slopes
  library(ggplot2)
  library(sjmisc)    # theme_sjplot2()
  
  # 1) Filter once (same as your plotting data)
  df_plot <- data %>%
    filter(!is.na(podtyp_nemoci_zjednoduseny),
           !is.na(odpoved_na_terapii_m0_vs_m6))
  
  # 2) Fit a single model with interactions so slopes can vary by panel
  #    (include interactions of var_indep_value with both factors;
  #     use the 3-way to allow fully cell-specific slopes)
  m <- lmer(
    var_dep_value ~ var_indep_value * podtyp_nemoci_zjednoduseny * poradie_vysetrenia +
      denna_davka_gk_mg_kg_bw + (1 | projekt_id),
    data = df_plot
  )
  
  # 3) Get per-panel trend (slope of var_indep_value) with p-values
  slopes <- emtrends(
    m,
    ~ podtyp_nemoci_zjednoduseny * poradie_vysetrenia,
    var = "var_indep_value"
  ) %>%
    summary(infer = c(TRUE, TRUE)) %>%   # gives SE, df, t-ratio, p.value, CI
    mutate(p_label = paste0("p = ", signif(p.value, 3)))
  
  # 4) Plot and annotate p-values in each facet
  fig <- ggplot(
    df_plot,
    aes(x = var_indep_value,
        y = var_dep_value,
        col  = podtyp_nemoci_zjednoduseny,
        fill = podtyp_nemoci_zjednoduseny)
  ) +
    scale_shape_manual(values = c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25),
                       name = "Treatment Effect") +
    theme_bw() +
    geom_point(alpha = 0.2, size = 3, aes(shape = odpoved_na_terapii_m0_vs_m6)) +
    stat_smooth(method = "lm", se = FALSE) +
    facet_grid(podtyp_nemoci_zjednoduseny ~ poradie_vysetrenia) +
    
    theme(axis.title = element_text(face = "bold", size = 20),
          axis.text = element_text(size = 18, color = "grey20"),
          strip.text = element_text(face = "bold", size = 20)
    ) +
    
    # annotate: match facet vars, pin to top-left
    geom_text(
      data = slopes,
      aes(label = p_label),
      x = Inf, y = -Inf,
      hjust = 1.1,   # push left from right edge
      vjust = -0.5,  # push up from bottom
      inherit.aes = FALSE,
      size = 8
    )
  
  if (!is.null(xlab) || !is.null(ylab)) {
    fig <- fig + labs(
      x = xlab %||% waiver(),
      y = ylab %||% waiver()
    )
  }
  
  return(fig)
  
}

plot_lmer_grid_basic_02 <- function(data,
                                    xlab = NULL, 
                                    ylab = NULL,
                                    covar_spec,
                                    covar = NULL,
                                    ...) {
  
  library(dplyr)
  library(lmerTest)  # lmer with p-values
  library(emmeans)   # emtrends for simple slopes
  library(ggplot2)
  library(sjmisc)    # theme_sjplot2()
  
  # 1) Filter once (same as your plotting data)
  df_plot <-data 
  
  
  ## formula specification
  if (is.null(covar)){
    rhs_terms <- paste(c("poradie_vysetrenia", "var_indep_value", covar_spec), 
                       collapse = " * ")
  } else{
    rhs_terms <- paste(paste(c("poradie_vysetrenia", "var_indep_value", covar_spec), 
                             collapse = " * "), " + ", covar)
  }
  
  fml <- as.formula(paste("var_dep_value ~", paste(rhs_terms, collapse = " + "), " + (1 | projekt_id)"))
  
  # 2) Fit a single model with interactions so slopes can vary by panel
  #    (include interactions of var_indep_value with both factors;
  #     use the 3-way to allow fully cell-specific slopes)
  
  m <- lmer(
    fml,
    data = df_plot
  )
  
  # 3) Get per-panel trend (slope of var_indep_value) with p-values
  spec <- as.formula(paste("~ ", covar_spec, "* poradie_vysetrenia"))
  spec2 <- as.formula(paste("~ ", covar_spec))
  
  
  slopes <- emtrends(
    m,
    specs = spec,
    var = "var_indep_value"
  ) |> 
    summary(infer = c(TRUE, TRUE)) |>    # gives SE, df, t-ratio, p.value, CI
    mutate(p_label = paste0("p = ", signif(p.value, 3)))
  
  slopes2 <- emtrends(
    m,
    specs = spec2,
    var = "var_indep_value"
  ) |> 
    summary(infer = c(TRUE, TRUE)) |>    # gives SE, df, t-ratio, p.value, CI
    mutate(p_label = paste0("p = ", signif(p.value, 3)),
           f_label = paste0(!! sym(covar_spec), " (",p_label, ")"))
  
  labs_vec  <- setNames(slopes2$f_label, as.character(slopes2[[covar_spec]]))
  labs_list <- list(); labs_list[[covar_spec]] <- labs_vec
  
  # 4) Plot and annotate p-values in each facet
  fig <- ggplot(
    df_plot |> 
      filter(!is.na(!! sym(covar_spec))),
    aes(x = var_indep_value,
        y = var_dep_value,
        colour = !! sym(covar_spec),
        fill   = !! sym(covar_spec))
  ) +
    scale_shape_manual(values = c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25),
                       name = "Treatment Effect") +
    theme_bw() +
    geom_point(alpha = 0.2, size = 3) +
    stat_smooth(method = "lm", se = FALSE) +
    facet_grid(
      rows = vars(!! sym(covar_spec)),
      cols = vars(poradie_vysetrenia),
      labeller = do.call(labeller, labs_list),
      drop = TRUE) +
    
    theme(axis.title = element_text(face = "bold", size = 20),
          axis.text = element_text(size = 18, color = "grey20"),
          strip.text = element_text(face = "bold", size = 20)
    ) +
    
    # annotate: match facet vars, pin to top-left
    geom_text(
      data = slopes,
      aes(label = p_label),
      x = Inf, y = -Inf,
      hjust = 1.1,   # push left from right edge
      vjust = -0.5,  # push up from bottom
      inherit.aes = FALSE,
      size = 8
    )
  
  if (!is.null(xlab) || !is.null(ylab)) {
    fig <- fig + labs(
      x = xlab %||% waiver(),
      y = ylab %||% waiver()
    )
  }
  
  return(fig)
  
}

# ggally
lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.15) +
    geom_smooth(method = method, se = FALSE, ...)
  p
}

# plot umap figures
plot_umap <- function(df, var) {
  ggplot(df, aes(x = UMAP1, y = UMAP2, color = .data[[var]])) +
    geom_point(size = 2, alpha = 0.9) +
    theme_minimal(base_size = 14) +
    labs(color = var, title = paste("UMAP –", var))
}

run_umap_plot <- function(model_famd, data, variable_to_plot) {
  
  # 1) Extract FAMD scores
  famd_scores <- model_famd$ind$coord
  
  # 2) Run UMAP
  umap_res <- umap::umap(famd_scores, n_neighbors = 15, min_dist = 0.1)
  
  # 3) Prepare dataframe
  umap_df <- as.data.frame(umap_res$layout) |>
    dplyr::rename(UMAP1 = V1, UMAP2 = V2) |>
    dplyr::mutate(
      !!variable_to_plot := data[[variable_to_plot]]
    )
  
  # 4) Create plot
  p <- plot_umap(umap_df, variable_to_plot)
  
  # 5) Save figure
  tiff(
    filename = paste0("output/figures/cluster analyses/",
                      format(Sys.Date(), "%y%m%d"),
                      "_umap_",
                      variable_to_plot,
                      ".tiff"),
    compression = "lzw",
    res = 300,
    width = 1600,
    height = 1200
  )
  print(p)
  dev.off()
  
  message("Saved: ", variable_to_plot)
}

## correlatoin kendall
make_corr_heatmap_kendall <- function(
    d_join,
    x_cols,
    y_cols,
    title_suffix = "ALL",
    legend_df = NULL
) {
  X <- d_join |>
    dplyr::select(dplyr::all_of(x_cols)) |>
    as.matrix() |>
    drop_const_cols()
  
  Y <- d_join |>
    dplyr::select(dplyr::all_of(y_cols)) |>
    as.matrix() |>
    drop_const_cols()
  
  cor_test_df <- expand.grid(
    activity_var = colnames(X),
    clinical_var = colnames(Y),
    stringsAsFactors = FALSE
  ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      test = list(
        stats::cor.test(
          x = X[, activity_var],
          y = Y[, clinical_var],
          method = "kendall",
          use = "pairwise.complete.obs",
          exact = FALSE
        )
      ),
      correlation = unname(test$estimate),
      p_value     = test$p.value
    ) |>
    dplyr::ungroup() |>
    dplyr::select(activity_var, clinical_var, correlation, p_value) |>
    dplyr::mutate(
      signif = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ ""
      )
    )
  
  # --- ORDER clinical variables by d11_legend$order (1..5), NAs last
  if (!is.null(legend_df)) {
    leg <- legend_df |>
      dplyr::transmute(
        clinical_var = names_new,
        ord = order
      ) |>
      dplyr::distinct()
    
    # keep only vars present in the plot
    leg_present <- leg |>
      dplyr::filter(clinical_var %in% unique(cor_test_df$clinical_var))
    
    # define order: ord 1..5 first, then ord NA, and within each group by name
    lev <- leg_present |>
      dplyr::mutate(ord2 = dplyr::if_else(is.na(ord), Inf, ord)) |>
      dplyr::arrange(ord2, clinical_var) |>
      dplyr::pull(clinical_var) |>
      unique()
    
    # append any clinical vars missing from legend at the end
    lev_missing <- setdiff(unique(cor_test_df$clinical_var), lev)
    lev <- c(lev, lev_missing)
    
    cor_test_df <- cor_test_df |>
      dplyr::mutate(clinical_var = factor(clinical_var, levels = lev))
  }
  
  fig <- ggplot2::ggplot(cor_test_df, ggplot2::aes(activity_var, clinical_var, fill = correlation)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = signif), color = "black", size = 3) +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title  = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      fill = "Kendall τ",
      title = paste0("Correlation: Activity vs Clinical (", title_suffix, ")"),
      subtitle = "* p < 0.05, ** p < 0.01, *** p < 0.001 (raw p-values)"
    )
  
  list(X = X, Y = Y, cor_df = cor_test_df, fig = fig)
}


# bout stacked chart
# 1) doplnit rest sloupce (na úrovni pacienta)
add_rest_cols <- function(df, ranges, bouts) {
  for (r in ranges) {
    total_col <- paste0(r, "_time")
    bout_cols <- paste0(r, "_", bouts, "_time")
    rest_col  <- paste0(r, "_rest_time")
    
    df <- df |>
      mutate(
        "{rest_col}" := pmax(
          .data[[total_col]] - rowSums(across(all_of(bout_cols)), na.rm = TRUE),
          0
        )
      )
  }
  df
}

# bezpečné uložení TIFF (ať nemáš 30× stejný blok)
save_tiff <- function(p, fname, res = 150, w = 1600, h = 1200) {
  tiff(filename = fname, compression = "lzw", res = res, width = w, height = h)
  print(p)
  dev.off()
}



# radar figure
plot_percentile_radar <- function(data, disease_type = "all") {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  
  # FILTR DAT

  if (!identical(disease_type, "all")) {
    
    data <- data |>
      dplyr::filter(subtype == as.numeric(disease_type))
    
    title_text <- paste0("subtype: ", disease_type)
    
  } else {
    
    title_text <- "subtype: all"
    
  }

  # RADIAL LABELS
  lab_data <- expand.grid(
    exam = unique(data$exam),
    y = c(50, 100, 200, 300, 400)
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      x = 0,
      label = as.character(y)
    )
  
  
  # MAPPING NAMES
  name_map <- c(
    "p9931_enmo" = "M010",
    "p9861_enmo" = "M020",
    "p9792_enmo" = "M030",
    "p9688_enmo" = "M045",
    "p9583_enmo" = "M060",
    "p9375_enmo" = "M090",
    "p9167_enmo" = "M120"
  )
  
  # MAIN PIPELINE
  data |> 
    dplyr::select(dplyr::starts_with("p9"), exam) |>
    tidyr::pivot_longer(
      cols = -exam,
      names_to = "var_name",
      values_to = "var_val"
    ) |>
    mutate(
      var_name = dplyr::recode(var_name, !!!name_map)
    ) |>
    dplyr::filter(!is.na(var_val)) |>
    dplyr::group_by(exam, var_name) |>
    dplyr::summarise(
      p05 = quantile(var_val, 0.05, na.rm = TRUE),
      p25 = quantile(var_val, 0.25, na.rm = TRUE),
      p50 = quantile(var_val, 0.50, na.rm = TRUE),
      p75 = quantile(var_val, 0.75, na.rm = TRUE),
      p95 = quantile(var_val, 0.95, na.rm = TRUE),
      p100 = max(var_val, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      `0-5`    = p05,
      `5-25`   = p25 - p05,
      `25-50`  = p50 - p25,
      `50-75`  = p75 - p50,
      `75-95`  = p95 - p75,
      `95-100` = p100 - p95
    ) |>
    dplyr::select(exam, var_name, `0-5`, `5-25`, `25-50`, `50-75`, `75-95`, `95-100`) |>
    tidyr::pivot_longer(
      cols = -c(exam, var_name),
      names_to = "percentil",
      values_to = "rozpeti"
    ) |>
    dplyr::mutate(
      percentil = factor(
        percentil,
        levels = c("95-100", "75-95", "50-75", "25-50", "5-25", "0-5")
      )
    ) |>
    ggplot(aes(x = rev(var_name), y = rozpeti, fill = percentil)) +
    geom_col(width = 1, linewidth = 0.1) +
    scale_fill_brewer(palette = "Reds", direction = -1) +
    coord_polar() +
    geom_abline(slope = 0, intercept = 100, col = "black", lty = 2) +
    geom_abline(slope = 0, intercept = 40, col = "black", lty = 2) +
    geom_text(
      data = lab_data,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE
    ) +
    facet_wrap(~exam) +
    labs(title = title_text) +
    theme_minimal()+
    theme(
      axis.title = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
}

plot_percentile_radar_02 <- function(data, tis_response_type = "all") {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  
  # příprava dat
  data <- data |>
    mutate(
      tis_response = as.factor(tis_response)
    ) |>
    filter(exam != "M0")
  
  # FILTR DAT
  if (!identical(tis_response_type, "all")) {
    
    data <- data |>
      filter(.data$tis_response == tis_response_type)
    
    title_text <- paste0("tis_response: ", tis_response_type)
    
  } else {
    
    title_text <- "tis_response: all"
    
  }
  
  # RADIAL LABELS
  lab_data <- expand.grid(
    exam = unique(data$exam),
    y = c(50, 100, 200, 300, 400)
  ) |>
    tibble::as_tibble() |>
    mutate(
      x = 0,
      label = as.character(y)
    )
  
  # MAPPING NAMES
  name_map <- c(
    "p9931_enmo" = "M010",
    "p9861_enmo" = "M020",
    "p9792_enmo" = "M030",
    "p9688_enmo" = "M045",
    "p9583_enmo" = "M060",
    "p9375_enmo" = "M090",
    "p9167_enmo" = "M120"
  )
  
  # MAIN PIPELINE
  data |> 
    select(starts_with("p9"), exam) |>
    pivot_longer(
      cols = -exam,
      names_to = "var_name",
      values_to = "var_val"
    ) |>
    mutate(
      var_name = dplyr::recode(var_name, !!!name_map)
    ) |>
    filter(!is.na(var_val)) |>
    group_by(exam, var_name) |>
    summarise(
      p05 = quantile(var_val, 0.05, na.rm = TRUE),
      p25 = quantile(var_val, 0.25, na.rm = TRUE),
      p50 = quantile(var_val, 0.50, na.rm = TRUE),
      p75 = quantile(var_val, 0.75, na.rm = TRUE),
      p95 = quantile(var_val, 0.95, na.rm = TRUE),
      p100 = max(var_val, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      `0-5`    = p05,
      `5-25`   = p25 - p05,
      `25-50`  = p50 - p25,
      `50-75`  = p75 - p50,
      `75-95`  = p95 - p75,
      `95-100` = p100 - p95
    ) |>
    select(exam, var_name, `0-5`, `5-25`, `25-50`, `50-75`, `75-95`, `95-100`) |>
    pivot_longer(
      cols = -c(exam, var_name),
      names_to = "percentil",
      values_to = "rozpeti"
    ) |>
    mutate(
      percentil = factor(
        percentil,
        levels = c("95-100", "75-95", "50-75", "25-50", "5-25", "0-5")
      )
    ) |>
    ggplot(aes(x = rev(var_name), y = rozpeti, fill = percentil)) +
    geom_col(width = 1, linewidth = 0.1) +
    scale_fill_brewer(palette = "Reds", direction = -1) +
    coord_polar() +
    geom_abline(slope = 0, intercept = 100, col = "black", lty = 2) +
    geom_abline(slope = 0, intercept = 40, col = "black", lty = 2) +
    geom_text(
      data = lab_data,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE
    ) +
    facet_wrap(~exam) +
    labs(title = title_text) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
}

plot_percentile_radar_tis_m6 <- function(data, tis_response_type = "all") {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  
  # vytvoření tis_response_2 (hodnota z M6 pro všechny exam pacienta)
  data <- data |>
    mutate(tis_response = as.factor(tis_response)) |>
    group_by(immet_id) |>
    mutate(
      tis_response_2 = tis_response[exam == "M6"][1]
    ) |>
    ungroup()
  
  # FILTR DAT
  if (!identical(tis_response_type, "all")) {
    
    data <- data |>
      filter(.data$tis_response_2 == tis_response_type)
    
    title_text <- paste0("tis_response (from M6): ", tis_response_type)
    
  } else {
    
    title_text <- "tis_response (from M6): all"
    
  }
  
  # RADIAL LABELS
  lab_data <- expand.grid(
    exam = unique(data$exam),
    y = c(50, 100, 200, 300, 400)
  ) |>
    tibble::as_tibble() |>
    mutate(
      x = 0,
      label = as.character(y)
    )
  
  # MAPPING NAMES
  name_map <- c(
    "p9931_enmo" = "M010",
    "p9861_enmo" = "M020",
    "p9792_enmo" = "M030",
    "p9688_enmo" = "M045",
    "p9583_enmo" = "M060",
    "p9375_enmo" = "M090",
    "p9167_enmo" = "M120"
  )
  
  # MAIN PIPELINE
  data |> 
    select(starts_with("p9"), exam) |>
    pivot_longer(
      cols = -exam,
      names_to = "var_name",
      values_to = "var_val"
    ) |>
    mutate(
      var_name = dplyr::recode(var_name, !!!name_map)
    ) |>
    filter(!is.na(var_val)) |>
    group_by(exam, var_name) |>
    summarise(
      p05 = quantile(var_val, 0.05, na.rm = TRUE),
      p25 = quantile(var_val, 0.25, na.rm = TRUE),
      p50 = quantile(var_val, 0.50, na.rm = TRUE),
      p75 = quantile(var_val, 0.75, na.rm = TRUE),
      p95 = quantile(var_val, 0.95, na.rm = TRUE),
      p100 = max(var_val, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      `0-5`    = p05,
      `5-25`   = p25 - p05,
      `25-50`  = p50 - p25,
      `50-75`  = p75 - p50,
      `75-95`  = p95 - p75,
      `95-100` = p100 - p95
    ) |>
    select(exam, var_name, `0-5`, `5-25`, `25-50`, `50-75`, `75-95`, `95-100`) |>
    pivot_longer(
      cols = -c(exam, var_name),
      names_to = "percentil",
      values_to = "rozpeti"
    ) |>
    mutate(
      percentil = factor(
        percentil,
        levels = c("95-100", "75-95", "50-75", "25-50", "5-25", "0-5")
      )
    ) |>
    ggplot(aes(x = rec(var_name), y = rozpeti, fill = percentil)) +
    geom_col(width = 1, linewidth = 0.1) +
    scale_fill_brewer(palette = "Reds", direction = -1) +
    coord_polar() +
    geom_abline(slope = 0, intercept = 100, col = "black", lty = 2) +
    geom_abline(slope = 0, intercept = 40, col = "black", lty = 2) +
    geom_text(
      data = lab_data,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE
    ) +
    facet_wrap(~exam) +
    labs(title = title_text) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
}

plot_percentile_radar_tis_m18 <- function(data, tis_response_type = "all") {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  
  # vytvoření tis_response_2 (hodnota z M18 pro všechny exam pacienta)
  data <- data |>
    mutate(
      tis_response = as.factor(tis_response)
    ) |>
    group_by(immet_id) |>
    mutate(
      tis_response_2 = tis_response[match("M18", exam)]
    ) |>
    ungroup()
  
  # FILTR DAT
  if (!identical(tis_response_type, "all")) {
    
    data <- data |>
      filter(.data$tis_response_2 == tis_response_type)
    
    title_text <- paste0("tis_response (from M18): ", tis_response_type)
    
  } else {
    
    title_text <- "tis_response (from M18): all"
    
  }
  
  # RADIAL LABELS
  lab_data <- expand.grid(
    exam = unique(data$exam),
    y = c(50, 100, 200, 300, 400)
  ) |>
    tibble::as_tibble() |>
    mutate(
      x = 0,
      label = as.character(y)
    )
  
  # MAPPING NAMES
  name_map <- c(
    "p9931_enmo" = "M010",
    "p9861_enmo" = "M020",
    "p9792_enmo" = "M030",
    "p9688_enmo" = "M045",
    "p9583_enmo" = "M060",
    "p9375_enmo" = "M090",
    "p9167_enmo" = "M120"
  )
  
  # MAIN PIPELINE
  data |> 
    select(starts_with("p9"), exam) |>
    pivot_longer(
      cols = -exam,
      names_to = "var_name",
      values_to = "var_val"
    ) |>
    mutate(
      var_name = dplyr::recode(var_name, !!!name_map)
    ) |>
    filter(!is.na(var_val)) |>
    group_by(exam, var_name) |>
    summarise(
      p05 = quantile(var_val, 0.05, na.rm = TRUE),
      p25 = quantile(var_val, 0.25, na.rm = TRUE),
      p50 = quantile(var_val, 0.50, na.rm = TRUE),
      p75 = quantile(var_val, 0.75, na.rm = TRUE),
      p95 = quantile(var_val, 0.95, na.rm = TRUE),
      p100 = max(var_val, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      `0-5`    = p05,
      `5-25`   = p25 - p05,
      `25-50`  = p50 - p25,
      `50-75`  = p75 - p50,
      `75-95`  = p95 - p75,
      `95-100` = p100 - p95
    ) |>
    select(exam, var_name, `0-5`, `5-25`, `25-50`, `50-75`, `75-95`, `95-100`) |>
    pivot_longer(
      cols = -c(exam, var_name),
      names_to = "percentil",
      values_to = "rozpeti"
    ) |>
    mutate(
      percentil = factor(
        percentil,
        levels = c("95-100", "75-95", "50-75", "25-50", "5-25", "0-5")
      )
    ) |>
    ggplot(aes(x = rev(var_name), y = rozpeti, fill = percentil)) +
    geom_col(width = 1, linewidth = 0.1) +
    scale_fill_brewer(palette = "Reds", direction = -1) +
    coord_polar() +
    geom_abline(slope = 0, intercept = 100, col = "black", lty = 2) +
    geom_abline(slope = 0, intercept = 40, col = "black", lty = 2) +
    geom_text(
      data = lab_data,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE
    ) +
    facet_wrap(~exam) +
    labs(title = title_text) +
    theme_minimal()
}


# tables ----
## column with 0, 1, and NA ----
is_01_col <- function(x) {
  all(unique(x) %in% c(0, 1, NA))
}

## drop zero-inflated columns ----
drop_const_cols <- function(M) {
  s <- apply(M, 2, sd, na.rm = TRUE)
  keep <- is.finite(s) & s > 0
  M[, keep, drop = FALSE]
}

## drop zero-inflated features - mafo ----
drop_const_features <- function(M) {
  s <- apply(M, 1, sd, na.rm = TRUE)
  keep <- is.finite(s) & s > 0
  M[keep, , drop = FALSE]
}


## model comparison ----
model_comp <- function(data, covar) {
  
  data <- data |> 
    select(var_dep_value, poradie_vysetrenia, var_indep_value, all_of(covar), projekt_id) |> 
    na.omit()
  
  rhs_terms <- c("poradie_vysetrenia", "var_indep_value", covar, "(1 | projekt_id)")
  fml <- as.formula(paste("var_dep_value ~", paste(rhs_terms, collapse = " + ")))
  
  ## tidymodels
  lmer_mod <- 
    linear_reg() |> 
    set_engine("lmer", REML = FALSE) |> 
    set_mode("regression")
  
  mod1 <- fit(
    lmer_mod,
    formula = var_dep_value ~ poradie_vysetrenia + var_indep_value + (1 | projekt_id),
    data = data
  )
  
  mod2 <- fit(
    lmer_mod,
    formula = fml,
    data = data
  )
  
  p_val_comp <- stats::anova(mod1$fit, mod2$fit)$`Pr(>Chisq)`[[2]]
  r2_mod <- performance(mod2, estimator = "ML")$R2_conditional
  p_val_covar <- model_parameters(mod2$fit)$p[[6]]
  
  return(list(table = tibble(p_val_comp = p_val_comp,
                             r2_mod = r2_mod,
                             p_val_covar = p_val_covar),
              model = mod2,
              model_noCovar = mod1)
  )
  
}

model_comp_int <- function(data, covar) {
  
  data <- data |> 
    dplyr::select(var_dep_value, poradie_vysetrenia, var_indep_value, 
                  all_of(covar), projekt_id) |> 
    na.omit()
  
  rhs_terms <- c("poradie_vysetrenia", paste0("var_indep_value*", covar), "(1 | projekt_id)")
  fml <- as.formula(paste("var_dep_value ~", paste(rhs_terms, collapse = " + ")))
  
  ## tidymodels
  lmer_mod <- 
    linear_reg() |> 
    set_engine("lmer", REML = FALSE) |> 
    set_mode("regression")
  
  mod1 <- fit(
    lmer_mod,
    formula = var_dep_value ~ poradie_vysetrenia + var_indep_value + (1 | projekt_id),
    data = data
  )
  
  mod2 <- fit(
    lmer_mod,
    formula = fml,
    data = data
  )
  
  p_val_comp <- stats::anova(mod1$fit, mod2$fit)$`Pr(>Chisq)`[[2]]
  r2_mod <- performance(mod2, estimator = "ML")$R2_conditional
  p_val_covar <- model_parameters(mod2$fit)$p[[7]]
  
  return(list(table = tibble(p_val_comp = p_val_comp,
                             r2_mod = r2_mod,
                             p_val_covar = p_val_covar),
              model = mod2,
              model_noCovar = mod1)
  )
  
}

model_min <- function(data, filter_na = TRUE) {
  
  if (filter_na){
    data <- data |> 
      select(var_dep_value, var_indep_value, poradie_vysetrenia, projekt_id) |> 
      na.omit()
  }
  
  ## tidymodels
  lmer_mod <- 
    linear_reg() |> 
    set_engine("lmer", REML = FALSE) |> 
    set_mode("regression")
  
  mod1 <- fit(
    lmer_mod,
    formula = var_dep_value ~ poradie_vysetrenia + var_indep_value + (1 | projekt_id),
    data = data
  )
  
  omega <- sjstats::anova_stats(mod1$fit)$partial.omegasq[[2]]
  r2_tbl <- performance::r2_nakagawa(mod1$fit)
  r2_mod <- r2_tbl$R2_conditional
  p_val_covar <- model_parameters(mod1$fit)$p[[5]]
  
  return(list(table = tibble(omega_sq_p = omega,
                             r2_mod = r2_mod,
                             p_val_covar = p_val_covar),
              model = mod1)
  )
  
  
}

model_min_2 <- function(data, filter_na = TRUE) {
  
  if (filter_na){
    data <- data |> 
      select(var_dep_value, var_indep_value, poradie_vysetrenia, projekt_id) |> 
      na.omit()
  }
  
  ## tidymodels
  lmer_mod <- 
    linear_reg() |> 
    set_engine("lmer", REML = FALSE) |> 
    set_mode("regression")
  
  mod1 <- fit(
    lmer_mod,
    formula = var_dep_value ~ poradie_vysetrenia + var_indep_value + (1 | projekt_id),
    data = data
  )
  
  omega <- sjstats::anova_stats(mod1$fit)$partial.omegasq[[2]]
  r2_tbl <- performance::r2_nakagawa(mod1$fit)
  r2_mod <- r2_tbl$R2_conditional
  p_val_covar <- model_parameters(mod1$fit)$p[[3]]
  
  return(list(table = tibble(omega_sq_p = omega,
                             r2_mod = r2_mod,
                             p_val_covar = p_val_covar),
              model = mod1)
  )
  
  
}

model_min_3 <- function(data, filter_na = TRUE) {
  
  if (filter_na){
    data <- data |> 
      select(var_dep_value, var_indep_value, poradie_vysetrenia, projekt_id) |> 
      na.omit()
  }
  
  ## tidymodels
  lmer_mod <- 
    linear_reg() |> 
    set_engine("lmer", REML = FALSE) |> 
    set_mode("regression")
  
  mod1 <- fit(
    lmer_mod,
    formula = var_dep_value ~ poradie_vysetrenia + var_indep_value + (1 | projekt_id),
    data = data
  )
  
  omega <- sjstats::anova_stats(mod1$fit)$partial.omegasq[[2]]
  r2_tbl <- performance::r2_nakagawa(mod1$fit)
  r2_mod <- r2_tbl$R2_conditional
  p_val_covar <- model_parameters(mod1$fit)$p[[4]]
  
  return(list(table = tibble(omega_sq_p = omega,
                             r2_mod = r2_mod,
                             p_val_covar = p_val_covar),
              model = mod1)
  )
  
  
}


model_covar <- function(data, covar, filter_na = TRUE) {
  
  if (filter_na){
    data <- data |> 
      select(var_dep_value, var_indep_value, poradie_vysetrenia, projekt_id) |> 
      na.omit()
  }
  
  ## formula specification
  rhs_terms <- c("poradie_vysetrenia", "var_indep_value", covar, "(1 | projekt_id)")
  fml <- as.formula(paste("var_dep_value ~", paste(rhs_terms, collapse = " + ")))
  
  ## tidymodels
  lmer_mod <- 
    linear_reg() |> 
    set_engine("lmer", REML = FALSE) |> 
    set_mode("regression")
  
  mod1 <- fit(
    lmer_mod,
    formula = fml,
    data = data
  )
  
  
  omega <- sjstats::anova_stats(mod1$fit)$partial.omegasq[[2]]
  r2_tbl <- performance::r2_nakagawa(mod1$fit)
  r2_mod <- r2_tbl$R2_conditional
  p_val_covar <- model_parameters(mod1$fit)$p[[5]]
  
  return(list(table = tibble(omega_sq_p = omega,
                             r2_mod = r2_mod,
                             p_val_covar = p_val_covar),
              model = mod1)
  )
  
}

model_covar_2 <- function(data, covar, filter_na = TRUE) {
  
  if (filter_na){
    data <- data |> 
      select(var_dep_value, var_indep_value, poradie_vysetrenia, projekt_id) |> 
      na.omit()
  }
  
  ## formula specification
  rhs_terms <- c("poradie_vysetrenia", "var_indep_value", covar, "(1 | projekt_id)")
  fml <- as.formula(paste("var_dep_value ~", paste(rhs_terms, collapse = " + ")))
  
  ## tidymodels
  lmer_mod <- 
    linear_reg() |> 
    set_engine("lmer", REML = FALSE) |> 
    set_mode("regression")
  
  mod1 <- fit(
    lmer_mod,
    formula = fml,
    data = data
  )
  
  
  omega <- sjstats::anova_stats(mod1$fit)$partial.omegasq[[2]]
  r2_tbl <- performance::r2_nakagawa(mod1$fit)
  r2_mod <- r2_tbl$R2_conditional
  p_val_covar <- model_parameters(mod1$fit)$p[[3]]
  
  return(list(table = tibble(omega_sq_p = omega,
                             r2_mod = r2_mod,
                             p_val_covar = p_val_covar),
              model = mod1)
  )
  
}


model_covar_3 <- function(data, covar, filter_na = TRUE) {
  
  if (filter_na){
    data <- data |> 
      select(var_dep_value, var_indep_value, poradie_vysetrenia, projekt_id) |> 
      na.omit()
  }
  
  ## formula specification
  rhs_terms <- c("poradie_vysetrenia", "var_indep_value", covar, "(1 | projekt_id)")
  fml <- as.formula(paste("var_dep_value ~", paste(rhs_terms, collapse = " + ")))
  
  ## tidymodels
  lmer_mod <- 
    linear_reg() |> 
    set_engine("lmer", REML = FALSE) |> 
    set_mode("regression")
  
  mod1 <- fit(
    lmer_mod,
    formula = fml,
    data = data
  )
  
  
  omega <- sjstats::anova_stats(mod1$fit)$partial.omegasq[[2]]
  r2_tbl <- performance::r2_nakagawa(mod1$fit)
  r2_mod <- r2_tbl$R2_conditional
  p_val_covar <- model_parameters(mod1$fit)$p[[4]]
  
  return(list(table = tibble(omega_sq_p = omega,
                             r2_mod = r2_mod,
                             p_val_covar = p_val_covar),
              model = mod1)
  )
  
}

futuremice_safe <- function(data, method = "pmm", m = 20, maxit = 15, workers = 4, parallel = FALSE) {
  # Zruš zbytky starých paralelních procesů
  future::plan(future::sequential)
  gc()
  
  # Nastav paralelní plán (volitelně)
  if (isTRUE(parallel)) {
    future::plan(future::multisession, workers = workers)
  }
  
  # Spusť imputaci
  imp <- tryCatch(
    mice::futuremice(
      data,
      method = method,
      m = m,
      maxit = maxit,
      ridge = 1e-02,
      remove.collinear = TRUE,
      remove.constant = TRUE
    ),
    error = function(e) {
      message("futuremice_safe: retry with higher ridge (1e-01). Original error: ", e$message)
      mice::futuremice(
        data,
        method = method,
        m = m,
        maxit = maxit,
        ridge = 1e-01,
        remove.collinear = TRUE,
        remove.constant = TRUE
      )
    }
  )
  
  # Po skončení vypni workery
  future::plan(future::sequential)
  
  # Vrací kompletní data (defaultně první set, můžeš změnit)
  mice::complete(imp)
}

# **********************************************************************
### rCCA -> export tabulek do XLSX (multi-sheet) ----
# **********************************************************************

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
  
  # ---- efektivní počet komponent (když by něco nesedělo)
  ncomp_eff <- min(ncomp, ncol(loadX), ncol(loadY), ncol(fit$variates$X), ncol(fit$variates$Y))
  
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
  
  # ---- loadings long (všechny komponenty, indexově)
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
  
  # ---- výběr top proměnných pro cross-correlation (COMP1 = sloupec 1)
  # top dle loadings comp1 (ne podle názvu)
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
  
  # cor() někdy vrátí vektor -> vynucení matice
  if (is.null(dim(cor_raw))) {
    cor_mat <- matrix(cor_raw, nrow = length(top_x_vars), ncol = length(top_y_vars),
                      dimnames = list(top_x_vars, top_y_vars))
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
    c("subtype", "response", "sex", "disease_subtype", "response_m0_m6"),
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

## export ----
format_corr_table <- function(cor_df) {
  cor_df |>
    dplyr::mutate(
      `Activity variable`  = activity_var,
      `Clinical variable`  = clinical_var,
      `Kendall tau`        = round(correlation, 3),
      `P-value`            = signif(p_value, 3),
      `Significance`       = signif
    ) |>
    dplyr::select(
      `Activity variable`,
      `Clinical variable`,
      `Kendall tau`,
      `P-value`,
      `Significance`
    ) |>
    dplyr::arrange(desc(abs(`Kendall tau`)))
}


## Possibly ----
possLMER <- possibly(.f=lmer, otherwise = NA_real_)
possmixMODmmt <- possibly(mod_mixMODmmt, otherwise = NULL)
# possMODELCOMP <- possibly(.f=model_comp, otherwise = NA_real_)
possCLMM <- possibly(.f = ordinal::clmm, otherwise = NA_real_)

## tidymodels ----
lmer_mod <- 
  linear_reg() |> 
  set_engine("lmer") |> 
  set_mode("regression")

lmer_tidier <- function(mod){
  mod_tidier <- 
    as_lmerModLmerTest(mod$fit) |> broom.mixed::tidy(conf.int = TRUE)
  return(mod_tidier)
} 



