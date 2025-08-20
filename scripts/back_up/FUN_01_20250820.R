# Libraries & functions ----
## libraries ----
pacman::p_load(update = T,  
               equatiomatic, beepr, tictoc, 
               tidyverse, purrr, furrr, easystats, rio, janitor, ggthemes, car,
               gtsummary, skimr, sjPlot, flextable, ggpubr, rstatix, tidymodels,
               kableExtra, skimr, GGally, testthat, factoextra, gplots, uwot,
               lmerTest, dlookr, multilevelmod,furrr,ggforce, lazyWeave, paletteer,
               emmeans, openxlsx, GLMMadaptive
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
      data = data_new,
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
                                xlab = NULL,
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
      x = xlab %||% waiver(),
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

# ggally
lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.15) +
    geom_smooth(method = method, se = FALSE, ...)
  p
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

## tidymodels
lmer_mod <- 
  linear_reg() |> 
  set_engine("lmer") |> 
  set_mode("regression")

lmer_tidier <- function(mod){
  mod_tidier <- 
    as_lmerModLmerTest(mod$fit) |> broom.mixed::tidy(conf.int = TRUE)
  return(mod_tidier)
} 
