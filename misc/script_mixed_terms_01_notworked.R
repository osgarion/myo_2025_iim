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
  # pval_df <- ct %>% 
  #   # strip off the common prefix so that it matches your 'group' levels
  #   mutate(group = str_extract(contrast, "\\d+"),
  #          p.label = paste0("p = ", signif(p.value, 3))) %>%
  #   select(group, p.label)
  
  pval_df <- ct %>%
    mutate(
      group   = as.character(.data[[var_indep2]]),         # vezmeme hodnoty ve sloupci "pohlavi"
      p.label = paste0("p = ", signif(p.value, 3))         # naformátujeme p‑hodnotu
    ) %>%
    select(group, p.label)
  
  
  data_new <- data
  dep_var <- deparse(extract_fit_engine(model)@call$formula)[1] |> 
    str_extract("^[^ ]+")
  
  
  # fml_new <- extract_fit_engine(model)@call$formula |> 
  #   deparse() |> 
  #   paste(collapse = " ") |> 
  #   str_squish() |> 
  #   str_replace("(\\+\\s*)podtyp_nemoci_zjednoduseny",
  #               "* podtyp_nemoci_zjednoduseny") |> 
  #   str_replace("(\\+\\s*)var_indep_value",
  #               "* var_indep_value") |>
  #   as.formula()
  # 
  data_new <- data
  # 
  # preds3 <- fit(lmer_mod, formula = fml_new, data = data_new)
  # 
  # predictor2 <- paste0(predictor, " [all]")
  # 
  # preds_int <- ggpredict(
  #   preds3,
  #   terms = c(predictor2, "podtyp_nemoci_zjednoduseny")
  # )   
  
  dep_var <- deparse(fml_new)[1] |> str_extract("^[^ ]+")
  
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