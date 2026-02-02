# **********************************************************************
# 4. rCCA ----   (MYOKINY jako Y blok místo clinics; jen M0 + M6)
# **********************************************************************


# 0) Helpers / config

iim_myo_rcca_myoX_prefix <- paste0(format(Sys.Date(), "%y%m%d"), "_iim_myo_rcca_myoX")
iim_myo_rcca_myoX_outdir <- file.path("output/figures/final")
dir.create(iim_myo_rcca_myoX_outdir, recursive = TRUE, showWarnings = FALSE)

norm_exam_m0m6 <- function(x) {
  x <- tolower(as.character(x))
  dplyr::case_when(
    x %in% c("m0", "0") ~ "m0",
    x %in% c("m6", "6") ~ "m6",
    TRUE ~ x
  )
}

# robustní sjednocení immet_id:
# - pokud je numeric -> vytvoří IMMETxx
# - pokud je character a vypadá jako "IMMET01" -> nechá
# - pokud je character a je jen číslo ("1") -> vytvoří IMMETxx
norm_immet_id <- function(x) {
  if (is.numeric(x)) {
    return(paste0("IMMET", stringr::str_pad(as.integer(x), 2, pad = "0")))
  }
  x_chr <- as.character(x)
  # pokud už začíná IMMET (libovolná velikost písmen), nech
  is_immet <- stringr::str_detect(toupper(x_chr), "^IMMET\\d+")
  out <- x_chr
  # pokud není IMMET a je to čistě číslo -> vytvoř IMMETxx
  is_num_string <- !is_immet & stringr::str_detect(x_chr, "^\\d+$")
  out[is_num_string] <- paste0("IMMET", stringr::str_pad(as.integer(out[is_num_string]), 2, pad = "0"))
  out
}


# 1) MYOKINY (Y block) preparation
#    - bez bbm_id
#    - jen M0/M6
#    - jen numerické myokiny

iim_myo_rcca_myoX_myok_df <- d01_myokiny_2 |>
  dplyr::mutate(
    immet_id   = norm_immet_id(immet_id),
    exam_order = factor(norm_exam_m0m6(exam_order), levels = c("m0", "m6"))
  ) |>
  dplyr::filter(exam_order %in% c("m0", "m6")) |>
  dplyr::select(-bbm_id)

iim_myo_rcca_myoX_myok_cols <- iim_myo_rcca_myoX_myok_df |>
  dplyr::select(-immet_id, -exam_order) |>
  dplyr::select(where(is.numeric)) |>
  names()

# ochrana: kdyby v myokinech nebyl žádný numeric sloupec
stopifnot(length(iim_myo_rcca_myoX_myok_cols) > 0)


# 2) AKTIVITA: X = numeric columns from d10_aktivita
#    Y = myokiny aligned by immet_id + exam_order

iim_myo_rcca_myoX_xcols_activ <- d10_aktivita |>
  dplyr::select(where(is.numeric)) |>
  names()

iim_myo_rcca_myoX_d_activ <- d12_join |>
  dplyr::mutate(
    immet_id   = norm_immet_id(immet_id),
    exam_order = factor(norm_exam_m0m6(exam_order), levels = c("m0", "m6"))
  ) |>
  dplyr::filter(exam_order %in% c("m0", "m6")) |>
  dplyr::inner_join(
    iim_myo_rcca_myoX_myok_df |>
      dplyr::select(immet_id, exam_order, dplyr::all_of(iim_myo_rcca_myoX_myok_cols)),
    by = c("immet_id", "exam_order")
  )

iim_myo_rcca_myoX_X_activ <- iim_myo_rcca_myoX_d_activ |>
  dplyr::select(dplyr::all_of(iim_myo_rcca_myoX_xcols_activ)) |>
  as.matrix()

iim_myo_rcca_myoX_Y_activ <- iim_myo_rcca_myoX_d_activ |>
  dplyr::select(dplyr::all_of(iim_myo_rcca_myoX_myok_cols)) |>
  as.matrix()

stopifnot(nrow(iim_myo_rcca_myoX_X_activ) == nrow(iim_myo_rcca_myoX_Y_activ))

# drop constant cols + scale
iim_myo_rcca_myoX_X_activ <- drop_const_cols(iim_myo_rcca_myoX_X_activ)
iim_myo_rcca_myoX_Y_activ <- drop_const_cols(iim_myo_rcca_myoX_Y_activ)

iim_myo_rcca_myoX_X_activ <- scale(iim_myo_rcca_myoX_X_activ)
iim_myo_rcca_myoX_Y_activ <- scale(iim_myo_rcca_myoX_Y_activ)

# rcc neumí NA -> vyřadit řádky s NA
iim_myo_rcca_myoX_keep_activ <- stats::complete.cases(iim_myo_rcca_myoX_X_activ) &
  stats::complete.cases(iim_myo_rcca_myoX_Y_activ)

iim_myo_rcca_myoX_X_activ <- iim_myo_rcca_myoX_X_activ[iim_myo_rcca_myoX_keep_activ, , drop = FALSE]
iim_myo_rcca_myoX_Y_activ <- iim_myo_rcca_myoX_Y_activ[iim_myo_rcca_myoX_keep_activ, , drop = FALSE]
iim_myo_rcca_myoX_d_activ <- iim_myo_rcca_myoX_d_activ[iim_myo_rcca_myoX_keep_activ, , drop = FALSE]

# po vyřazení řádků ještě jednou konstanty + scale (bezpečné)
iim_myo_rcca_myoX_X_activ <- drop_const_cols(iim_myo_rcca_myoX_X_activ) |> scale()
iim_myo_rcca_myoX_Y_activ <- drop_const_cols(iim_myo_rcca_myoX_Y_activ) |> scale()


# 3) ČAS: X = numeric columns from d10_cas
#    Y = myokiny aligned by immet_id + exam_order

iim_myo_rcca_myoX_xcols_time <- d10_cas |>
  dplyr::select(where(is.numeric)) |>
  names()

iim_myo_rcca_myoX_d_time <- d13_join |>
  dplyr::mutate(
    immet_id   = norm_immet_id(immet_id),
    exam_order = factor(norm_exam_m0m6(exam_order), levels = c("m0", "m6"))
  ) |>
  dplyr::filter(exam_order %in% c("m0", "m6")) |>
  dplyr::inner_join(
    iim_myo_rcca_myoX_myok_df |>
      dplyr::select(immet_id, exam_order, dplyr::all_of(iim_myo_rcca_myoX_myok_cols)),
    by = c("immet_id", "exam_order")
  )

iim_myo_rcca_myoX_X_time <- iim_myo_rcca_myoX_d_time |>
  dplyr::select(dplyr::all_of(iim_myo_rcca_myoX_xcols_time)) |>
  as.matrix()

iim_myo_rcca_myoX_Y_time <- iim_myo_rcca_myoX_d_time |>
  dplyr::select(dplyr::all_of(iim_myo_rcca_myoX_myok_cols)) |>
  as.matrix()

stopifnot(nrow(iim_myo_rcca_myoX_X_time) == nrow(iim_myo_rcca_myoX_Y_time))

# drop constant cols + scale
iim_myo_rcca_myoX_X_time <- drop_const_cols(iim_myo_rcca_myoX_X_time)
iim_myo_rcca_myoX_Y_time <- drop_const_cols(iim_myo_rcca_myoX_Y_time)

iim_myo_rcca_myoX_X_time <- scale(iim_myo_rcca_myoX_X_time)
iim_myo_rcca_myoX_Y_time <- scale(iim_myo_rcca_myoX_Y_time)

# rcc neumí NA -> vyřadit řádky s NA
iim_myo_rcca_myoX_keep_time <- stats::complete.cases(iim_myo_rcca_myoX_X_time) &
  stats::complete.cases(iim_myo_rcca_myoX_Y_time)

iim_myo_rcca_myoX_X_time <- iim_myo_rcca_myoX_X_time[iim_myo_rcca_myoX_keep_time, , drop = FALSE]
iim_myo_rcca_myoX_Y_time <- iim_myo_rcca_myoX_Y_time[iim_myo_rcca_myoX_keep_time, , drop = FALSE]
iim_myo_rcca_myoX_d_time <- iim_myo_rcca_myoX_d_time[iim_myo_rcca_myoX_keep_time, , drop = FALSE]

# po vyřazení řádků ještě jednou konstanty + scale
iim_myo_rcca_myoX_X_time <- drop_const_cols(iim_myo_rcca_myoX_X_time) |> scale()
iim_myo_rcca_myoX_Y_time <- drop_const_cols(iim_myo_rcca_myoX_Y_time) |> scale()

# **********************************************************************
# 4) Analyses (rCCA fits) ----
# **********************************************************************
iim_myo_rcca_myoX_fit_activ <- mixOmics::rcc(
  X = iim_myo_rcca_myoX_X_activ,
  Y = iim_myo_rcca_myoX_Y_activ,
  ncomp  = 2,
  method = "shrinkage"
)

iim_myo_rcca_myoX_fit_time <- mixOmics::rcc(
  X = iim_myo_rcca_myoX_X_time,
  Y = iim_myo_rcca_myoX_Y_time,
  ncomp  = 2,
  method = "shrinkage"
)

# optional: quick look at top loadings
sort(iim_myo_rcca_myoX_fit_activ$loadings$X[, 1], decreasing = TRUE)[1:10]
sort(iim_myo_rcca_myoX_fit_activ$loadings$Y[, 1], decreasing = TRUE)[1:10]

sort(iim_myo_rcca_myoX_fit_time$loadings$X[, 1], decreasing = TRUE)[1:10]
sort(iim_myo_rcca_myoX_fit_time$loadings$Y[, 1], decreasing = TRUE)[1:10]

# **********************************************************************
# 5) AKTIVITA outputs ---- (Y = MYOKINY)
# **********************************************************************

# --- Canonical scores scatter (comp1)
iim_myo_rcca_myoX_scores_activ <- tibble::tibble(
  immet_id   = iim_myo_rcca_myoX_d_activ$immet_id,
  exam_order = iim_myo_rcca_myoX_d_activ$exam_order,
  X1 = iim_myo_rcca_myoX_fit_activ$variates$X[, 1],
  Y1 = iim_myo_rcca_myoX_fit_activ$variates$Y[, 1],
  X2 = iim_myo_rcca_myoX_fit_activ$variates$X[, 2],
  Y2 = iim_myo_rcca_myoX_fit_activ$variates$Y[, 2]
)

iim_myo_rcca_myoX_fig_can_activ_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_scores_activ,
  ggplot2::aes(X1, Y1, color = exam_order)
) +
  ggplot2::geom_point(alpha = 0.8) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, linetype = 2) +
  ggplot2::labs(
    title = "rCCA: Canonical variates (Component 1)",
    x = "X variate 1 (activity block)",
    y = "Y variate 1 (myokine block)",
    color = "Exam"
  ) +
  sjPlot::theme_sjplot2()

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_canon_comp_activ_01.tiff")),
  compression = "lzw", res = 300, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_can_activ_01)
grDevices::dev.off()

# --- Heatmap: cross-correlation top loadings
iim_myo_rcca_myoX_topX_activ <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_X_activ),
  loading = iim_myo_rcca_myoX_fit_activ$loadings$X[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:15)

iim_myo_rcca_myoX_topY_activ <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_Y_activ),
  loading = iim_myo_rcca_myoX_fit_activ$loadings$Y[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:20)

iim_myo_rcca_myoX_cor_activ <- stats::cor(
  iim_myo_rcca_myoX_X_activ[, iim_myo_rcca_myoX_topX_activ$var, drop = FALSE],
  iim_myo_rcca_myoX_Y_activ[, iim_myo_rcca_myoX_topY_activ$var, drop = FALSE]
)

iim_myo_rcca_myoX_cor_df_activ <- as.data.frame(as.table(iim_myo_rcca_myoX_cor_activ)) |>
  stats::setNames(c("x_var", "y_var", "cor"))

iim_myo_rcca_myoX_fig_hm_activ_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_cor_df_activ,
  ggplot2::aes(x_var, y_var, fill = cor)
) +
  ggplot2::geom_tile() +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Cross-correlation heatmap (top activity vs top myokines)",
    x = "Activity (top loadings)",
    y = "Myokines (top loadings)",
    fill = "cor"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", space = "Lab")

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_heatmap_loading_activ_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_hm_activ_01)
grDevices::dev.off()

# --- Latent axis loadings (barplots)
iim_myo_rcca_myoX_topX_tab_activ <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_X_activ),
  loading = iim_myo_rcca_myoX_fit_activ$loadings$X[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:20)

iim_myo_rcca_myoX_topY_tab_activ <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_Y_activ),
  loading = iim_myo_rcca_myoX_fit_activ$loadings$Y[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:30)

iim_myo_rcca_myoX_load_df_activ <- dplyr::bind_rows(
  iim_myo_rcca_myoX_topX_tab_activ |> dplyr::mutate(type = "Activity"),
  iim_myo_rcca_myoX_topY_tab_activ |> dplyr::mutate(type = "Myokines")
)

iim_myo_rcca_myoX_fig_axis_activ_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_load_df_activ,
  ggplot2::aes(x = loading, y = reorder(var, loading), fill = type)
) +
  ggplot2::geom_col(alpha = 0.8) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::facet_wrap(~type, scales = "free_y") +
  ggplot2::labs(
    x = "Loading (Canonical component 1)",
    y = NULL,
    title = "Opposing poles of the latent activity–myokine axis",
    subtitle = "Variables on opposite sides represent contrasting patterns"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_activity_myokine_axis_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_axis_activ_01)
grDevices::dev.off()

# --- Correlation with canonical axis (tiles)
iim_myo_rcca_myoX_cor_axis_activ <- dplyr::bind_rows(
  tibble::tibble(
    var  = colnames(iim_myo_rcca_myoX_X_activ),
    cor  = stats::cor(iim_myo_rcca_myoX_X_activ, iim_myo_rcca_myoX_fit_activ$variates$X[, 1]),
    type = "Activity"
  ),
  tibble::tibble(
    var  = colnames(iim_myo_rcca_myoX_Y_activ),
    cor  = stats::cor(iim_myo_rcca_myoX_Y_activ, iim_myo_rcca_myoX_fit_activ$variates$Y[, 1]),
    type = "Myokines"
  )
) |>
  dplyr::filter(abs(cor) > 0.15)

iim_myo_rcca_myoX_fig_cor_activ_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_cor_axis_activ,
  ggplot2::aes(x = type, y = reorder(var, cor), fill = cor)
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  ggplot2::labs(
    fill = "Correlation",
    title = "Correlation of individual variables with the latent CCA axis",
    subtitle = "Red vs blue indicate opposing patterns"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_canon_cor_activ_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_cor_activ_01)
grDevices::dev.off()

# --- Subtype and response visualization
iim_myo_rcca_myoX_plot_activ <- iim_myo_rcca_myoX_d_activ |>
  dplyr::transmute(
    immet_id, exam_order,
    subtype  = as.factor(disease_subtype),
    response = as.factor(response_m0_m6),
    compX1   = iim_myo_rcca_myoX_fit_activ$variates$X[, 1],
    compY1   = iim_myo_rcca_myoX_fit_activ$variates$Y[, 1]
  )

iim_myo_rcca_myoX_fig_subtype_activ_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_plot_activ,
  ggplot2::aes(compX1, compY1, color = subtype)
) +
  ggplot2::geom_point(size = 2.8, alpha = 0.85) +
  ggplot2::stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
  ggplot2::labs(
    x = "rCCA component 1 (Activity)",
    y = "rCCA component 1 (Myokines)",
    title = "Patients projected onto the shared activity–myokine space",
    subtitle = "Ellipses highlight subtype clustering (80% t-ellipse)"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_rcca_subtype_activ_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_subtype_activ_01)
grDevices::dev.off()

iim_myo_rcca_myoX_fig_response_activ_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_plot_activ,
  ggplot2::aes(compX1, compY1, color = response)
) +
  ggplot2::geom_point(size = 2.8, alpha = 0.85) +
  ggplot2::stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
  ggplot2::labs(
    x = "rCCA component 1 (Activity)",
    y = "rCCA component 1 (Myokines)",
    title = "Patients projected onto the shared activity–myokine space",
    subtitle = "Ellipses highlight response clustering (80% t-ellipse)"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_rcca_response_activ_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_response_activ_01)
grDevices::dev.off()

# **********************************************************************
# 6) ČAS outputs ---- (Y = MYOKINY)
# **********************************************************************

# --- Canonical scores scatter (comp1)
iim_myo_rcca_myoX_scores_time <- tibble::tibble(
  immet_id   = iim_myo_rcca_myoX_d_time$immet_id,
  exam_order = iim_myo_rcca_myoX_d_time$exam_order,
  X1 = iim_myo_rcca_myoX_fit_time$variates$X[, 1],
  Y1 = iim_myo_rcca_myoX_fit_time$variates$Y[, 1],
  X2 = iim_myo_rcca_myoX_fit_time$variates$X[, 2],
  Y2 = iim_myo_rcca_myoX_fit_time$variates$Y[, 2]
)

iim_myo_rcca_myoX_fig_can_time_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_scores_time,
  ggplot2::aes(X1, Y1, color = exam_order)
) +
  ggplot2::geom_point(alpha = 0.8) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, linetype = 2) +
  ggplot2::labs(
    title = "rCCA: Canonical variates (Component 1)",
    x = "X variate 1 (time block)",
    y = "Y variate 1 (myokine block)",
    color = "Exam"
  ) +
  sjPlot::theme_sjplot2()

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_canon_comp_time_01.tiff")),
  compression = "lzw", res = 300, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_can_time_01)
grDevices::dev.off()

# --- Heatmap: cross-correlation top loadings
iim_myo_rcca_myoX_topX_time <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_X_time),
  loading = iim_myo_rcca_myoX_fit_time$loadings$X[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:15)

iim_myo_rcca_myoX_topY_time <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_Y_time),
  loading = iim_myo_rcca_myoX_fit_time$loadings$Y[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:20)

iim_myo_rcca_myoX_cor_time <- stats::cor(
  iim_myo_rcca_myoX_X_time[, iim_myo_rcca_myoX_topX_time$var, drop = FALSE],
  iim_myo_rcca_myoX_Y_time[, iim_myo_rcca_myoX_topY_time$var, drop = FALSE]
)

iim_myo_rcca_myoX_cor_df_time <- as.data.frame(as.table(iim_myo_rcca_myoX_cor_time)) |>
  stats::setNames(c("x_var", "y_var", "cor"))

iim_myo_rcca_myoX_fig_hm_time_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_cor_df_time,
  ggplot2::aes(x_var, y_var, fill = cor)
) +
  ggplot2::geom_tile() +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Cross-correlation heatmap (top time vs top myokines)",
    x = "Time (top loadings)",
    y = "Myokines (top loadings)",
    fill = "cor"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", space = "Lab")

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_heatmap_loading_time_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_hm_time_01)
grDevices::dev.off()

# --- Latent axis loadings (barplots)
iim_myo_rcca_myoX_topX_tab_time <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_X_time),
  loading = iim_myo_rcca_myoX_fit_time$loadings$X[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:20)

iim_myo_rcca_myoX_topY_tab_time <- tibble::tibble(
  var = colnames(iim_myo_rcca_myoX_Y_time),
  loading = iim_myo_rcca_myoX_fit_time$loadings$Y[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::slice(1:30)

iim_myo_rcca_myoX_load_df_time <- dplyr::bind_rows(
  iim_myo_rcca_myoX_topX_tab_time |> dplyr::mutate(type = "Time"),
  iim_myo_rcca_myoX_topY_tab_time |> dplyr::mutate(type = "Myokines")
)

iim_myo_rcca_myoX_fig_axis_time_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_load_df_time,
  ggplot2::aes(x = loading, y = reorder(var, loading), fill = type)
) +
  ggplot2::geom_col(alpha = 0.8) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::facet_wrap(~type, scales = "free_y") +
  ggplot2::labs(
    x = "Loading (Canonical component 1)",
    y = NULL,
    title = "Opposing poles of the latent time–myokine axis",
    subtitle = "Variables on opposite sides represent contrasting patterns"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_time_myokine_axis_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_axis_time_01)
grDevices::dev.off()

# --- Correlation with canonical axis (tiles)
iim_myo_rcca_myoX_cor_axis_time <- dplyr::bind_rows(
  tibble::tibble(
    var  = colnames(iim_myo_rcca_myoX_X_time),
    cor  = stats::cor(iim_myo_rcca_myoX_X_time, iim_myo_rcca_myoX_fit_time$variates$X[, 1]),
    type = "Time"
  ),
  tibble::tibble(
    var  = colnames(iim_myo_rcca_myoX_Y_time),
    cor  = stats::cor(iim_myo_rcca_myoX_Y_time, iim_myo_rcca_myoX_fit_time$variates$Y[, 1]),
    type = "Myokines"
  )
) |>
  dplyr::filter(abs(cor) > 0.15)

iim_myo_rcca_myoX_fig_cor_time_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_cor_axis_time,
  ggplot2::aes(x = type, y = reorder(var, cor), fill = cor)
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  ggplot2::labs(
    fill = "Correlation",
    title = "Correlation of individual variables with the latent CCA axis",
    subtitle = "Red vs blue indicate opposing patterns"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_canon_cor_time_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_cor_time_01)
grDevices::dev.off()

# --- Subtype and response visualization
iim_myo_rcca_myoX_plot_time <- iim_myo_rcca_myoX_d_time |>
  dplyr::transmute(
    immet_id, exam_order,
    subtype  = as.factor(disease_subtype),
    response = as.factor(response_m0_m6),
    compX1   = iim_myo_rcca_myoX_fit_time$variates$X[, 1],
    compY1   = iim_myo_rcca_myoX_fit_time$variates$Y[, 1]
  )

iim_myo_rcca_myoX_fig_response_time_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_plot_time,
  ggplot2::aes(compX1, compY1, color = response)
) +
  ggplot2::geom_point(size = 2.8, alpha = 0.85) +
  ggplot2::stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
  ggplot2::labs(
    x = "rCCA component 1 (Time)",
    y = "rCCA component 1 (Myokines)",
    title = "Patients projected onto the shared time–myokine space",
    subtitle = "Ellipses highlight response clustering (80% t-ellipse)"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_rcca_response_time_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_response_time_01)
grDevices::dev.off()

iim_myo_rcca_myoX_fig_subtype_time_01 <- ggplot2::ggplot(
  iim_myo_rcca_myoX_plot_time,
  ggplot2::aes(compX1, compY1, color = subtype)
) +
  ggplot2::geom_point(size = 2.8, alpha = 0.85) +
  ggplot2::stat_ellipse(type = "t", linetype = "dashed", level = 0.8, alpha = .5) +
  ggplot2::labs(
    x = "rCCA component 1 (Time)",
    y = "rCCA component 1 (Myokines)",
    title = "Patients projected onto the shared time–myokine space",
    subtitle = "Ellipses highlight subtype clustering (80% t-ellipse)"
  ) +
  ggplot2::theme_minimal(base_size = 13)

tiff(
  filename = file.path(iim_myo_rcca_myoX_outdir, paste0(iim_myo_rcca_myoX_prefix, "_rcca_subtype_time_01.tiff")),
  compression = "lzw", res = 150, width = 1600, height = 1200
)
print(iim_myo_rcca_myoX_fig_subtype_time_01)
grDevices::dev.off()


# **********************************************************************
# 10. výstup pro vyhodnocení ----
# **********************************************************************

# ==========================================================
# rCCA -> export tabulek do XLSX (multi-sheet) ----
# ==========================================================

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
  
  # ---- zajistit jména proměnných (mixOmics někdy nedá rownames)
  if (is.null(rownames(loadX))) rownames(loadX) <- colnames(X)
  if (is.null(rownames(loadY))) rownames(loadY) <- colnames(Y)
  
  # ---- efektivní počet komponent
  ncomp_eff <- min(
    ncomp,
    ncol(loadX), ncol(loadY),
    ncol(fit$variates$X), ncol(fit$variates$Y)
  )
  
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
  
  # ---- loadings long
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
  
  # ---- top proměnné pro cross-correlation (COMP1)
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
  
  if (is.null(dim(cor_raw))) {
    cor_mat <- matrix(
      cor_raw,
      nrow = length(top_x_vars),
      ncol = length(top_y_vars),
      dimnames = list(top_x_vars, top_y_vars)
    )
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
    c("subtype", "response", "sex", "disease_subtype", "response_m0_m6", "exam_order"),
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

# ==========================================================
# Export tabulek pro tento projekt (nekolidující názvy) ----
# ==========================================================

# META musí odpovídat řádkům X/Y (po complete.cases filtraci už sedí iim_myo_rcca_myoX_d_*)

iim_myo_rcca_myoX_meta_activ <- iim_myo_rcca_myoX_d_activ |>
  dplyr::select(immet_id, exam_order, dplyr::any_of(c("disease_subtype","response_m0_m6","sex"))) |>
  dplyr::mutate(
    exam_order = factor(tolower(as.character(exam_order)), levels = c("m0","m6"))
  )

iim_myo_rcca_myoX_meta_time <- iim_myo_rcca_myoX_d_time |>
  dplyr::select(immet_id, exam_order, dplyr::any_of(c("disease_subtype","response_m0_m6","sex"))) |>
  dplyr::mutate(
    exam_order = factor(tolower(as.character(exam_order)), levels = c("m0","m6"))
  )

# 1) AKTIVITA × MYOKINY
iim_myo_rcca_myoX_tabs_activ <- build_rcca_tables(
  fit = iim_myo_rcca_myoX_fit_activ,
  X   = iim_myo_rcca_myoX_X_activ,
  Y   = iim_myo_rcca_myoX_Y_activ,
  meta = iim_myo_rcca_myoX_meta_activ,
  analysis_name = "activity_vs_myokines_M0M6",
  block_x = "Activity",
  block_y = "Myokines",
  ncomp = 2,
  top_n_load_x = 25,
  top_n_load_y = 25,
  top_n_pairs_x = 15,
  top_n_pairs_y = 20,
  axis_cor_thr = 0.15
)

export_rcca_tables(
  iim_myo_rcca_myoX_tabs_activ,
  file.path("output/tables", paste0(format(Sys.Date(), "%y%m%d"), "_iim_myo_rcca_activity_vs_myokines_M0M6.xlsx"))
)

# 2) ČAS × MYOKINY
iim_myo_rcca_myoX_tabs_time <- build_rcca_tables(
  fit = iim_myo_rcca_myoX_fit_time,
  X   = iim_myo_rcca_myoX_X_time,
  Y   = iim_myo_rcca_myoX_Y_time,
  meta = iim_myo_rcca_myoX_meta_time,
  analysis_name = "time_vs_myokines_M0M6",
  block_x = "Time",
  block_y = "Myokines",
  ncomp = 2,
  top_n_load_x = 25,
  top_n_load_y = 25,
  top_n_pairs_x = 15,
  top_n_pairs_y = 20,
  axis_cor_thr = 0.15
)

export_rcca_tables(
  iim_myo_rcca_myoX_tabs_time,
  file.path("output/tables", paste0(format(Sys.Date(), "%y%m%d"), "_iim_myo_rcca_time_vs_myokines_M0M6.xlsx"))
)
