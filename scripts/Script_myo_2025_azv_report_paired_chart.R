# ************************************************************
# Requirements
# ************************************************************
# install.packages(c("dplyr","tidyr","purrr","ggplot2","sjPlot","gridExtra"))
# (grid is base R)

# ************************************************************
# Variable sets
# ************************************************************
fig_sel_01 <- c("enmo_full_recording_mean", "inactivity", "light",
                "moderate", "t400_time_min")

fig_sel_02 <- c("ig_gradient_enmo", "ig_intercept_enmo", "m5_enmo")

fig_sel_03 <- c("t100_b1m80_time_min", "t100_b5m80_time_min", "t100_b10m80_time_min")

fig_sel_04 <- c("p9931_enmo", "p9861_enmo", "p9792_enmo", "p9583_enmo", "p9167_enmo")

# ************************************************************
# 1) Data prep + FILTER to complete pairs (M0 & M6) per variable
# ************************************************************

fig_line_m0m6_respcol_build_df <- function(
    data,
    vars_sel,
    id_col = "immet_id",
    response_col = "response_m0_m6",
    exam_col = "exam_order",
    exam_levels = c("m0", "m6"),
    out_var_name = "var_name",
    out_var_value = "var_value",
    keep_complete_pairs = TRUE
) {
  
  df <- data |>
    dplyr::select(all_of(c(vars_sel, id_col, response_col, exam_col))) |>
    dplyr::mutate(
      "{exam_col}" := factor(tolower(as.character(.data[[exam_col]])), levels = exam_levels),
      "{response_col}" := factor(.data[[response_col]])
    ) |>
    tidyr::pivot_longer(
      cols = all_of(vars_sel),
      names_to = out_var_name,
      values_to = out_var_value
    ) |>
    dplyr::mutate("{out_var_name}" := factor(.data[[out_var_name]], levels = vars_sel))
  
  if (!keep_complete_pairs) return(df)
  
  # Keep only IDs that have BOTH M0 and M6 and NON-NA value, within each var_name
  df |>
    dplyr::group_by(.data[[out_var_name]], .data[[id_col]]) |>
    dplyr::filter(
      dplyr::n_distinct(.data[[exam_col]][!is.na(.data[[out_var_value]])]) == length(exam_levels),
      all(exam_levels %in% as.character(.data[[exam_col]][!is.na(.data[[out_var_value]])]))
    ) |>
    dplyr::ungroup()
}

# ************************************************************
# 2) One variable = one column (facet rows = response)
#    FIX: connect jittered points -> same jitter for line + point
# ************************************************************

fig_line_m0m6_respcol_make_var_column <- function(
    df_var,
    var_label,
    id_col = "immet_id",
    response_col = "response_m0_m6",
    exam_col = "exam_order",
    value_col = "var_value",
    color_values = c(
      "No Improvement" = "blue",
      "Minimal"        = "red",
      "Moderate"       = "darkgreen",
      "Major"          = "#9900FF"
    ),
    jitter_width = 0.08,
    jitter_seed  = 123
) {
  
  pos_jit <- ggplot2::position_jitter(
    width = jitter_width,
    height = 0,
    seed = jitter_seed
  )
  
  ggplot2::ggplot(
    df_var,
    ggplot2::aes(
      x = .data[[exam_col]],
      y = .data[[value_col]],
      col   = .data[[response_col]],
      shape = .data[[response_col]],
      group = .data[[id_col]]
    )
  ) +
    ggplot2::geom_line(
      position = pos_jit,
      linewidth = 0.9,
      alpha = 0.5
    ) +
    ggplot2::geom_point(
      position = pos_jit,
      size = 2,
      alpha = 0.5
    ) +
    ggplot2::scale_color_manual(values = color_values) +
    ggplot2::labs(title = var_label) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data[[response_col]]), scales = "fixed") +
    sjPlot::theme_sjplot2() +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title   = ggplot2::element_text(size = 11, face = "bold"),
      axis.text    = ggplot2::element_text(size = 10, face = "bold"),
      strip.text.y = ggplot2::element_text(size = 11, face = "bold"),
      legend.position = "none"
    )
}

# ************************************************************
# 3) Build all columns and assemble into one-row figure via gridExtra
# ************************************************************

fig_line_m0m6_respcol_grid_cols <- function(
    data,
    vars_sel,
    id_col = "immet_id",
    response_col = "response_m0_m6",
    exam_col = "exam_order",
    exam_levels = c("m0", "m6"),
    color_values = c(
      "No Improvement" = "blue",
      "Minimal"        = "red",
      "Moderate"       = "darkgreen",
      "Major"          = "#9900FF"
    ),
    jitter_width = 0.08,
    jitter_seed  = 123,
    keep_complete_pairs = TRUE
) {
  
  df <- fig_line_m0m6_respcol_build_df(
    data = data,
    vars_sel = vars_sel,
    id_col = id_col,
    response_col = response_col,
    exam_col = exam_col,
    exam_levels = exam_levels,
    keep_complete_pairs = keep_complete_pairs
  )
  
  plots_split <- df |>
    dplyr::group_by(var_name) |>
    dplyr::group_split()
  
  plots <- purrr::imap(plots_split, ~ fig_line_m0m6_respcol_make_var_column(
    df_var = .x,
    var_label = as.character(unique(.x$var_name)),
    id_col = id_col,
    response_col = response_col,
    exam_col = exam_col,
    value_col = "var_value",
    color_values = color_values,
    jitter_width = jitter_width,
    jitter_seed  = jitter_seed
  ))
  
  grob_row <- gridExtra::arrangeGrob(
    grobs = plots,
    nrow = 1
  )
  
  list(
    fig_line_m0m6_respcol_df    = df,
    fig_line_m0m6_respcol_plots = plots,
    fig_line_m0m6_respcol_grob  = grob_row
  )
}

# ************************************************************
# 4) Export helper for grob to TIFF (smoother, higher DPI)
# ************************************************************

save_line_grob_tiff <- function(
    grob,
    filename,
    ncols,
    res = 300,
    col_width_px = 520,
    height_px = 1400,
    compression = "lzw",
    use_cairo = TRUE
) {
  w <- max(1400, ncols * col_width_px)
  
  # "cairo" tends to look smoother (if available on your system)
  tiff(
    filename = filename,
    compression = compression,
    res = res,
    width = w,
    height = height_px,
    type = if (isTRUE(use_cairo)) "cairo" else "native"
  )
  
  grid::grid.newpage()
  grid::grid.draw(grob)
  dev.off()
}

# ************************************************************
# 5) Run for 4 sets + save 4 figures
# ************************************************************

fig_sets <- list(
  fig_sel_01 = fig_sel_01,
  fig_sel_02 = fig_sel_02,
  fig_sel_03 = fig_sel_03,
  fig_sel_04 = fig_sel_04
)

dir.create("output/figures/final", recursive = TRUE, showWarnings = FALSE)

out_list <- purrr::imap(fig_sets, ~{
  out <- fig_line_m0m6_respcol_grid_cols(
    data = d14_join,
    vars_sel = .x,
    jitter_width = 0.08,
    jitter_seed  = 123,
    keep_complete_pairs = TRUE
  )
  
  save_line_grob_tiff(
    grob = out$fig_line_m0m6_respcol_grob,
    filename = file.path(
      "output/figures/final",
      paste0(format(Sys.Date(), "%y%m%d"), "_slope_line_chart_", .y, ".tiff")
    ),
    ncols = length(.x),
    res = 300,
    col_width_px = 520,
    height_px = 1400,
    use_cairo = TRUE
  )
  
  out
})

# Optional: grobs to draw interactively
fig_line_01 <- out_list$fig_sel_01$fig_line_m0m6_respcol_grob
fig_line_02 <- out_list$fig_sel_02$fig_line_m0m6_respcol_grob
fig_line_03 <- out_list$fig_sel_03$fig_line_m0m6_respcol_grob
fig_line_04 <- out_list$fig_sel_04$fig_line_m0m6_respcol_grob

# Example viewer:
grid::grid.newpage(); grid::grid.draw(fig_line_01)


# Hezčí grafy s trendem -----
# ************************************************************
# Requirements
# ************************************************************
# install.packages(c("dplyr","tidyr","purrr","ggplot2","sjPlot","gridExtra"))

# ************************************************************
# Variable sets
# ************************************************************
fig_sel_01 <- c("enmo_full_recording_mean", "inactivity", "light",
                "moderate", "t400_time_min")

fig_sel_02 <- c("ig_gradient_enmo", "ig_intercept_enmo", "m5_enmo")

fig_sel_03 <- c("t100_b1m80_time_min", "t100_b5m80_time_min", "t100_b10m80_time_min")

fig_sel_04 <- c("p9931_enmo", "p9861_enmo", "p9792_enmo", "p9583_enmo", "p9167_enmo")

# ************************************************************
# 1) Data prep + keep only COMPLETE M0–M6 pairs
# ************************************************************

fig_line_m0m6_respcol_build_df <- function(
    data,
    vars_sel,
    id_col = "immet_id",
    response_col = "response_m0_m6",
    exam_col = "exam_order",
    exam_levels = c("m0", "m6"),
    out_var_name = "var_name",
    out_var_value = "var_value",
    keep_complete_pairs = TRUE
) {
  
  df <- data |>
    dplyr::select(all_of(c(vars_sel, id_col, response_col, exam_col))) |>
    dplyr::mutate(
      "{exam_col}" := factor(tolower(as.character(.data[[exam_col]])),
                             levels = exam_levels),
      "{response_col}" := factor(.data[[response_col]])
    ) |>
    tidyr::pivot_longer(
      cols = all_of(vars_sel),
      names_to = out_var_name,
      values_to = out_var_value
    ) |>
    dplyr::mutate("{out_var_name}" :=
                    factor(.data[[out_var_name]], levels = vars_sel))
  
  if (!keep_complete_pairs) return(df)
  
  df |>
    dplyr::group_by(.data[[out_var_name]], .data[[id_col]]) |>
    dplyr::filter(
      all(exam_levels %in% .data[[exam_col]][!is.na(.data[[out_var_value]])])
    ) |>
    dplyr::ungroup()
}

# ************************************************************
# 2) One variable = one column
#    + individual trajectories (jittered)
#    + LM line per response facet (NOT jittered)
# ************************************************************

fig_line_m0m6_respcol_make_var_column <- function(
    df_var,
    var_label,
    id_col = "immet_id",
    response_col = "response_m0_m6",
    exam_col = "exam_order",
    value_col = "var_value",
    color_values = c(
      "No Improvement" = "blue",
      "Minimal"        = "red",
      "Moderate"       = "darkgreen",
      "Major"          = "#9900FF"
    ),
    jitter_width = 0.08,
    jitter_seed  = 123
) {
  
  pos_jit <- ggplot2::position_jitter(
    width = jitter_width,
    height = 0,
    seed = jitter_seed
  )
  
  ggplot2::ggplot(
    df_var,
    ggplot2::aes(
      x = .data[[exam_col]],
      y = .data[[value_col]],
      group = .data[[id_col]]
    )
  ) +
    # ---- individual trajectories (visual only)
    ggplot2::geom_line(
      ggplot2::aes(color = .data[[response_col]]),
      position = pos_jit,
      linewidth = 0.8,
      alpha = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[[response_col]],
                   shape = .data[[response_col]]),
      position = pos_jit,
      size = 2,
      alpha = 0.5
    ) +
    # ---- LM trend (REAL DATA, no jitter)
    ggplot2::geom_smooth(
      ggplot2::aes(
        x = as.numeric(.data[[exam_col]]),
        y = .data[[value_col]],
        color = .data[[response_col]],
        group = .data[[response_col]]
      ),
      method = "lm",
      se = TRUE,
      linewidth = 1.2,
      alpha = 0.18,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_color_manual(values = color_values) +
    ggplot2::labs(title = var_label) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data[[response_col]]),
      scales = "fixed"
    ) +
    sjPlot::theme_sjplot2() +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title   = ggplot2::element_text(size = 11, face = "bold"),
      axis.text    = ggplot2::element_text(size = 10, face = "bold"),
      strip.text.y = ggplot2::element_text(size = 11, face = "bold"),
      legend.position = "none"
    )
}

# ************************************************************
# 3) Assemble columns with gridExtra
# ************************************************************

fig_line_m0m6_respcol_grid_cols <- function(
    data,
    vars_sel,
    jitter_width = 0.08,
    jitter_seed  = 123
) {
  
  df <- fig_line_m0m6_respcol_build_df(
    data = data,
    vars_sel = vars_sel,
    keep_complete_pairs = TRUE
  )
  
  plots <- df |>
    dplyr::group_by(var_name) |>
    dplyr::group_split() |>
    purrr::imap(~ fig_line_m0m6_respcol_make_var_column(
      df_var = .x,
      var_label = as.character(unique(.x$var_name)),
      jitter_width = jitter_width,
      jitter_seed  = jitter_seed
    ))
  
  gridExtra::arrangeGrob(grobs = plots, nrow = 1)
}

# ************************************************************
# 4) Export helper (high-quality TIFF)
# ************************************************************

save_line_grob_tiff <- function(
    grob,
    filename,
    ncols,
    res = 300,
    col_width_px = 520,
    height_px = 1400
) {
  w <- max(1400, ncols * col_width_px)
  
  tiff(
    filename = filename,
    compression = "lzw",
    res = res,
    width = w,
    height = height_px,
    type = "cairo"
  )
  grid::grid.newpage()
  grid::grid.draw(grob)
  dev.off()
}

# ************************************************************
# 5) Run & save all 4 figures
# ************************************************************

fig_sets <- list(
  fig_sel_01 = fig_sel_01,
  fig_sel_02 = fig_sel_02,
  fig_sel_03 = fig_sel_03,
  fig_sel_04 = fig_sel_04
)

dir.create("output/figures/final", recursive = TRUE, showWarnings = FALSE)

purrr::imap(fig_sets, ~{
  g <- fig_line_m0m6_respcol_grid_cols(
    data = d14_join,
    vars_sel = .x
  )
  
  save_line_grob_tiff(
    grob = g,
    filename = file.path(
      "output/figures/final",
      paste0(format(Sys.Date(), "%y%m%d"), "_slope_line_chart_", .y, ".tiff")
    ),
    ncols = length(.x)
  )
})
