









# 
# 
# mod_mixMODmmt(d04_sel1_nested$data[[5]], 100)
# 
# data <- d04_sel1_nested$data[[5]]
# 
# 
# mod_mixMODmmt(d04_sel1_nested$data[[1]], 80)
# 
# 
# 
# data <- d04_sel1_nested$data[[5]]
# upper = 100







# fig_mixMODmmt <- function(data, mod, var_indep_name, var_dep_name, upper) {
#   data_set <- data |> 
#     mutate(
#       var_indep_value_scl = log(var_indep_value),
#       y_capped = pmin(var_dep_value, upper),
#       ind = as.integer(var_dep_value >= upper)
#     )
#   
#   pred_df <- GLMMadaptive::effectPlotData(mod, newdata = na.omit(data_set))
#   
#   fig <- ggplot(pred_df, aes(x = var_indep_value, y = var_dep_value)) +
#     geom_point(aes(color = poradie_vysetrenia), alpha = 0.6, show.legend = FALSE) +
#     geom_line(aes(y = pred, color = poradie_vysetrenia), size = 1, show.legend = FALSE) +
#     geom_ribbon(aes(ymax = upp, ymin = low, fill = poradie_vysetrenia), 
#                 alpha = 0.3, linetype = 0) +
#     facet_wrap(~ poradie_vysetrenia, scales = "free_y") +
#     scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
#     labs(
#       x = var_indep_name,
#       y = var_dep_name,
#       fill = "Time point"
#     ) +
#     theme_sjplot2() +
#     theme(
#       axis.title = element_text(size = 14, face = "bold"),
#       strip.text = element_text(face = "bold"),
#       strip.background = element_rect(fill = "gray90", color = NA)
#     )
#   
#   return(fig)
# }
# 
