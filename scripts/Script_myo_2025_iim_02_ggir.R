# 1. Libraries & functions ----
# *******************************************************
# required
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, purrr, conflicted, renv)
# Functions and libraries uploading
list.files("scripts/functions/", pattern="*.*", full.names=TRUE) |> map(~source(.))
renv::status()

# 2. Environment reproducibility ----
# *******************************************************
# RENV
# renv::init()        # once
# renv::status()      # sync?
# renv::snapshot()    # after changing deps
# renv::restore()     # on a new machine

# Back-up
back_up("scripts/functions/FUN_01.R") # the the destination subdirectory specify using 'path_dest'
back_up("scripts/functions/OBJ_01.R") # the the destination subdirectory specify using 'path_dest'
back_up("scripts/Script_myo_2025_iim_01.R") # the the destination subdirectory specify using 'path_dest'

# Save to .RData
save(xx1,
     file = "reports/markD_01.RData")
save.image("reports/markD_01.RData")
# Save the tables into data/tables.RData by listing them individually
cgwtools::resave(xx2,xx3, file = "reports/markD_01.RData") # resave a list of tables that I'll use in the .Rmd file.
# Save the tables into data/tables.RData using "patterns" 
cgwtools::resave(list=ls(pattern="tbl"), file = "reports/markD_01.RData")

# *******************************************************
# 3. EDA ----
# *******************************************************

## figures ----
### ostatni ----
d10_ostatni |>
  ggplot(aes(`AD_ig_intercept_ENMO_0-24hr`, `AD_ig_gradient_ENMO_0-24hr`)) +
  geom_point(shape = 15, col = "darkblue") +
  geom_smooth(
    method = "loess", se = FALSE,
    linetype = "dashed", col = "darkblue"
  ) +
  geom_point(
    data = d10_ostatni |> filter(`AD_ig_rsquared_ENMO_0-24hr` > 0.9),
    shape = 17, col = "darkgreen"
  ) +
  geom_smooth(
    data = d10_ostatni |> filter(`AD_ig_rsquared_ENMO_0-24hr` > 0.9),
    method = "loess", se = FALSE,
    linetype = "dashed", col = "darkgreen"
  ) +
  labs(
    title = "Intercept vs. gradient (ENMO 0–24 h)",
    subtitle = "Dashed LOESS fit for all samples (blue) and high-fit subset R² > 0.9 (green)",
    x = "Intercept (ENMO, 0–24 h)",
    y = "Gradient (ENMO, 0–24 h)"
  ) +
  sjPlot::theme_sjplot2() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(size = 15, face = "italic"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, face = "bold")
  )

## q-q plots ----
# *******************************************************
### čas ----
# *******************************************************
d10_cas |> 
  filter(exam %in% c("M0", "M6")) |>
  pivot_longer(cols = where(is.numeric),
               names_to = "var",
               values_to = "values") |> 
  ggqqplot(x = "values",
           color = "exam",
           palette = c("#00AFBB", "#E7B800"),
           facet.by = "var", scales = "free")





## ggally ----
# *******************************************************
### čas ----
# *******************************************************
d10_cas |> 
  filter(exam %in% c("M0", "M6")) |>
  (\(df)
   ggpairs(
     df,
     columns = 3:ncol(df),
     mapping = ggplot2::aes(color = exam),
     lower = list(continuous = my_fn),
     upper = list(continuous = GGally::wrap("cor", method = "pearson"))
   )
  )()

### aktivita ----
# *******************************************************
d10_aktivita |> 
  filter(exam %in% c("M0", "M6")) |>
  (\(df)
   ggpairs(
     df,
     columns = 3:ncol(df),
     mapping = ggplot2::aes(color = exam),
     lower = list(continuous = my_fn),
     upper = list(continuous = GGally::wrap("cor", method = "pearson"))
   )
  )()

# ***********************************************************************
d10_cas |> 
  filter(exam %in% c("M0", "M6")) |>
  mutate(
    across(
      where(is.numeric),
      ~ as.numeric(scale(.x)))) |> 
  
  ## missingness ----
# *******************************************************
vis_dat(d11_clinics)
vis_miss(d11_clinics |> 
           filter(exam_order %in% c("M0", "M6")) |> 
           select(
             where(~ mean(is.na(.x)) <= 0.20)
           ))





