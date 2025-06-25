# Header 1 ----
## Header 2 ----
### Header 3 ----

# References ----
# model diagnosis
## https://easystats.github.io/see/articles/performance.html
## https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
# statistical data and model report
## report(model) # report(model)
## https://github.com/easystats/report
# Multivariate analyses
## https://cran.r-project.org/web/packages/MVN/vignettes/MVN.pdf
# Regression model summary
## https://strengejacke.github.io/sjPlot/articles/tab_model_estimates.html
## https://cran.r-project.org/web/packages/jtools/vignettes/summ.html
# Regression model plots
## https://cran.r-project.org/web/packages/sjPlot/vignettes/plot_model_estimates.html
# Markdown session info
## https://stackoverflow.com/questions/33156946/how-can-i-format-sessioninfo-in-rmarkdown

# Directories ----
dir.create("data/raw", recursive = TRUE)
dir.create("data/processed", recursive = TRUE)
dir.create("scripts/functions", recursive = TRUE)
dir.create("output/figures", recursive = TRUE)
dir.create("output/tables", recursive = TRUE)
dir.create("misc")
dir.create("reports")



# Libraries & functions ----
# required
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, purrr,conflicted)

# load libraries
pacman::p_load(
  beepr, tictoc,
  tidyverse, readxl, data.table, magrittr, Rmisc, purrr, easystats, rio, janitor,
  ggthemes, car, labelled, testthat
)

# Missing values, multivariate analyses
pacman::p_load(naniar, MVN) 

# Functions and libraries uploading
list.files("scripts/functions/", pattern="*.*", full.names=TRUE) |> map(~source(.))

# Back-up
back_up("scripts/functions/FUN_01.R") # the the destination subdirectory specify using 'path_dest'
back_up("scripts/functions/OBJ_01.R") # the the destination subdirectory specify using 'path_dest'
back_up("scripts/Script_myo_2025_iim_01.R") # the the destination subdirectory specify using 'path_dest'

# Markdown ----
options(knitr.kable.NA = '',     # empty space in cells with NAs in kable
        scipen = 999)            # non-academic format of numbers
# Save to .RData
save(xx1,
     file = "reports/markD_01.RData")
save.image("reports/markD_01.RData")
# Save the tables into data/tables.RData by listing them individually
cgwtools::resave(xx2,xx3, file = "reports/markD_01.RData") # resave a list of tables that I'll use in the .Rmd file.
# Save the tables into data/tables.RData using "patterns" 
cgwtools::resave(list=ls(pattern="tbl"), file = "reports/markD_01.RData")

# copy-to-clipboard botton
# ```{r xaringanExtra-clipboard, echo=FALSE}
# xaringanExtra::use_clipboard()
# ```

# EDA ----
## Duplicates ----
eval_dup_01 <- d03 |> get_dupes()  |>  capture_messages()

## Missing values ----
vis_miss(d03)
gg_miss_var(d03)
gg_miss_var(d03 ,facet=poradie_vysetrenia)
gg_miss_case(d03)
gg_miss_case(d03 ,facet = poradie_vysetrenia)
(miss_tab01 <- miss_var_summary(d03))
(miss_tab02 <- miss_case_summary(d03))

gg_miss_upset(d03) # give an overall pattern of missingness
gg_miss_fct(x = d03, fct = poradie_vysetrenia)

## overview ----
d03 |>
  select(-projekt_id) |> 
  my_skim()

## correlation ----
### ggally ----
ggpairs(d03,
        columns = d03 |> select(where(is.numeric)) |> names(), 
        lower = list(continuous = my_fn)) +
  theme_sjplot2()

# Statistics ----



