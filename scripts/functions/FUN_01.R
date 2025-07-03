# Libraries & functions ----
## libraries ----
pacman::p_load(update = T,  
               equatiomatic, beepr, tictoc, 
               tidyverse, purrr, furrr, easystats, rio, janitor, ggthemes, car,
               gtsummary, skimr, sjPlot, flextable, ggpubr, rstatix, tidymodels,
               kableExtra, skimr, GGally, testthat, factoextra, gplots, uwot
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
  lmerTest::lmer
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

## Possibly ----
possLMER <- possibly(.f=lmer, otherwise = NA_real_)


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









# tables ----
## column with 0, 1, and NA ----
is_01_col <- function(x) {
  all(unique(x) %in% c(0, 1, NA))
}

## skimr ----
my_skim <- skimr::skim_with(numeric = sfl(median = ~ median(., na.rm = TRUE),
                                          mad = ~ mad(., na.rm = TRUE)), 
                            append = T)


file_to_load <- function(folder_path, phrase) {
  # Získání všech souborů (včetně cesty)
  files <- list.files(folder_path, recursive = TRUE, full.names = TRUE)
  
  # Filtrování souborů podle fráze v názvu
  matching_files <- files[grepl(phrase, basename(files), ignore.case = TRUE)]
  
  return(matching_files)
}

