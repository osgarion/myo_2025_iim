##############################################
# PLS-DA PIPELINE WITH AUTOMATIC ELLIPSE CHECK
# AND SAFE DROPLEVELS HANDLING
# Author: Jiří Baloun
# Date: 13. 11. 2025
##############################################

library(mixOmics)
library(dplyr)
library(purrr)

################################################
# Y variables checked by PLS-DA
################################################

grp_vars <- c(
  "podtyp_nemoci_zjednoduseny",
  "mmt8_categories",
  "mda_categories",
  "odpoved_na_terapii_m0_vs_m6"
)

################################################
# AUTOMATIC ELLIPSE VALIDATION
################################################

ellipse_allowed <- function(Y) {
  group_sizes <- table(Y)
  
  if (any(group_sizes <= 2)) {
    message("⚠️  Ellipse disabled: some groups have ≤ 2 samples (",
            paste(names(group_sizes[group_sizes <= 2]), collapse = ", "),
            ")")
    return(FALSE)
  } else {
    return(TRUE)
  }
}

################################################
# Run ONE PLS-DA model
################################################

run_plsda <- function(df,
                      yvar,
                      prefix,
                      outdir = "output/figures/cluster analyses/",
                      ncomp = 2) {
  
  message("\n==============================")
  message("Running PLS-DA for: ", yvar)
  message("==============================\n")
  
  # STEP 1 — keep rows with non-NA Y
  df2 <- df |> filter(!is.na(.data[[yvar]]))
  
  # STEP 2 — remove Y and _total columns from predictors
  predictors <- df2 |> 
    select(-all_of(yvar)) |>
    select(-contains("_total"))
  
  # STEP 3 — remove incomplete predictor rows
  complete_idx <- complete.cases(predictors)
  predictors <- predictors[complete_idx, ]
  Y <- df2[[yvar]][complete_idx]
  
  # NEW FIX: remove unused factor levels
  Y <- droplevels(Y)
  
  # STEP 4 — create model matrix
  X <- model.matrix(~ . - 1, data = predictors)
  
  # SAFETY CHECK
  if (nrow(X) != length(Y)) {
    stop("X/Y mismatch: X=", nrow(X), " Y=", length(Y))
  }
  
  # STEP 5 — PLSDA
  pls_model <- plsda(X, Y, ncomp = ncomp, scale = TRUE)
  
  # STEP 6 — check ellipse validity
  ellipse_flag <- ellipse_allowed(Y)
  
  # STEP 7 — output paths
  group_file  <- paste0(outdir, prefix, "_", yvar, "_group.tiff")
  load_file   <- paste0(outdir, prefix, "_", yvar, "_loadings.tiff")
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # STEP 8 — group plot
  tiff(group_file, compression="lzw", res=300, width=1600, height=1200)
  
  plotIndiv(
    pls_model,
    comp = c(1, 2),
    group = Y,
    ellipse = ellipse_flag,
    legend = TRUE,
    pch = 19,                  # FIX: disable shape scale → avoid ggplot bug
    ind.names = FALSE,
    title = paste0("PLS-DA — ", yvar)
  )
  
  dev.off()
  
  # STEP 9 — loadings plot
  tiff(load_file, compression="lzw", res=300, width=1600, height=1200)
  
  plotLoadings(
    pls_model,
    contrib = "max",
    method = "median",
    title = paste0("Loadings — ", yvar)
  )
  
  dev.off()
  
  return(pls_model)
}

################################################
# Run ALL PLS-DA models
################################################

run_all_plsda <- function(df,
                          prefix = "pls_m0",
                          outdir = "output/figures/cluster analyses/",
                          ncomp = 2) {
  
  message("\n>>> Starting full PLS-DA pipeline...\n")
  
  results <- map(
    grp_vars,
    ~ run_plsda(df, .x, prefix, outdir, ncomp)
  )
  
  names(results) <- grp_vars
  
  message("\n>>> Pipeline finished.\n")
  
  return(results)
}

################################################
# END OF SCRIPT
################################################

# Usage:
# pls_m0_results    <- run_all_plsda(d08_m0_pls,    prefix="251113_m0_pls")
# pls_delta_results <- run_all_plsda(d08_delta_pls, prefix="251113_delta_pls")
