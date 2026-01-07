##############################################
# PLS-DA PIPELINE — AUTO Y-VAR DETECTION,
# VECTOR GRAPHICS AND VARIABLE RENAME
##############################################

library(mixOmics)
library(dplyr)
library(purrr)

################################################
# POSSIBLE Y VARIABLES
################################################

grp_vars_all <- c(
  "mmt8_categories",
  "mda_categories"
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
  
  # STEP 2 — remove Y and apply rename_vec
  predictors <- df2 |> 
    select(-all_of(yvar)) |>
    rename(any_of(rename_vec))   # must exist in environment
  
  # STEP 3 — remove incomplete predictor rows
  complete_idx <- complete.cases(predictors)
  predictors <- predictors[complete_idx, ]
  Y <- droplevels(df2[[yvar]][complete_idx])
  
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
  
  # STEP 7 — output paths (SVG — vector graphics)
  group_file  <- paste0(outdir, prefix, "_", yvar, "_group.svg")
  load_file   <- paste0(outdir, prefix, "_", yvar, "_loadings.svg")
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # STEP 8 — group plot (SVG)
  svg(group_file, width = 12, height = 9)
  plotIndiv(
    pls_model,
    comp = c(1, 2),
    group = Y,
    ellipse = ellipse_flag,
    legend = TRUE,
    pch = 19,
    ind.names = FALSE,
    title = paste0("PLS-DA — ", yvar)
  )
  dev.off()
  
  # STEP 9 — loadings plot (SVG)
  svg(load_file, width = 12, height = 9)
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
# Run ALL PLS-DA models — ONLY EXISTING Y-VARS
################################################

run_all_plsda <- function(df,
                          prefix = "pls_m0",
                          outdir = "output/figures/cluster analyses/",
                          ncomp = 2) {
  
  message("\n>>> Starting full PLS-DA pipeline...\n")
  
  # Detect which Y-vars exist in the table
  available_y <- grp_vars_all[grp_vars_all %in% colnames(df)]
  
  if (length(available_y) == 0) {
    stop("❌ No Y variables found (expected mmt8_categories or mda_categories).")
  }
  
  message("Found Y variables: ", paste(available_y, collapse = ", "))
  
  results <- map(
    available_y,
    ~ run_plsda(df, .x, prefix, outdir, ncomp)
  )
  
  names(results) <- available_y
  
  message("\n>>> Pipeline finished.\n")
  
  return(results)
}

################################################
# END OF SCRIPT
################################################
