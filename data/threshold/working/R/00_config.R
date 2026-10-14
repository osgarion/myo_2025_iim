# ======================================================================
# Script:  00_config.R
# Project: MYO-II / Paper A — ENMO intensity-spectrum (thresholds) analysis
# Author:  Veronika Balounová (prepared with Claude Opus 5.5)
# Date:    2026-10-13
# ======================================================================
#
# PURPOSE
# ---
# Shared configuration for every script in `working/`: locates the working
# directory and the master RDS on any machine, checks packages, defines
# locked analysis constants and output folders. Sourced, never run alone.
#
# INPUTS
# ---
#   MYOII_ACC_OUTCOMES_MASTER_V1.rds — canonical master (see README.md for lookup order)
#
# ANALYSES
# ---
#   none (configuration only)
#
# OUTPUTS
# ---
#   Objects: WORKING_DIR, INPUT_RDS, PATHS, constants (MAIN_VISITS, B_BOOT, SEED, ...)
#
# ======================================================================

# **********************************************************************
# 1. Working directory ----
# **********************************************************************
# Lookup order: env MYOII_WORKING_DIR -> folder of the active RStudio document
# (scripts live in working/, helpers in working/R/) -> current working directory.
locate_working_dir <- function() {
  env <- Sys.getenv("MYOII_WORKING_DIR")
  if (nzchar(env)) return(normalizePath(env, winslash = "/", mustWork = TRUE))
  is_working <- function(p) file.exists(file.path(p, "R", "00_config.R"))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) {
      cand <- c(dirname(p), dirname(dirname(p)))
      hit  <- cand[vapply(cand, is_working, logical(1))]
      if (length(hit)) return(normalizePath(hit[1], winslash = "/"))
    }
  }
  cand <- c(getwd(), file.path(getwd(), "working"), file.path(getwd(), "data", "threshold", "working"), dirname(getwd()))
  hit  <- cand[vapply(cand, is_working, logical(1))]
  if (!length(hit)) stop("Cannot locate the 'working' folder. Set MYOII_WORKING_DIR ",
                         "or setwd() to the folder that contains R/00_config.R.")
  normalizePath(hit[1], winslash = "/")
}
if (!exists("WORKING_DIR") || !file.exists(file.path(WORKING_DIR, "R", "00_config.R"))) WORKING_DIR <- locate_working_dir()

# **********************************************************************
# 2. Packages ----
# **********************************************************************
# No silent installation: missing packages stop with an explicit install command.
PKG_ANALYSIS <- c("dplyr", "tidyr", "purrr", "tibble", "stringr", "ggplot2", "scales",
                  "patchwork", "lme4", "glmmTMB", "future", "furrr", "rio", "writexl")
PKG_REPORT   <- c("rmarkdown", "knitr", "kableExtra", "sessioninfo")

check_packages <- function(pkgs = PKG_ANALYSIS) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages. Install with:\n  install.packages(c(\"",
                         paste(miss, collapse = "\", \""), "\"))", call. = FALSE)
  invisible(TRUE)
}
check_packages()
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(ggplot2)
})

# **********************************************************************
# 3. Master RDS ----
# **********************************************************************
# Lookup order: env MYOII_MASTER_RDS -> working/data/ -> parent folder of working/
# (= data/threshold/ in the myo_2025_iim project).
RDS_NAME <- "MYOII_ACC_OUTCOMES_MASTER_V1.rds"
locate_master <- function() {
  env  <- Sys.getenv("MYOII_MASTER_RDS")
  cand <- c(if (nzchar(env)) env,
            file.path(WORKING_DIR, "data", RDS_NAME),
            file.path(dirname(WORKING_DIR), RDS_NAME))
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("Master RDS not found. Tried:\n  ", paste(cand, collapse = "\n  "),
                         "\nSet Sys.setenv(MYOII_MASTER_RDS = '<path>') before sourcing.", call. = FALSE)
  normalizePath(hit[1], winslash = "/")
}
INPUT_RDS <- locate_master()
INPUT_MD5 <- unname(tools::md5sum(INPUT_RDS))

# **********************************************************************
# 4. Locked analysis constants ----
# **********************************************************************
MAIN_VISITS        <- c("M0", "M3", "M6", "M18")
FOLLOWUP_VISITS    <- c("M3", "M6", "M18")
EXPLORATORY_VISITS <- c("M30", "M48", "M60")
THRESHOLD_GRID_MG  <- c(20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 140, 160, 180,
                        200, 220, 240, 260, 280, 300, 325, 350, 375, 400)

# Subtype coding (confirmed by user 2026-10-13): 0 and 4 = DM, 2 and 3 = IMNM, 5 = ASyS.
# PM (1) and "UCTD/SSc?" are excluded from subtype analyses only.
# All patients stay in the main analyses.
SUBTYPE_LABELS <- c("0" = "DM", "4" = "DM", "2" = "IMNM", "3" = "IMNM", "5" = "ASyS")

# Coverage limits (min/day per 10 mg for bins, min/day above T for cumulative).
# Descriptive and pragmatic, not clinically validated.
COVERAGE_LIMITS  <- c(0, 0.1, 0.25, 0.5, 1)
COVERAGE_PRIMARY <- 0.25

# Bootstrap (patient-level cluster bootstrap).
# B can be lowered for a quick test run: Sys.setenv(MYOII_B = 50) before sourcing.
B_BOOT <- as.integer(Sys.getenv("MYOII_B", unset = "2000"))
SEED   <- 20261013L
N_WORKERS <- local({
  n <- suppressWarnings(parallel::detectCores(logical = TRUE))
  n <- if (is.na(n)) 1L else max(1L, n - 1L)
  cap <- getOption("myo_ncores", n)          # project-level cap (.Rprofile), if set
  as.integer(max(1L, min(n, cap)))
})

# Visit encoding for all figures (B/W printable).
VISIT_LINETYPES <- c(M0 = "dotdash", M3 = "dotted", M6 = "dashed", M18 = "solid")
VISIT_SHAPES    <- c(M0 = 1, M3 = 2, M6 = 0, M18 = 16)

# **********************************************************************
# 5. Output folders ----
# **********************************************************************
PATHS <- list(
  rdata   = file.path(WORKING_DIR, "output", "RData"),
  tables  = file.path(WORKING_DIR, "output", "tables"),
  figures = file.path(WORKING_DIR, "output", "figures"),
  logs    = file.path(WORKING_DIR, "output", "logs")
)
invisible(lapply(PATHS, dir.create, recursive = TRUE, showWarnings = FALSE))

date_prefix <- function() format(Sys.Date(), "%y%m%d")

# **********************************************************************
# 6. Saving result objects without patient-level rows ----
# **********************************************************************
# ggplot objects (plot_env, aes() quosures, layer internals) keep references to
# the environment they were built in. Scripts run via source(local = new.env()),
# so that environment -- with the patient-level data -- would be serialised into
# every RData holding a figure. save_public() saves a private copy in which all
# reachable non-package environments are stripped of bindings that hold rows
# with a patient_id column; it refuses to write if any such rows remain.
holds_patient_rows <- function(v, depth = 0L) {
  if (is.data.frame(v)) return("patient_id" %in% names(v))
  if (is.list(v) && !is.environment(v) && depth < 6L) {
    return(any(vapply(v, holds_patient_rows, logical(1), depth = depth + 1L)))
  }
  FALSE
}
walk_envs <- function(x, on_env, seen = new.env(), depth = 0L) {
  if (depth > 60L || is.null(x)) return(invisible())
  if (is.environment(x)) {
    if (isNamespace(x) || identical(x, globalenv()) || identical(x, baseenv()) ||
        identical(x, emptyenv()) || grepl("^package:", environmentName(x))) return(invisible())
    key <- format(x)
    if (!is.null(seen[[key]])) return(invisible())
    seen[[key]] <- TRUE
    on_env(x)
    for (n in ls(x, all.names = TRUE)) {
      tryCatch(walk_envs(get(n, envir = x), on_env, seen, depth + 1L), error = function(e) NULL)
    }
    return(walk_envs(parent.env(x), on_env, seen, depth + 1L))
  }
  if (is.function(x)) return(walk_envs(environment(x), on_env, seen, depth + 1L))
  if (is.list(x)) for (el in x) walk_envs(el, on_env, seen, depth + 1L)
  at <- attributes(x)
  for (a in setdiff(names(at), c("names", "row.names", "class", "dim", "dimnames", "levels"))) {
    walk_envs(at[[a]], on_env, seen, depth + 1L)
  }
  invisible()
}
scrub_patient_envs <- function(x) {
  x <- unserialize(serialize(x, NULL))   # private copy; the live session is untouched
  walk_envs(x, function(e) {
    drop <- Filter(function(n) tryCatch(holds_patient_rows(get(n, envir = e)), error = function(err) FALSE),
                   ls(e, all.names = TRUE))
    if (length(drop)) rm(list = drop, envir = e)
  })
  x
}
check_no_patient_envs <- function(x) {
  left <- character()
  walk_envs(x, function(e) {
    hit <- Filter(function(n) tryCatch(holds_patient_rows(get(n, envir = e)), error = function(err) FALSE),
                  ls(e, all.names = TRUE))
    left <<- c(left, hit)
  })
  unique(left)
}
save_public <- function(obj_name, file, envir = parent.frame()) {
  x <- scrub_patient_envs(get(obj_name, envir = envir))
  left <- c(if (holds_patient_rows(x)) obj_name, check_no_patient_envs(x))
  if (length(left)) stop("save_public(): patient-level objects still reachable: ", paste(left, collapse = ", "), call. = FALSE)
  e <- new.env(); assign(obj_name, x, envir = e)
  save(list = obj_name, file = file, envir = e)
  invisible(file)
}

message("[config] WORKING_DIR = ", WORKING_DIR)
message("[config] INPUT_RDS   = ", INPUT_RDS, " (md5 ", INPUT_MD5, ")")
message("[config] workers     = ", N_WORKERS, ", B = ", B_BOOT, ", seed = ", SEED)
