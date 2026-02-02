# **********************************************************************
# 1. preparation ----
# **********************************************************************
# **********************************************************************
## libraries ----
# **********************************************************************
library(lme4)
library(broom.mixed)
library(recipes)


# **********************************************************************
## data ----
# **********************************************************************
d14_join <- d10_aktivita |>
  transmute(
    immet_id   = as.character(sample),
    exam_order = tolower(as.character(exam)),
    across(where(is.numeric), identity)
  ) |>
  inner_join(
    d10_cas |>
      transmute(
        immet_id   = as.character(sample),
        exam_order = tolower(as.character(exam)),
        across(where(is.numeric), identity)
      ),
    by = c("immet_id", "exam_order"),
    suffix = c("_activ", "_time")
  ) |>
  inner_join(
    d11_clinics_imp |>
      mutate(exam_order = tolower(as.character(exam_order))) |>
      mutate(immet_id = as.character(immet_id)),
    by = c("immet_id", "exam_order")
  ) |>
  filter(exam_order %in% c("m0", "m6")) |>
  distinct()

age_m0 <- d14_join |>
  filter(exam_order == "m0") |>
  select(immet_id, age) |>
  rename(age_m0 = age)

d14_join <- d14_join |>
  left_join(age_m0, by = "immet_id")

# **********************************************************************
# 2. Analyses ----
# **********************************************************************
## Clinics effect on accelerometry ----
# **********************************************************************
### Latent score ----
# **********************************************************************

x_activ <- c(
  "enmo_full_recording_mean","acc_day_spt_wei","m5_enmo","m1_enmo",
  "p9167_enmo","p9375_enmo","p9583_enmo","p9688_enmo","p9792_enmo",
  "p9861_enmo","p9931_enmo","ig_gradient_enmo","ig_intercept_enmo"
)

x_time <- c(
  "m5hr_time","t20_time_min","t20_b1m80_time_min","t20_b5m80_time_min","t20_b10m80_time_min",
  "t30_time_min","t30_b1m80_time_min","t30_b5m80_time_min","t30_b10m80_time_min",
  "t40_time_min","t40_b1m80_time_min","t40_b5m80_time_min","t40_b10m80_time_min",
  "t100_time_min","t100_b1m80_time_min","t100_b5m80_time_min","t100_b10m80_time_min",
  "moderate", "light", "inactivity"
)

prep_scores <- function(df, cols, prefix){
  rec <- recipe(~ ., data = df |> select(all_of(cols))) |>
    step_impute_knn(all_predictors()) |>
    step_zv(all_predictors()) |>
    step_corr(all_predictors(), threshold = 0.9) |>
    step_normalize(all_predictors()) |>
    step_pca(all_predictors(), num_comp = 2) |>
    prep()
  
  scores <- bake(rec, new_data = NULL) |>
    transmute(
      !!paste0(prefix, "_pc1") := PC1,
      !!paste0(prefix, "_pc2") := PC2
    )
  list(recipe = rec, scores = scores)
}

tmp_activ <- prep_scores(d14_join, x_activ, "activ")
tmp_time  <- prep_scores(d14_join, x_time,  "time")

d15_scores <- d14_join |>
  bind_cols(tmp_activ$scores, tmp_time$scores)


clin_num <- clin_cols |> keep(~ is.numeric(d11_clinics_imp[[.x]]))
clin_fac <- clin_cols |> keep(~ is.factor(d11_clinics_imp[[.x]]) || is.character(d11_clinics_imp[[.x]]) || is.logical(d11_clinics_imp[[.x]]))

rec_scaled <- recipe(~ ., data = d15_scores) |>
  step_impute_knn(all_of(clin_num), "age_m0") |>
  step_normalize(all_of(clin_num), "age_m0") |>
  prep()

d15_scaled <- bake(rec_scaled, new_data = d15_scores)

# příklad: klinická proměnná physician_vas, adjustace sex + bmi + age_m0
# fit1 <- lmer(
#   activ_pc1 ~ exam_order + physician_vas + exam_order:physician_vas +
#     sex + bmi + age_m0 + (1 | immet_id),
#   data = d15_scaled
# )
# 
# tidy(fit1, effects = "fixed", conf.int = TRUE)

pc_cols <- c("activ_pc1", "time_pc1")
clin_exclude <- c("immet_id", "exam_order")  
clin_cols <- setdiff(names(d11_clinics_imp), c(clin_exclude, pc_cols))
base_covars <- c("bmi", "age_m0") |> intersect(names(d15_scores))

d15_scaled <- bake(rec_scaled, new_data = d15_scores) |> 
  mutate(exam_order = factor(tolower(as.character(exam_order)), levels = c("m0", "m6"))) |>
  mutate(immet_id = as.character(immet_id))

pc_long <- d15_scaled |>
  select(immet_id, exam_order, all_of(pc_cols)) |>
  pivot_longer(cols = all_of(pc_cols), names_to = "pc_name", values_to = "pc_value")


clin_long_num <- d15_scaled |>
  select(immet_id, exam_order, disease_subtype, all_of(clin_num), all_of(base_covars)) |>
  pivot_longer(cols = setdiff(clin_num, base_covars), 
               names_to = "clin_name", values_to = "clin_value") |>
  mutate(clin_type = "numeric",
         clin_value = as.numeric(clin_value)) |> 
  full_join(pc_long, by = c("immet_id", "exam_order")) |> 
  group_by(clin_name, pc_name) |> 
  nest() |> 
  ungroup()

clin_long_fac <- d15_scaled |>
  select(immet_id, exam_order, disease_subtype, response_m0_m6, all_of(clin_fac), all_of(base_covars)) |>
  mutate(
    response_m0_m6   = factor(response_m0_m6),
    disease_subtype  = factor(disease_subtype)
  ) |>
  pivot_longer(
    cols = -c(immet_id, exam_order, disease_subtype, response_m0_m6, all_of(base_covars)),
    names_to = "clin_name",
    values_to = "clin_value"
  ) |>
  mutate(
    clin_type  = "factor",
    clin_value = factor(clin_value)
  ) |>
  full_join(pc_long, by = c("immet_id", "exam_order")) |>
  group_by(clin_name, pc_name) |>
  nest() |>
  ungroup()


clin_long <- bind_rows(clin_long_num, clin_long_fac)


res_clin_act_time_01 <- clin_long |> 
  mutate(mod = map(data, ~possLMER(data = .x,
                               pc_value ~ exam_order*clin_value + bmi + age_m0 + 
                                 (1 | immet_id))),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE)))

res_clin_act_time_01_tab <- res_clin_act_time_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "clin_value")) |> 
  mutate(term = str_replace(term, "clin_value", clin_name),
         pc_name = str_remove(pc_name, "_pc1")) |> 
  select(-data, -mod, -x, -std.error, -statistic, -df, -clin_name, -effect)

 
res_clin_act_sel_01 <- res_clin_act_time_01_tab |> 
  filter(p.value < 0.06 & pc_name == "activ") |> 
  pull(term) |> 
  str_remove("exam_orderm6:") |> 
  unique() |> 
  str_replace("response_m0_m6Minimal", "response_m0_m6") |> 
  str_replace("mtx1", "mtx") |> 
  str_replace("cpa1", "cpa") |> 
  str_replace("arthritis1", "arthritis") |> 
  str_replace("raynaud1", "raynaud") |> 
  str_replace("oesophageal_motility_disorder1", "oesophageal_motility_disorder")

res_clin_time_sel_01 <- res_clin_act_time_01_tab |> 
  filter(p.value < 0.06 & pc_name == "time") |> 
  pull(term) |> 
  str_remove("exam_orderm6:") |> 
  unique() |> 
  str_replace("disease_subtype2", "disease_subtype") |> 
  str_replace("mtx1", "mtx") |> 
  str_replace("cpa1", "cpa") |> 
  str_replace("arthritis1", "arthritis") |> 
  str_replace("raynaud1", "raynaud") |> 
  str_replace("ild1", "ild") |> 
  str_replace("statins1", "statins")


# **********************************************************************
## Confirmatory analyses ----
# **********************************************************************
### LMER ----
# **********************************************************************

clin_long_num_02 <- d15_scaled |> 
  select(immet_id, exam_order, all_of(x_activ), all_of(x_time), all_of(base_covars),
         intersect(res_clin_act_sel_01, clin_num)) |> 
  pivot_longer(cols = intersect(res_clin_act_sel_01, clin_num), 
               names_to = "clin_name", values_to = "clin_value") |> 
  pivot_longer(cols = all_of(c(x_activ, x_time)), 
               names_to = "var_name", values_to = "var_value") |> 
  mutate(clin_type = "numeric",
         clin_value = as.numeric(clin_value)) |> 
  group_by(clin_name, var_name) |> 
  nest() |> 
  ungroup()


clin_long_fac_02 <- d15_scaled |> 
  select(immet_id, exam_order, all_of(x_activ), all_of(x_time), all_of(base_covars),
         intersect(res_clin_act_sel_01, clin_fac)) |> 
  pivot_longer(cols = intersect(res_clin_act_sel_01, clin_fac), 
               names_to = "clin_name", values_to = "clin_value") |> 
  pivot_longer(cols = all_of(c(x_activ, x_time)), 
               names_to = "var_name", values_to = "var_value") |> 
  mutate(clin_type = "factor",
         clin_value = as.factor(clin_value)) |> 
  group_by(clin_name, var_name) |> 
  nest() |> 
  ungroup()

clin_long_02 <- bind_rows(clin_long_num_02, clin_long_fac_02)


res_clin_act_time_02 <- clin_long_02 |> 
  mutate(mod = map(data, ~possLMER(data = .x,
                                   var_value ~ exam_order*clin_value + bmi + age_m0 + 
                                     (1 | immet_id))),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE)))


res_clin_act_time_02_tab <- res_clin_act_time_02 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "clin_value")) |> 
  mutate(term = str_replace(term, "clin_value", clin_name),
         term = str_replace(term, "exam_orderm6", "exam")) |> 
  select(-data, -mod, -std.error, -statistic, -df,  -effect) |> 
  mutate(ggir_type = case_when(
    var_name %in% x_activ ~ "active",
    var_name %in% x_time ~ "time"
  ))

res_clin_act_time_02_tab_sig <- res_clin_act_time_02_tab |>
  dplyr::group_by(var_name, clin_name) |>
  dplyr::filter(
    any(
      (p.value < 0.06) & (sign(conf.low) == sign(conf.high)),
      na.rm = TRUE)) |>
  dplyr::ungroup()

export(res_clin_act_time_02_tab_sig,
       "output/tables/final/260126_assoc_clin_activ_time_signif_01.xlsx")

plot_emm_clin_exam <- function(model, clin_name, var_name) {
  sjPlot::plot_model(
    model = model,
    type  = "emm",
    term  = c("clin_value", "exam_order")
  ) +
    sjPlot::theme_sjplot2() +
    ggplot2::labs(
      title    = "Association between variables",
      x = clin_name,
      y = var_name
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
      axis.title    = ggplot2::element_text(size = 12, face = "bold"),
      axis.text     = ggplot2::element_text(size = 11, face = "bold")
    )
}


res_clin_act_time_02_fig <- res_clin_act_time_02 |>
  dplyr::semi_join(
    res_clin_act_time_02_tab_sig |>
      dplyr::select(clin_name, var_name) |>
      dplyr::distinct(),
    by = c("clin_name", "var_name")) |> 
  mutate(fig = pmap(list(mod, clin_name, var_name),
                    ~plot_emm_clin_exam(model = ..1, clin_name = ..2, var_name = ..3)))


pdf(file.path("output/figures", 
              paste0(format(Sys.Date(), "%y%m%d"),
              "_", 
              "signif_assoc_clin_time_activ_01.pdf")), 
    width = 10, height = 7)
walk(res_clin_act_time_02_fig$fig, print)
dev.off()



# **********************************************************************
### sPLS ----
# **********************************************************************

clin_sel <- intersect(res_clin_act_time_02_tab_sig$clin_name |> unique(), 
                      names(d15_scaled))
y_activ  <- intersect(x_activ,  names(d15_scaled))
y_time   <- intersect(x_time,   names(d15_scaled))

df <- d15_scaled


# ******************************
# 1) Helper: prepare X/Y matrices for a given timepoint
# ******************************
prep_xy <- function(df, exam_level, x_vars, y_vars, id_var = "immet_id") {
  
  dat <- df |>
    filter(exam_order == exam_level) |>
    select(all_of(id_var), all_of(x_vars), all_of(y_vars)) |>
    distinct()
  
  # drop columns with 0 variance (mixOmics will complain)
  nzv <- function(v) sd(v, na.rm = TRUE) > 0
  x_keep <- x_vars |> keep(~ is.numeric(dat[[.x]]) && nzv(dat[[.x]]))
  y_keep <- y_vars |> keep(~ is.numeric(dat[[.x]]) && nzv(dat[[.x]]))
  
  # NOTE: sPLS works on numeric matrices; factors need to be encoded beforehand if you want them in X
  # Here we keep only numeric clinical predictors. If you want binary factors, convert 0/1 to numeric in advance.
  X <- dat |> select(all_of(x_keep)) |> as.matrix()
  Y <- dat |> select(all_of(y_keep)) |> as.matrix()
  
  # simple median imputation (recommended)
  for (j in seq_len(ncol(X))) {
    if (anyNA(X[, j])) X[is.na(X[, j]), j] <- median(X[, j], na.rm = TRUE)
  }
  for (j in seq_len(ncol(Y))) {
    if (anyNA(Y[, j])) Y[is.na(Y[, j]), j] <- median(Y[, j], na.rm = TRUE)
  }
  
  list(X = X, Y = Y, x_keep = x_keep, y_keep = y_keep)
}

# ******************************
# 2) Helper: tune sPLS
# ******************************
tune_spls_xy <- function(X, Y, ncomp = 2,
                         test_keepX = NULL,
                         test_keepY = NULL,
                         folds = 5, nrepeat = 30, seed = 1) {
  
  set.seed(seed)
  
  if (is.null(test_keepX)) {
    # if your X is already a shortlist, you can also set keepX = ncol(X) and skip tuning on X
    test_keepX <- unique(pmin(c(3, 5, 8, 10, 15), ncol(X)))
  } else {
    test_keepX <- unique(pmin(test_keepX, ncol(X)))
  }
  
  if (is.null(test_keepY)) {
    test_keepY <- unique(pmin(c(3, 5, 8, 10, 15), ncol(Y)))
  } else {
    test_keepY <- unique(pmin(test_keepY, ncol(Y)))
  }
  
  tune <- tune.spls(
    X, Y,
    ncomp = ncomp,
    test.keepX = test_keepX,
    test.keepY = test_keepY,
    validation = "Mfold",
    folds = folds,
    nrepeat = nrepeat,
    progressBar = TRUE,
    scale = TRUE   # important unless you are 100% sure everything is already comparable
  )
  
  # choice: one keepX/keepY per component
  # tune$choice.keepX and $choice.keepY are vectors of length ncomp
  list(
    tune = tune,
    keepX = tune$choice.keepX,
    keepY = tune$choice.keepY
  )
}

# ******************************
# 3) Fit + summarise
# ******************************
fit_spls_xy <- function(X, Y, keepX, keepY, ncomp = 2) {
  fit <- spls(X, Y, ncomp = ncomp, keepX = keepX, keepY = keepY, scale = TRUE)
  
  # selected variables
  selX <- lapply(seq_len(ncomp), function(h) selectVar(fit, comp = h)$name)
  selY <- lapply(seq_len(ncomp), function(h) selectVar(fit, comp = h, block = "Y")$name) # sometimes not needed
  names(selX) <- paste0("comp", seq_len(ncomp))
  names(selY) <- paste0("comp", seq_len(ncomp))
  
  list(fit = fit, selX = selX, selY = selY)
}

# ******************************
# 4) Run: AKTIVITA and CAS, separately for M0 and M6
# ******************************

run_one <- function(exam_level, y_vars, label) {
  
  xy <- prep_xy(df, exam_level = exam_level, x_vars = clin_sel, y_vars = y_vars)
  
  # if your clin_sel is already small, you can force keepX = ncol(X) (no selection on X)
  # and only tune Y:
  # keepX_fixed <- rep(ncol(xy$X), 2)
  # keepY_tune <- tune_spls_xy(xy$X, xy$Y, ncomp = 2, test_keepX = keepX_fixed)
  
  tuned <- tune_spls_xy(xy$X, xy$Y, ncomp = 2, folds = 5, nrepeat = 30, seed = 1)
  
  fitted <- fit_spls_xy(xy$X, xy$Y, keepX = tuned$keepX, keepY = tuned$keepY, ncomp = 2)
  
  list(
    exam = exam_level,
    block = label,
    X_vars_used = xy$x_keep,
    Y_vars_used = xy$y_keep,
    tune = tuned$tune,
    keepX = tuned$keepX,
    keepY = tuned$keepY,
    fit = fitted$fit,
    selected_X = fitted$selX,
    selected_Y = fitted$selY
  )
}

# --- Activity ---
spls_activ_m0 <- run_one("m0", y_activ, "AKTIVITA")
spls_activ_m6 <- run_one("m6", y_activ, "AKTIVITA")

# --- Time ---
spls_time_m0  <- run_one("m0", y_time, "CAS")
spls_time_m6  <- run_one("m6", y_time, "CAS")

# ******************************
## Quick outputs / plots ----
# ******************************

# Example: correlation of variates (X vs Y) for M0 activity
plotIndiv(spls_activ_m0$fit, comp = c(1,2), group = NULL, legend = TRUE,
          title = "sPLS (M0) - AKTIVITA - X variates")
plotVar(spls_activ_m0$fit, comp = c(1,2), var.names = TRUE,
        title = "sPLS (M0) - AKTIVITA - selected variables")

plotVar(spls_time_m0$fit, comp = c(1,2), var.names = TRUE,
        title = "sPLS (M0) - TIME - selected variables")


# Network plot (often very useful)
network(spls_activ_m0$fit, comp = 1, cutoff = 0.3, title = "sPLS network (M0) - AKTIVITA - comp1")
network(spls_time_m0$fit,  comp = 1, cutoff = 0.3, title = "sPLS network (M0) - CAS - comp1")

# Selected variables (lists)
spls_activ_m0$selected_X
spls_activ_m0$selected_Y

# **********************************************************************
# 3. Accelerometry - subtype ----
# **********************************************************************




# **********************************************************************
## LMER ----
# **********************************************************************
### subtype ----
# **********************************************************************
res_subtype_act_time_01 <- d15_scores |> 
  select(immet_id, exam_order, all_of(pc_cols), disease_subtype, all_of(base_covars)) |>
  pivot_longer(cols = all_of(pc_cols), names_to = "pc_name", values_to = "pc_value") |> 
  group_by(pc_name) |> 
  nest() |> 
  mutate(mod = map(data, ~possLMER(data = .x,
                                   pc_value ~ exam_order + disease_subtype + bmi + age_m0 + 
                                     (1 | immet_id))),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE)))

res_subtype_act_time_01_tab <- res_subtype_act_time_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "disease_subtype")) |> 
  mutate(pc_name = str_remove(pc_name, "_pc1")) |> 
  select(-data, -mod, -std.error, -statistic, -df, -effect)

res_subtype_act_sel_01 <- res_subtype_act_time_01_tab |> 
  filter(p.value < 0.06 & pc_name == "activ") 

res_subtype_time_sel_01 <- res_subtype_act_time_01_tab |> 
  filter(p.value < 0.06 & pc_name == "time")

# **********************************************************************
### response ----
# **********************************************************************
res_response_act_time_01 <- d15_scores |> 
  select(immet_id, exam_order, all_of(pc_cols), response_m0_m6, all_of(base_covars)) |>
  pivot_longer(cols = all_of(pc_cols), names_to = "pc_name", values_to = "pc_value") |> 
  group_by(pc_name) |> 
  nest() |> 
  mutate(mod = map(data, ~possLMER(data = .x,
                                   pc_value ~ exam_order + response_m0_m6 + bmi + age_m0 + 
                                     (1 | immet_id))),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE)))

res_response_act_time_01_tab <- res_response_act_time_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "response")) |> 
  mutate(pc_name = str_remove(pc_name, "_pc1")) |> 
  select(-data, -mod, -std.error, -statistic, -df, -effect)

res_response_act_sel_01 <- res_response_act_time_01_tab |> 
  filter(p.value < 0.06 & pc_name == "activ") 

res_response_time_sel_01 <- res_response_act_time_01_tab |> 
  filter(p.value < 0.06 & pc_name == "time")

# **********************************************************************
## Confirmatory ----
# **********************************************************************
safe_emm   <- possibly(emmeans, otherwise = NULL)
safe_contr <- possibly(contrast, otherwise = NULL)
safe_tidy  <- possibly(broom::tidy, otherwise = tibble())

res_subtype_act_02 <- d15_scores |> 
  select(immet_id, exam_order, all_of(x_activ), all_of(base_covars), disease_subtype) |> 
  pivot_longer(cols = all_of(x_activ), 
               names_to = "var_name", values_to = "var_value") |> 
  group_by(var_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, ~possLMER(data = .x,
                                   var_value ~ exam_order + disease_subtype + bmi + age_m0 + 
                                     (1 | immet_id))),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE)))

res_subtype_act_02_contrast <- res_subtype_act_02 |> unnest(tidier) |> 
  filter(str_detect(term, "disease_subtype") &
           (p.value < 0.06) & (sign(conf.low) == sign(conf.high)) ) |> 
  mutate(emm_subtype_by_exam = map(mod,
    ~ safe_emm(.x, specs = ~ disease_subtype | exam_order)),
    contr_subtype_by_exam = map(emm_subtype_by_exam,
      ~ safe_contr(.x, method = "pairwise", adjust = "tukey")),
    contr_tidy = map(contr_subtype_by_exam,
      ~ if (is.null(.x)) tibble() else tidy(.x, conf.int = TRUE))) |> 
  select(var_name, contr_tidy) |> 
  unnest(contr_tidy)

res_response_act_02 <- d15_scores |> 
  select(immet_id, exam_order, all_of(x_activ), all_of(base_covars), response_m0_m6) |> 
  pivot_longer(cols = all_of(x_activ), 
               names_to = "var_name", values_to = "var_value") |> 
  group_by(var_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, ~possLMER(data = .x,
                                   var_value ~ exam_order + response_m0_m6 + bmi + age_m0 + 
                                     (1 | immet_id))),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE)))

res_response_act_02_contrast <- res_response_act_02 |> unnest(tidier) |> 
  filter(str_detect(term, "response") &
           (p.value < 0.06) & (sign(conf.low) == sign(conf.high)) ) |> 
  mutate(emm_subtype_by_exam = map(mod,
                                   ~ safe_emm(.x, specs = ~ disease_subtype | exam_order)),
         contr_subtype_by_exam = map(emm_subtype_by_exam,
                                     ~ safe_contr(.x, method = "pairwise", adjust = "tukey")),
         contr_tidy = map(contr_subtype_by_exam,
                          ~ if (is.null(.x)) tibble() else tidy(.x, conf.int = TRUE))) |> 
  select(var_name, contr_tidy) |> 
  unnest(contr_tidy)

# 
res_response_act_02b <- d15_scores |> 
  select(immet_id, exam_order, all_of(x_activ), all_of(base_covars), response_m0_m6) |> 
  mutate(response_m0_m6 = factor(response_m0_m6, ordered = FALSE)) |> 
  pivot_longer(cols = all_of(x_activ), 
               names_to = "var_name", values_to = "var_value") |> 
  group_by(var_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, ~possLMER(data = .x,
                                   var_value ~ exam_order + response_m0_m6 + bmi + age_m0 + 
                                     (1 | immet_id))),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE)))

res_response_act_02b_contrast <- res_response_act_02b |> unnest(tidier) |> 
  filter(str_detect(term, "response") &
           (p.value < 0.06) & (sign(conf.low) == sign(conf.high)) ) |> 
  mutate(emm_subtype_by_exam = map(mod,
                                   ~ safe_emm(.x, specs = ~ disease_subtype | exam_order)),
         contr_subtype_by_exam = map(emm_subtype_by_exam,
                                     ~ safe_contr(.x, method = "pairwise", adjust = "tukey")),
         contr_tidy = map(contr_subtype_by_exam,
                          ~ if (is.null(.x)) tibble() else tidy(.x, conf.int = TRUE))) |> 
  select(var_name, contr_tidy) |> 
  unnest(contr_tidy)

# **********************************************************************
### response ----
# **********************************************************************
res_response_time_01 <- d15_scores |>
  select(immet_id, exam_order, all_of(x_time), all_of(base_covars), response_m0_m6) |>
  mutate(
    exam_order     = factor(exam_order),
    response_m0_m6 = factor(response_m0_m6, ordered = TRUE)  # necháš jako ordered, ale níže testuješ přes emmeans
  ) |>
  pivot_longer(
    cols = all_of(x_time),
    names_to = "var_name",
    values_to = "var_value"
  ) |>
  group_by(var_name) |>
  nest() |>
  ungroup() |>
  mutate(
    mod = map(data, ~ possLMER(
      data = .x,
      var_value ~ exam_order + response_m0_m6 + bmi + age_m0 + (1 | immet_id)
    )),
    tidier = map(mod, ~ tidy(.x, effects = "fixed", conf.int = TRUE))
  )

# 1) Které time parametry mají významný efekt response (v modelu)?
sig_response_time_01 <- res_response_time_01 |>
  transmute(var_name, tidier) |>
  unnest(tidier) |>
  filter(str_detect(term, "^response_m0_m6")) |>
  group_by(var_name) |>
  summarise(
    min_p_response = min(p.value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(min_p_response < 0.05)

sig_response_time_01
# -> Tohle je seznam time proměnných, kde vyšel response efekt (alespoň jeden kontrast)

# 2) Pro významné proměnné: emmeans(response | exam_order) + pairwise kontrasty
res_response_time_01_contrast <- res_response_time_01 |>
  semi_join(sig_response_time_01, by = "var_name") |>
  mutate(
    emm_resp_by_exam = map(mod, ~ safe_emm(.x, specs = ~ response_m0_m6 | exam_order)),
    contr_resp_by_exam = map(emm_resp_by_exam, ~ safe_contr(.x, method = "pairwise", adjust = "tukey")),
    contr_tidy = map(contr_resp_by_exam, ~ if (is.null(.x)) tibble() else tidy(.x, conf.int = TRUE))
  ) |>
  select(var_name, contr_tidy) |>
  unnest(contr_tidy) |>
  # volitelně: jen "jasně významné" kontrasty (p<0.05 a CI mimo 0)
  filter(adj.p.value < 0.06, sign(conf.low) == sign(conf.high)) |>
  arrange(var_name, exam_order, adj.p.value)

res_response_time_01_contrast
# -> konkrétní významné rozdíly mezi response skupinami zvlášť pro M0 a M6


d14_join |> 
  select(immet_id, exam_order, response_m0_m6, disease_subtype, bmi, age_m0) |> 
  mutate(response_m0_m6 = response_m0_m6 |> 
           tolower() |> 
           factor(levels = c("no improvement", "minimal", 
                                "moderate", "major")),
         mod = ordinal::clmm(response_m0_m6 ~ exam_order * disease_subtype + 
                               bmi + age_m0 + (1|immet_id)))
# **********************************************************************
### response - ordered ----
# **********************************************************************

possCLMM <- function(data, formula, ...) {
  tryCatch(
    ordinal::clmm(formula, data = data, Hess = TRUE, ...),
    error = function(e) NULL
  )
}

library(dplyr)
library(tidyr)
library(purrr)
library(ordinal)
library(emmeans)
library(broom)

res_ord_01 <- d15_scores |>
  select(immet_id, exam_order, response_m0_m6, disease_subtype, bmi, age_m0) |>
  mutate(
    response_m0_m6 = response_m0_m6 |>
      tolower() |>
      factor(levels = c("no improvement", "minimal", "moderate", "major")) |>
      ordered(),
    exam_order = factor(exam_order),
    disease_subtype = factor(disease_subtype)
  ) |>
  mutate(.grp = "all") |>
  group_by(.grp) |>
  nest() |>
  ungroup() |>
  mutate(
    mod = map(data, ~ ordinal::clmm(
      data = .x,
      response_m0_m6 ~ exam_order + disease_subtype + bmi + age_m0 + (1 | immet_id),
      link = "logit"
    )),
    mod_m0 = map(data, ~ ordinal::clmm(
      data = .x,
      response_m0_m6 ~ exam_order + bmi + age_m0 + (1 | immet_id),
      link = "logit"
    )),
    emms_lat = map(mod, ~ emmeans(.x, ~ disease_subtype | exam_order, mode = "latent")),
    contrasts_lat = map(emms_lat, ~ pairs(.x, adjust = "tukey")),
    contrasts_lat_tbl = map(contrasts_lat, ~ broom::tidy(.x, conf.int = TRUE)),
    emms_prob = map(mod, ~ emmeans(.x, ~ disease_subtype | exam_order, mode = "prob")),
    emms_prob_tbl = map(emms_prob, as.data.frame)
  ) |>
  select(-.grp)

res_ord_01_tab_summary <- res_ord_01$mod[[1]] |> tidy(conf.int = TRUE)
res_ord_01_tab_contrast <- res_ord_01$contrasts_lat_tbl[[1]]

res_ord_01_anova <- anova(
  res_ord_01$mod_m0[[1]],
  res_ord_01$mod[[1]]
) |> as.data.frame()


export(list(model_sum = res_ord_01_tab_summary,
            model_cotrast = res_ord_01_tab_contrast,
            model_comp_m0 = res_ord_01_anova),
       "output/tables/final/260129_resp_diseas_type_output_01.xlsx")

# **********************************************************************
# 4. figures ----
# **********************************************************************
## line graf ----
# **********************************************************************


fig_sel_01 <- c("enmo_full_recording_mean", "inactivity", "light", 
                "moderate", "t400_time_min")

fig_sel_02 <- c("ig_gradient_enmo", "ig_intercept_enmo", "m5_enmo")

fig_sel_03 <- c("t100_b1m80_time_min", "t100_b5m80_time_min", "t100_b10m80_time_min")

fig_sel_04 <- c("p9931_enmo", "p9861_enmo", "p9792_enmo", "p9583_enmo","p9167_enmo")



# ************************************************************
# FIG: line plots M0 vs M6 by response, 4 vars, patchwork 
# ************************************************************

fig_line_m0m6_respcol_build_df <- function(
    data,
    vars_sel,
    id_col = "immet_id",
    response_col = "response_m0_m6",
    exam_col = "exam_order",
    exam_levels = c("m0", "m6"),
    out_var_name = "var_name",
    out_var_value = "var_value"
) {
  data |>
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
}

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
    )
) {
  ggplot(df_var, aes(
    x = .data[[exam_col]],
    y = .data[[value_col]],
    col = .data[[response_col]],
    shape = .data[[response_col]],
    group = .data[[id_col]]
  )) +
    geom_line(linewidth = 0.9, alpha = 0.5) +
    geom_point(
      position = position_jitter(width = 0.08, height = 0),
      size = 2, alpha = 0.5
    ) +
    scale_color_manual(values = color_values) +
    labs(title = var_label) +
    facet_grid(rows = vars(.data[[response_col]]), scales = "fixed") +
    sjPlot::theme_sjplot2() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text  = element_text(size = 10, face = "bold"),
      strip.text.y = element_text(size = 11, face = "bold"),
      legend.position = "none"
    )
}

fig_line_m0m6_respcol_patchwork_4cols <- function(
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
    )
) {
  fig_line_m0m6_respcol_df <- fig_line_m0m6_respcol_build_df(
    data = data,
    vars_sel = vars_sel,
    id_col = id_col,
    response_col = response_col,
    exam_col = exam_col,
    exam_levels = exam_levels
  )
  
  # 1 plot = 1 sloupec (1 var_name), uvnitř 4 řádky (response)
  fig_line_m0m6_respcol_plots <- fig_line_m0m6_respcol_df |>
    group_by(var_name) |>
    group_split() |>
    imap(~ fig_line_m0m6_respcol_make_var_column(
      df_var = .x,
      var_label = as.character(unique(.x$var_name)),
      id_col = id_col,
      response_col = response_col,
      exam_col = exam_col,
      value_col = "var_value",
      color_values = color_values
    ))
  
  # složit do 1 řady (4 sloupce); každý sloupec má vlastní y (protože je to samostatný ggplot)
  fig_line_m0m6_respcol_patch <- patchwork::wrap_plots(fig_line_m0m6_respcol_plots, nrow = 1)
  
  list(
    fig_line_m0m6_respcol_df    = fig_line_m0m6_respcol_df,
    fig_line_m0m6_respcol_plots = fig_line_m0m6_respcol_plots,
    fig_line_m0m6_respcol_patch = fig_line_m0m6_respcol_patch
  )
}

# ---- použít:
fig_line_m0m6_respcol_out <- fig_line_m0m6_respcol_patchwork_4cols(
  data = d14_join,
  vars_sel = fig_sel_01
)

fig_line_01 <- fig_line_m0m6_respcol_out$fig_line_m0m6_respcol_patch
tiff(
  filename = file.path("output/figures/final",
                       paste0(format(Sys.Date(), "%y%m%d"),
                       "_", 
                       "slope_line_chart_01.tiff")),
  compression = "lzw",
  res = "150",
  width = 1400,
  height = 1000
)
print(fig_line_01)
dev.off()

# **********************************************************************
# 5. myokines ----
# **********************************************************************
d17_myokines <- import(paste0(folder_path, "Treatment Response_sérum,expr_pre STAT .xlsx")) |> 
  select(!contains("pg")) |> 
  rename(immet_id = `projekt ID`) |> 
  mutate(immet_id = str_remove_all(immet_id, "_"))


# **********************************************************************
# Data: jen numerické proměnné
# **********************************************************************
# **********************************************************************
# PCA myokines (tidymodels-native)
# **********************************************************************

# 1) Recipe + prep
pca_myokines_rec <- recipe(~ ., data = d17_myokines) |>
  update_role(immet_id, new_role = "id") |>
  step_zv(all_numeric_predictors(), id = "zv") |>
  step_normalize(all_numeric_predictors(), id = "norm") |>
  step_pca(all_numeric_predictors(), threshold = 0.90, id = "pca")

pca_myokines_prep <- prep(pca_myokines_rec)

# 2) Skupiny
pca_myokines_grp <- d14_join |>
  select(immet_id, response_m0_m6) |>
  distinct()

# 3) Scores + skupiny
pca_myokines_scores <- juice(pca_myokines_prep) |>
  left_join(pca_myokines_grp, by = "immet_id")

# 4) Explained variance (%)
pca_myokines_var <- tidy(pca_myokines_prep, id = "pca", type = "variance") |>
  filter(terms == "variance") |>
  mutate(
    percent  = 100 * value / sum(value),
    cum_perc = cumsum(percent)
  )

pc1_perc <- pca_myokines_var |> filter(component == 1) |> pull(percent)
pc2_perc <- pca_myokines_var |> filter(component == 2) |> pull(percent)

# 5) Scree plot
fig_pca_myokines_scree <-
  ggplot(pca_myokines_var, aes(component, percent)) +
  geom_col() +
  geom_line(aes(group = 1)) +
  geom_point() +
  labs(y = "% explained variance", x = "PC") +
  theme_minimal(base_size = 12)

# 6) PC1 vs PC2 plot (filled ellipses)
fig_pca_myokines_pc1_pc2 <-
  ggplot(
    pca_myokines_scores,
    aes(PC1, PC2, color = response_m0_m6, fill = response_m0_m6)
  ) +
  geom_point(size = 2.6, alpha = 0.85) +
  stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.15, linewidth = 1) +
  labs(
    title = "PCA - myokines: PC1 vs PC2",
    x = paste0("PC1 (", round(pc1_perc, 1), "%)"),
    y = paste0("PC2 (", round(pc2_perc, 1), "%)"),
    color = "Response",
    fill  = "Response"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 20),
    axis.title   = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold")
  ) +
  paletteer::scale_color_paletteer_d("MetBrewer::Egypt") +
  paletteer::scale_fill_paletteer_d("MetBrewer::Egypt")

# 7) Export PC1–PC2
tiff(
  filename = file.path(
    "output/figures",
    paste0(format(Sys.Date(), "%y%m%d"), "_pca_myokines_pc1_pc2.tiff")
  ),
  compression = "lzw",
  res = 150,
  width = 1600,
  height = 1200
)
print(fig_pca_myokines_pc1_pc2)
dev.off()

# (volitelně) zobrazit v RStudio
fig_pca_myokines_scree
fig_pca_myokines_pc1_pc2

# **********************************************************************
# 6. lmer - association analyses - PtVAS covariate ----
# **********************************************************************
## clinics ----
# **********************************************************************
### time ----
# **********************************************************************
res_ass_time_clincs_01 <- d14_join |> 
  select(immet_id, exam_order, all_of(sel_time_01), all_of(sel_clinics_01), physician_vas, age_m0) |> 
  pivot_longer(cols = all_of(sel_time_01),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_clinics_01),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order*scale(indep_value) + scale(physician_vas) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_time_clincs_01_tab <- res_ass_time_clincs_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "scale\\(indep_value\\)", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_time_clincs_01_fig <- res_ass_time_clincs_01 |> 
  inner_join(res_ass_time_clincs_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Clinical and Time features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 11, face = "bold")
                      )
  ))


pdf(file.path("output/figures", 
              paste0(format(Sys.Date(), "%y%m%d"),
                     "_", 
                     "assoc_time_clin_01.pdf")), 
    width = 10, height = 7)
walk(res_ass_time_clincs_01_fig$fig, print)
dev.off()

export(res_ass_time_clincs_01_tab, "output/tables/260130_ass_time_clincs_01.xlsx")

# **********************************************************************
### activity ----
# **********************************************************************
res_ass_activ_clincs_01 <- d14_join |> 
  select(immet_id, exam_order, all_of(sel_activ_01), all_of(sel_clinics_01), physician_vas, age_m0) |> 
  pivot_longer(cols = all_of(sel_activ_01),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_clinics_01),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order*scale(indep_value) + scale(physician_vas) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_activ_clincs_01_tab <- res_ass_activ_clincs_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "scale\\(indep_value\\)", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_activ_clincs_01_fig <- res_ass_activ_clincs_01 |> 
  inner_join(res_ass_activ_clincs_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Clinical and Activity features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 11, face = "bold")
                      )
  ))


pdf(file.path("output/figures", 
              paste0(format(Sys.Date(), "%y%m%d"),
                     "_", 
                     "assoc_activ_clin_01.pdf")), 
    width = 10, height = 7)
walk(res_ass_activ_clincs_01_fig$fig, print)
dev.off()

export(res_ass_activ_clincs_01_tab, "output/tables/260130_ass_activ_clincs_01.xlsx")

# **********************************************************************
##  myokines ----
# **********************************************************************
### time ----
# **********************************************************************
res_ass_time_myokin_01 <- d14_join |> 
  inner_join(d17_myokines, by = "immet_id") |> 
  select(immet_id, exam_order, all_of(sel_time_01), all_of(sel_myokines_01),  physician_vas, age_m0) |> 
  pivot_longer(cols = all_of(sel_time_01),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_myokines_01),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order + indep_value + scale(physician_vas) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_time_myokin_01_tab <- res_ass_time_myokin_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "indep_value", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_time_myokin_01_fig <- res_ass_time_myokin_01 |> 
  inner_join(res_ass_time_myokin_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Myokines and Time features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 11, face = "bold")
                      )
  ))

# # No significant 
# pdf(file.path("output/figures", 
#               paste0(format(Sys.Date(), "%y%m%d"),
#                      "_", 
#                      "assoc_time_myokin_01.pdf")), 
#     width = 10, height = 7)
# walk(res_ass_time_myokin_01_fig$fig, print)
# dev.off()

export(res_ass_time_myokin_01_tab, "output/tables/260130_ass_time_myokin_01.xlsx")

# **********************************************************************
### activity ----
# **********************************************************************
res_ass_activ_myokin_01 <- d14_join |> 
  inner_join(d17_myokines, by = "immet_id") |> 
  select(immet_id, exam_order, all_of(sel_activ_01), all_of(sel_myokines_01), physician_vas, age_m0) |> 
  pivot_longer(cols = all_of(sel_activ_01),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_myokines_01),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order + indep_value + scale(physician_vas) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_activ_myokin_01_tab <- res_ass_activ_myokin_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "indep_value", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_activ_myokin_01_fig <- res_ass_activ_myokin_01 |> 
  inner_join(res_ass_activ_myokin_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Myokines and Activity features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 11, face = "bold")
                      )
  ))

# # No significant 
# pdf(file.path("output/figures", 
#               paste0(format(Sys.Date(), "%y%m%d"),
#                      "_", 
#                      "assoc_activ_myokin_01.pdf")), 
#     width = 10, height = 7)
# walk(res_ass_activ_myokin_01_fig$fig, print)
# dev.off()

export(res_ass_activ_myokin_01_tab, "output/tables/260130_ass_activ_myokin_01.xlsx")

# **********************************************************************
# 7. lmer - association analyses - bcm covariate ----
# **********************************************************************
# **********************************************************************
## clinics ----
# **********************************************************************
### time ----
# **********************************************************************
res_ass_time_clincs_bcmADJ_01 <- d14_join |> 
  mutate(creatinine = dplyr::if_else(creatinine > 120, NA_real_, creatinine)) |> 
  select(immet_id, exam_order, all_of(sel_time_02), all_of(sel_clinics_02), bcm, age_m0) |> 
  pivot_longer(cols = all_of(sel_time_02),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_clinics_02),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order*scale(indep_value) + scale(bcm) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_time_clincs_bcmADJ_01_tab <- res_ass_time_clincs_bcmADJ_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "scale\\(indep_value\\)", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_time_clincs_bcmADJ_01_fig <- res_ass_time_clincs_bcmADJ_01 |> 
  inner_join(res_ass_time_clincs_bcmADJ_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 show.data = TRUE,
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Clinical and Time features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 18, face = "bold")
                      )
  ))


pdf(file.path("output/figures", 
              paste0(format(Sys.Date(), "%y%m%d"),
                     "_", 
                     "assoc_time_clin_bcmADJ_01.pdf")), 
    width = 10, height = 7)
walk(res_ass_time_clincs_bcmADJ_01_fig$fig, print)
dev.off()

export(res_ass_time_clincs_bcmADJ_01_tab, "output/tables/260131_ass_time_clincs_bcmADJ_01.xlsx")

# **********************************************************************
### activity ----
# **********************************************************************
res_ass_activ_clincs_bcmADJ_01 <- d14_join |> 
  mutate(creatinine = dplyr::if_else(creatinine > 120, NA_real_, creatinine)) |> 
  select(immet_id, exam_order, all_of(sel_activ_02), all_of(sel_clinics_02), bcm, age_m0) |> 
  pivot_longer(cols = all_of(sel_activ_02),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_clinics_02),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order*scale(indep_value) + scale(bcm) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_activ_clincs_bcmADJ_01_tab <- res_ass_activ_clincs_bcmADJ_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "scale\\(indep_value\\)", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_activ_clincs_bcmADJ_01_fig <- res_ass_activ_clincs_bcmADJ_01 |> 
  inner_join(res_ass_activ_clincs_bcmADJ_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 show.data = TRUE,
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Clinical and Activity features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 18, face = "bold")
                      )
  ))


pdf(file.path("output/figures", 
              paste0(format(Sys.Date(), "%y%m%d"),
                     "_", 
                     "assoc_activ_clin_bcmADJ_01.pdf")), 
    width = 10, height = 7)
walk(res_ass_activ_clincs_bcmADJ_01_fig$fig, print)
dev.off()

export(res_ass_activ_clincs_bcmADJ_01_tab, "output/tables/260131_ass_activ_clincs_bcmADJ_01.xlsx")

# **********************************************************************
## myokines ----
# **********************************************************************
### time ----
# **********************************************************************
res_ass_time_myokin_bcmADJ_01 <- d14_join |> 
  mutate(creatinine = dplyr::if_else(creatinine > 120, NA_real_, creatinine)) |> 
  inner_join(d17_myokines, by = "immet_id") |> 
  select(immet_id, exam_order, all_of(sel_time_02), all_of(sel_myokines_01),  bcm, age_m0) |> 
  pivot_longer(cols = all_of(sel_time_02),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_myokines_01),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order + indep_value + scale(bcm) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_time_myokin_bcmADJ_01_tab <- res_ass_time_myokin_bcmADJ_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "indep_value", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_time_myokin_bcmADJ_01_fig <- res_ass_time_myokin_bcmADJ_01 |> 
  inner_join(res_ass_time_myokin_bcmADJ_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 show.data = TRUE,
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Myokines and Time features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 11, face = "bold")
                      )
  ))


# pdf(file.path("output/figures", 
#               paste0(format(Sys.Date(), "%y%m%d"),
#                      "_", 
#                      "assoc_time_myokin_bcmADJ_01.pdf")), 
#     width = 10, height = 7)
# walk(res_ass_time_myokin_bcmADJ_01_fig$fig, print)
# dev.off()

export(res_ass_time_myokin_bcmADJ_01_tab, "output/tables/260131_ass_time_myokin_bcmADJ_01.xlsx")

# **********************************************************************
### activity ----
# **********************************************************************
res_ass_activ_myokin_bcmADJ_01 <- d14_join |> 
  mutate(creatinine = dplyr::if_else(creatinine > 120, NA_real_, creatinine)) |> 
  inner_join(d17_myokines, by = "immet_id") |> 
  select(immet_id, exam_order, all_of(sel_activ_01), all_of(sel_myokines_01), bcm, age_m0) |> 
  pivot_longer(cols = all_of(sel_activ_01),
               names_to = "dep_name",
               values_to = "dep_value") |> 
  pivot_longer(cols = all_of(sel_myokines_01),
               names_to = "indep_name",
               values_to = "indep_value") |> 
  group_by(dep_name, indep_name) |> 
  nest() |> 
  ungroup() |> 
  mutate(mod = map(data, 
                   ~lmer(dep_value ~ exam_order + indep_value + scale(bcm) + scale(age_m0) +
                           (1|immet_id), data = .x)),
         tidier = map(mod, ~tidy(.x, effects = "fixed", conf.int = TRUE))
  )

res_ass_activ_myokin_bcmADJ_01_tab <- res_ass_activ_myokin_bcmADJ_01 |> 
  unnest(tidier) |> 
  filter(str_detect(term, "indep_value")) |> 
  mutate(term = str_replace(term, "indep_value", indep_name),
         term = str_replace(term, "exam_orderm6", "examination")) |> 
  select(-data, -mod, -effect, - std.error, -statistic, -df) 

res_ass_activ_myokin_bcmADJ_01_fig <- res_ass_activ_myokin_bcmADJ_01 |> 
  inner_join(res_ass_activ_myokin_bcmADJ_01_tab |> 
               filter((p.value < 0.06) & (sign(conf.low) == sign(conf.high))) |> 
               select(dep_name,indep_name) |> 
               distinct(),
             by = c("dep_name","indep_name")) |> 
  mutate(fig = pmap(list(mod, dep_name, indep_name), ~plot_model(model = ..1,
                                                                 type = "emm",
                                                                 show.data = TRUE,
                                                                 term  =  c("indep_value", "exam_order")) +
                      labs(title = "Association betwen Myokines and Activity features",
                           y = ..2,
                           x = ..3) +
                      sjPlot::theme_sjplot2() +
                      ggplot2::theme(
                        plot.title    = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
                        plot.subtitle = ggplot2::element_text(size = 15, face = "italic"),
                        axis.title    = ggplot2::element_text(size = 12, face = "bold"),
                        axis.text     = ggplot2::element_text(size = 11, face = "bold")
                      )
  ))


# pdf(file.path("output/figures", 
#               paste0(format(Sys.Date(), "%y%m%d"),
#                      "_", 
#                      "assoc_activ_myokin_bcmADJ_01.pdf")), 
#     width = 10, height = 7)
# walk(res_ass_activ_myokin_bcmADJ_01_fig$fig, print)
# dev.off()

export(res_ass_activ_myokin_bcmADJ_01_tab, "output/tables/260131_ass_activ_myokin_bcmADJ_01.xlsx")

# **********************************************************************
# 8. response to trt - ordinal analysis ----
# **********************************************************************
res_resp_ordinal_01 <- d16_statistika_fi |> 
  filter(exam %in% c("m0", "m6")) |> 
  select(where(~ mean(is.na(.x)) <= 0.20)) |> 
  pivot_longer(cols = -c(immet_id, exam),
               names_to = "var_name",
               values_to = "var_value") |> 
  left_join(d11_clinics |> select(immet_id, disease_subtype, response_m0_m6, bmi), by = "immet_id") |> 
  left_join(age_m0, by = "immet_id") |> 
  group_by(var_name) |> 
  nest() |> 
  mutate(
    mod = map(data, ~ possCLMM(
      data = .x,
      response_m0_m6 ~ exam + var_value + bmi + age_m0 + (1 | immet_id),
      link = "logit"
    )),
    mod_m0 = map(data, ~ possCLMM(
      data = .x,
      response_m0_m6 ~ exam + bmi + age_m0 + (1 | immet_id),
      link = "logit"
    )),
    emms_lat = map(mod, ~ emtrends(.x, ~ var_value | exam, mode = "latent")),
    contrasts_lat = map(emms_lat, ~ pairs(.x, adjust = "tukey")),
    contrasts_lat_tbl = map(contrasts_lat, ~ broom::tidy(.x, conf.int = TRUE)),
    emms_prob = map(mod, ~ emmeans(.x, ~ var_value | exam, mode = "prob")),
    emms_prob_tbl = map(emms_prob, as.data.frame), 
    mod_diff = map2(mod_m0, mod, anova)
  )
