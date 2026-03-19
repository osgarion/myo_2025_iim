# Manual AI Predictors Evaluation 01

## Purpose

This manual describes a reproducible workflow for comparing repeated-measures mixed models with `exam` included as a fixed effect. The workflow evaluates whether activity-related GGIR parameters can:

- reflect a clinical outcome,
- complement an existing clinical measure,
- or serve as candidate surrogate markers.

The workflow is designed for longitudinal data with repeated patient visits and always includes `exam` as a fixed effect.

The manual is meant for an AI agent that should be able to recreate the task end-to-end in R.

## Important instruction before any replication

Before replicating or extending this workflow, the AI agent must ask the user:

1. Which statistical method should be used.
2. Which predictors should be treated as tested predictors.

This is mandatory, because the workflow structure is reusable, but the exact method and tested variables may change between tasks.

This manual does not describe any non-exam-adjusted workflow.

## Core modeling principle

The central goal is to compare multiple highly similar models where:

- the model family is the same across models,
- the covariate structure is the same,
- the random-effect structure is the same,
- the tested predictor changes one by one,
- the outcome is fixed within a given comparison block,
- `exam` is included as a fixed effect because the data are longitudinal and values visibly change over time,
- the patient identifier is included as a random intercept via `(1 | immet_id)`.

The standard mixed-model template is:

```r
outcome ~ tested_variable + exam + age_m0 + sex + (1 | immet_id)
```

For adjusted models with an additional base clinical predictor:

```r
outcome ~ base_predictor + tested_variable + exam + age_m0 + sex + (1 | immet_id)
```

## Why `exam` must be included

`exam` represents scheduled repeated visits such as `M0`, `M3`, `M6`, `M18`.

If `exam` is omitted:

- a tested predictor may partly behave as a proxy for time,
- beta estimates may absorb longitudinal trend instead of a true independent association,
- semi-partial R2 and marginal R2 may be artificially inflated,
- ranking of predictors may reflect shared time trend rather than independent explanatory value.

Therefore, the exam-adjusted workflows use:

```r
+ exam + age_m0 + sex + (1 | immet_id)
```

and not just:

```r
+ age_m0 + sex + (1 | immet_id)
```

## Why interaction with `exam` is not used by default

Do not add `tested_variable * exam` as the default screening model.

Reason:

- the main screening question is whether a tested predictor has an overall association with the outcome after adjustment for visit,
- adding an interaction changes the question to whether the association differs by visit,
- interaction terms produce several visit-specific coefficients instead of one focal effect,
- this breaks the simple ranking logic based on a single tested predictor term,
- the models become less stable and harder to compare,
- confidence intervals typically widen,
- interpretation becomes much harder for routine ranking workflows.

Therefore, the standard workflow uses:

```r
outcome ~ tested_variable + exam + age_m0 + sex + (1 | immet_id)
```

and not:

```r
outcome ~ tested_variable * exam + age_m0 + sex + (1 | immet_id)
```

Interaction with `exam` should only be considered as a secondary, focused analysis for a few top candidate predictors if the user explicitly wants to test whether the association changes over visits.

## Why `(exam | immet_id)` is not used by default

Do not use `(exam | immet_id)` as the default random-effects structure in these screening scripts.

Reason:

- the current workflow is intended as a stable ranking screen across many predictors,
- `(1 | immet_id)` is more robust and easier to compare across many models,
- `(exam | immet_id)` adds patient-specific visit slopes and is substantially more complex,
- with sparse or unbalanced repeated visits it can lead to singular fits or convergence problems,
- if `exam` is a factor with several levels, the random-effects structure becomes much heavier.

Therefore, the default structure remains:

```r
(1 | immet_id)
```

Use `(exam | immet_id)` only in a later, more focused confirmatory step if the user explicitly requests patient-specific visit trajectories and the data are sufficient for that model complexity.

## Data objects used in the project

### Main longitudinal data object

`d18_ggir_data`

Expected columns typically include:

- `immet_id`
- `exam`
- `age_m0` or `age`
- `sex`
- `ph_vas`
- `mmt8`
- `fi2`
- `sf36_pcs`
- `sf36_mcs`
- `sf36pf`
- `sf36rp`
- `sf36bp`
- `sf36gh`
- `sf36vt`
- `sf36sf`
- `sf36re`
- `sf36mf`
- GGIR activity parameters such as:
  - `ig_gradient_enmo`
  - `ig_intercept_enmo`
  - `acc_day_spt_wei`
  - `m5_enmo`
  - `mvpa_t100_time_min`
  - `p9931_enmo`

### TIS data object

`d18_tis_comp_01`

Expected columns used for merge:

- `immet_id`
- `exam`
- `tis_total`

For TIS workflows, export or select:

```r
select(immet_id, exam, tis_total)
```

and merge it into `d18_ggir_data` by:

```r
by = c("immet_id", "exam")
```

## Standard preprocessing rules

### Column resolution

The scripts should tolerate small naming differences and resolve the first available matching column, for example:

- `exam_order` or `exam`
- `age_m0` or `age`
- `fi2` or `fi_2`
- `mmt8` or `mmt8_total`
- `ph_vas`, `physician_vas`, or `PhVAS`

### Harmonized working dataset

The working dataset should standardize the columns used in modeling:

- `immet_id` as character
- `exam` as ordered factor
- `sex` as factor
- `age_m0` as numeric
- outcome as numeric
- tested predictors as numeric

If `exam` contains labels like `M0`, `M3`, `M6`, `M18`, order them numerically by the visit number.

## Tested predictor groups

### GGIR activity parameters

The standard activity list used in the workflows is:

```r
c(
  "ig_gradient_enmo",
  "ig_intercept_enmo",
  "acc_day_spt_wei",
  "m5_enmo",
  "mvpa_t100_time_min",
  "p9931_enmo"
)
```

### Clinical predictors used in GGIR-outcome models

```r
c(
  "mmt8",
  "fi2",
  "ph_vas",
  "sf36_pcs",
  "sf36_mcs",
  "sf36pf",
  "sf36rp",
  "sf36bp",
  "sf36gh",
  "sf36vt",
  "sf36sf",
  "sf36re",
  "sf36mf"
)
```

### Single-predictor comparison set

For workflows comparing standalone predictors, the standard list is:

```r
c(
  "fi2",
  "mmt8",
  "ig_gradient_enmo",
  "ig_intercept_enmo",
  "acc_day_spt_wei",
  "m5_enmo",
  "mvpa_t100_time_min",
  "p9931_enmo"
)
```

## Supported workflow variants

### Variant A: PhVAS workflows

This combines the ideas of the original PhVAS comparison workflows into one exam-adjusted implementation.

#### A1. Base-adjusted comparison

Outcome:

- `ph_vas`

Models:

```r
ph_vas ~ fi2 + tested_variable + exam + age_m0 + sex + (1 | immet_id)
ph_vas ~ mmt8 + tested_variable + exam + age_m0 + sex + (1 | immet_id)
```

Tested variables:

- the GGIR activity parameters

#### A2. Single-predictor comparison

Outcome:

- `ph_vas`

Models:

```r
ph_vas ~ tested_variable + exam + age_m0 + sex + (1 | immet_id)
```

Tested variables:

- `fi2`
- `mmt8`
- the GGIR activity parameters

Output folder:

- `output/exam_incl/ph_vas`

### Variant B: SF-36 workflows

For each SF-36 outcome, run three exam-adjusted blocks:

- `fi2`-adjusted GGIR comparison
- `mmt8`-adjusted GGIR comparison
- single-predictor comparison

Outcomes:

```r
c(
  "sf36_pcs", "sf36_mcs", "sf36pf", "sf36rp", "sf36bp",
  "sf36gh", "sf36vt", "sf36sf", "sf36re", "sf36mf"
)
```

Models:

```r
sf36_xxx ~ fi2 + tested_variable + exam + age_m0 + sex + (1 | immet_id)
sf36_xxx ~ mmt8 + tested_variable + exam + age_m0 + sex + (1 | immet_id)
sf36_xxx ~ tested_variable + exam + age_m0 + sex + (1 | immet_id)
```

Output folder:

- `output/exam_incl/sf36`

### Variant C: TIS workflows

Merge `tis_total` into `d18_ggir_data`, then run:

- `fi2`-adjusted GGIR comparison
- `mmt8`-adjusted GGIR comparison
- single-predictor comparison

Outcome:

- `tis_total`

Models:

```r
tis_total ~ fi2 + tested_variable + exam + age_m0 + sex + (1 | immet_id)
tis_total ~ mmt8 + tested_variable + exam + age_m0 + sex + (1 | immet_id)
tis_total ~ tested_variable + exam + age_m0 + sex + (1 | immet_id)
```

Output folder:

- `output/exam_incl/tis`

### Variant D: GGIR-as-outcome workflows

For each GGIR parameter used as an outcome, compare the clinical predictors one by one.

Outcomes:

```r
c(
  "ig_gradient_enmo",
  "ig_intercept_enmo",
  "acc_day_spt_wei",
  "m5_enmo",
  "mvpa_t100_time_min",
  "p9931_enmo"
)
```

Tested predictors:

- `mmt8`
- `fi2`
- `ph_vas`
- all SF-36 parameters

Model:

```r
ggir_parameter ~ tested_variable + exam + age_m0 + sex + (1 | immet_id)
```

Output folder:

- `output/exam_incl/ggir`

## Statistical method

The project uses linear mixed-effects models:

```r
lme4::lmer(..., REML = FALSE)
```

Reasons:

- same family across all compared models,
- stable AIC comparison,
- compatibility with repeated visits,
- compatibility with patient-level random intercept.

Do not use `REML = TRUE` for these ranking workflows because AIC comparison across models with different fixed effects should be based on `REML = FALSE`.

The default method for the current project is therefore:

```r
lme4::lmer(..., REML = FALSE)
```

However, the AI agent still must ask the user to confirm the desired statistical method before reproducing the task.

## Metric extraction

Each model comparison must evaluate the tested predictor only.

This is crucial: after adding `exam`, the extraction must still target the tested predictor, not the `exam` term.

### Correct extraction logic

Use:

```r
filter(.data$Parameter == focal_term)
```

where `focal_term` is explicitly the currently tested predictor name.

Do not extract:

- `exam`
- any `exam` contrast
- any covariate other than the tested predictor

This control is especially important after adding `exam`, because mixed-model summaries may contain several `exam` coefficients. None of those may replace the tested predictor in the evaluation table or radar graph.

### Metrics to extract

For each tested predictor:

1. Standardized beta
   - use `parameters::standardize_parameters(method = "refit")`
   - extract the focal term only
   - use absolute beta for scoring

2. Standardized beta 95% CI
   - `CI_low`
   - `CI_high`
   - `ci_width = abs(CI_high - CI_low)`

3. CI precision
   - use the non-binary support proportion approach
   - recommended formula:

```r
ci_precision_raw = (1 / ci_width) * support_prop
```

4. Semi-partial R2
   - for the tested predictor
   - in mixed models, compute it as:

```r
R2_marginal(full_model) - R2_marginal(reduced_model_without_tested_predictor)
```

5. R2 marginal
   - from `performance::r2_nakagawa()`

6. R2 conditional
   - from `performance::r2_nakagawa()`

7. AIC
   - lower is better

## Scoring and ranking

All metrics must be converted to 0-1 scale.

Use a helper similar to:

```r
scale01(x, higher_better = TRUE)
```

Rules:

- higher is always better after scaling,
- if all values are identical, return `0.5`,
- AIC must be reverse-scaled.

### Radar/composite metrics

Use these scoring axes:

- `score_beta`
- `score_ci`
- `score_semi_partial_r2`
- `score_r2_marginal`
- `score_r2_conditional`
- `score_aic`

### Composite score

Create:

- `composite_score`

as the mean of the radar metrics, unless the user explicitly requests weighted scoring.

Use `composite_score` for ranking.

The ranking objective is to compare tested predictors, not visit effects. Any implementation must therefore preserve focal extraction and focal scoring strictly on the tested predictor.

## Radar graph requirements

Use `fmsb::radarchart()`.

Each polygon represents one tested predictor.

### Aesthetic requirements used in the project

- stable color per variable name across graphs
- same variable must keep the same color even if order changes
- named color mapping must be used, not position-based mapping
- titles and legends should be large and readable
- exported image should use `png(..., type = "cairo")`

### Standard radar axes

- `|beta| std.`
- `CI precision`
- `R2 - semi-partial`
- `R2 marginal`
- `R2 conditional`
- `AIC (reverse)`

## Tables

Each workflow must produce:

### Composite table

Include at least:

- `rank`
- `parameter`
- `response`
- `beta_std`
- `ci_low`
- `ci_high`
- `ci_width`
- `ci_support_prop`
- `semi_partial_r2`
- `r2_marginal`
- `r2_conditional`
- `aic`
- `score_beta`
- `score_ci`
- `score_semi_partial_r2`
- `score_r2_marginal`
- `score_r2_conditional`
- `score_aic`
- `composite_score`

### Raw metrics table

Include the same metrics before export formatting, typically with slightly higher precision.

### Export formatting

Round numeric values for the export view, usually to 3 decimals for the composite table and 4 decimals for raw metrics if needed.

## Export structure

### Root folder

All exam-adjusted outputs go under:

```text
output/exam_incl
```

### Required subfolders

- `output/exam_incl/ph_vas/tables`
- `output/exam_incl/ph_vas/figures`
- `output/exam_incl/sf36/tables`
- `output/exam_incl/sf36/figures`
- `output/exam_incl/tis/tables`
- `output/exam_incl/tis/figures`
- `output/exam_incl/ggir/tables`
- `output/exam_incl/ggir/figures`

If a folder does not exist, create it automatically.

### File naming

Use names beginning with the date:

```text
YYMMDD_short_description_suffix_01
```

Examples:

- `260318_ph_vas_phvas_fi2adj_exam_01.xlsx`
- `260318_tis_total_tis_single_exam_01.png`
- `260318_sf36pf_sf36_mmt8adj_exam_01.xlsx`
- `260318_ig_gradient_enmo_ggir_exam_01.png`

## Required script outputs in R

Each script should:

- create a main workflow result object,
- print radar plot(s),
- print export table(s),
- export the files,
- print or message the output directory path.

Use distinct object names so workflows do not overwrite each other unnecessarily.

Examples:

- `ggir_vis01_exam_results`
- `ggir_vis03_exam_results`
- `ggir_vis04_exam_results`
- `ggir_vis05_exam_results`

The printed preview in RStudio should come from standard `print()` calls:

- radar plots should appear in the `Plots` pane,
- formatted tables should appear in the console,
- export should happen in parallel as file output.

## Quality control checklist

Before considering the task finished, the AI agent must verify:

1. `exam` is present in every intended mixed model.
2. The random structure remains `(1 | immet_id)`.
3. The tested predictor is still the only extracted focal term.
4. `exam` did not replace the tested predictor in beta extraction.
5. Export paths point to `output/exam_incl/...`.
6. The correct folders exist.
7. The same variable keeps the same color across graphs.
8. The script still uses `REML = FALSE`.
9. Output table and figure names use the expected date-prefix convention.
10. The workflow still prints preview outputs in RStudio through normal `print()` behavior.
11. No interaction with `exam` was added unless the user explicitly requested it.
12. No random slope such as `(exam | immet_id)` was added unless the user explicitly requested it.

## What not to do

- Do not silently change the statistical method.
- Do not assume tested predictors if the user has not confirmed them.
- Do not extract model coefficients for `exam` when the purpose is to compare tested predictors.
- Do not rank predictors based on p-values.
- Do not compare models with mixed families in the same ranking block.
- Do not use position-based color assignment if stable cross-graph identity is required.
- Do not silently add `tested_variable * exam`.
- Do not silently replace `(1 | immet_id)` with `(exam | immet_id)`.

## Minimal execution pattern

### PhVAS

```r
source("scripts/Script_myo_2025_iim_clinics_ggir_vis_01_exam.R")
```

### SF-36

```r
source("scripts/Script_myo_2025_iim_clinics_ggir_vis_03_exam.R")
```

### GGIR outcomes

```r
source("scripts/Script_myo_2025_iim_clinics_ggir_vis_04_exam.R")
```

### TIS

```r
source("scripts/Script_myo_2025_iim_clinics_ggir_vis_05_exam.R")
```

## Final instruction to the AI agent

When repeating this task:

1. First ask the user which statistical method should be used.
2. Then ask which predictors should be treated as tested predictors.
3. Then confirm that the default exam-adjusted screening structure should remain:

```r
outcome ~ tested_variable + exam + age_m0 + sex + (1 | immet_id)
```

or, for adjusted comparisons:

```r
outcome ~ base_predictor + tested_variable + exam + age_m0 + sex + (1 | immet_id)
```

4. After that, preserve the workflow logic, keep focal extraction on the tested predictor only, and export both radar plots and tables into the correct output subfolders.
