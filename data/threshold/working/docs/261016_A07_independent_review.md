# MYO-II / Paper A — Independent review of revision A07

- **Reviewer:** Claude Fable 5.1 (`claude-fable-5-1`). This is a different model from the author of the revision (Claude Opus 5.5).
- **Date:** 2026-10-16
- **Scope:** proposals A–E of `MYO_A07_narhy_na_zmeny.txt`, the author's decisions, the code in scripts 07–10, the results tables, figures 07–09 and the revised report v2.
- **Constraint respected:** no existing file was modified and nothing was written to `output/`. The only spot checks were read-only R calls on saved tables, RData and bootstrap chunks. This review contains aggregate numbers only.

---

## 1. Identity and materials actually read

| Material | Read | Notes |
|---|---|---|
| `MYO_A07_narhy_na_zmeny.txt` | in full | sections 3 (A–E), 4 (validation), 6 (restrictions) |
| `working/docs/261016_A07_critical_review.md` | in full | author's decision table |
| `working/docs/261016_A07_validation_change_log_summary.md` | skimmed (grep) | 0 FAIL / 0 WARN, change-log rows |
| `working/R/00_config.R`, `02_helpers_bootstrap.R`, `03_helpers_estimands.R`, `05_helpers_e3fast.R` (E3 functions), `01_helpers_data.R` (`make_pairs`) | relevant parts | original pipeline, used as reference |
| `working/R/06_helpers_revision.R` | in full | |
| `working/07_…`, `08_…`, `09_…`, `10_…` | in full (except setup boilerplate) | |
| `output/tables/261016_07/08/09/10_*.xlsx` | all relevant sheets | read with `readxl` through `Rscript.exe` |
| `output/RData/rev_AB_objects.RData` and `boot_chunks/rev_B_*` | spot check | bootstrap bias of `share_pos` |
| Figures `261016_07_B1`, `261016_08_C1` | viewed | other figures checked through their tables (R10) |
| `reports/261016_thresholds_report_v2.Rmd` | revised sections in full | lines 150–410, 576–610, 820–1310, 1464–1616 |
| `reports/261016_thresholds_report_v2.html` | text extracted | compared rendered numbers with the Rmd and the tables |
| `261013_phase3_disposition_and_locked_plan.md`, original report | headings and diff only | context |

---

## 2. Overall verdict: **MINOR REVISION**

- The code is correct in every place I checked. The original results are reproduced exactly and no population or estimand was changed without documentation.
- Accepting A, the B estimand and the C/D analyses is statistically justified. E is set up reasonably and labelled as exploratory.
- **No critical findings.**
- **Three major findings concern interpretation, not computation:**
  - the MMT8 "attenuation" claim (F1);
  - the meaning of the equal-share reference and the "most concentrated at M3 high intensity" claim (F2);
  - the "df = 3 is the best compromise / verdict holds for all calibrated specifications" claim (F3).
- Each can be fixed by rewording the report text. None needs a re-run, apart from an optional extra simulation cell for F3.

---

## 3. Findings

Severity: **critical** = wrong result or invalid conclusion; **major** = a conclusion stated in the report that the numbers do not support, or a material methodological gap; **minor** = wording, presentation, documentation or a secondary methodological point.

| # | Sev. | Location | Problem | Evidence (what I checked) | Recommended fix |
|---|---|---|---|---|---|
| F1 | **major** | Report: Key message 6; §"Subtypes adjusted…" Evaluation (bullet 3: "attenuated… Part of the smaller IMNM gain may therefore reflect their greater baseline muscle weakness"); "How to use this" ("its partial dependence on baseline muscle strength"); H5 ("partly explained by greater baseline muscle weakness"); change-log row "Subtypes"; `validation_change_log_summary.md` row 53 | The "attenuation by MMT8" is a cross-sample comparison. Against the ph_VAS-adjusted estimate (−20.2, n = 47) it looks like attenuation. **Within the MMT8 sample**, IMNM − DM ΔT100 moves only from −18.8 (95 % CI −36.5 to −1.3, unadjusted) to −17.5 (−36.6 to 3.0, adjusted for MMT8), a change of 1.3 min/day (≈ 7 %). MMT8 is not associated with ΔT100: the coefficient is 0.87 min/day per 5 points (−3.4 to 6.6) and the LRT p = 0.71. Crossing 0 is mostly widening of the CI, not a shift of the estimate. The design rows 1 vs 2 were meant to separate sample restriction from adjustment, yet the text compares across samples. | `09` sheet `models`: rows E_a/E_b, severity mmt8, mean M3–M18, IMNM − DM; sheet `lrt` row "severity main effect", dT100, mmt8 | Replace with: "In the MMT8 sample the difference is similar unadjusted (−18.8) and adjusted (−17.5); its CI includes 0 after adjustment. MMT8 itself is not associated with ΔT100. The data therefore do not show that baseline weakness explains the IMNM deficit." Drop "partly explained" from H5, KM6 and the next-study bullet, or reformulate it as an untested hypothesis. |
| F2 | **major** | Report: KM5; F6 Evaluation ("far more than equal sharing… **everywhere**", "clearly higher above ~240 mg at M3", "the early high-intensity gain is the most concentrated"); H3 ("2–4 times"); change log row F6 | (a) The reference `min(k,n₊)/n₊` is the **mathematical minimum** of `share_pos`. Any non-identical set of positive gains exceeds it. "Above equal sharing everywhere" / "unequal at every intensity" is therefore true by construction, not an empirical finding. It needs a realistic reference, e.g. the share expected under a stated continuous distribution, or a Gini/Lorenz summary. (b) "Clearly higher above ~240 mg at M3" and "most concentrated early" rest on point estimates with overlapping **pointwise** CIs (M3 325–350: 65 [45–77] vs M6: 50 [36–62]). No formal comparison was made, and at 240–260 mg M6 (53.7) ≥ M3 (52.8). (c) The text summarises only ≥ 60 mg. The highest shares in the whole table are at 20–40 mg (M6 20–30: 70 %; M18: 65 %; M3: 61 %), which is visible in the figure but not mentioned. (d) `share_pos`/ref actually ranges from 2.0 to 5.1 (M3 325–350: 65.3/12.8), not "2–4". (e) On the absolute-minutes scale, concentration at high intensity is partly driven by between-patient heterogeneity in baseline level, because minutes at > 250 mg are highly skewed. This alternative explanation is not mentioned. | `07` sheet `B_concentration` (all 75 rows); figure `07_B1` | (a) Present the grey line as a lower bound ("the minimum possible share"), not as "equal-sharing expectation", and delete "everywhere/at every intensity" as a finding. (b) Downgrade to "point estimates are highest at M3 above ~280 mg; CIs overlap and no formal comparison was made". If the claim is to be kept, bootstrap the M3 − M6 difference in mean `share_pos` over the high-intensity bins. (c) Mention the low-bin values and their higher reference. (d) Write "2–5×". (e) Add one sentence on scale dependence. Optionally report `share_pos` of log-ratio changes as a sensitivity analysis. |
| F3 | **major** | Report §"Is the E3 verdict robust…" ("df = 3 is the best compromise", "holds for every specification that is calibrated (df ≥ 3)"); Limitations 4; Scientific summary 2; change-log row "E3 robustness" | "Calibrated" is defined under a true null equal to the **df-3 fit itself**, so df 3 is favoured by design. Under the misspecified null (true smooth = df 5), any-bump rates are: df 3 = 5.4 / **18.2** / 7.6 % and df 4 = 1.6 / **23.8** / 1.0 % (M3/M6/M18). Neither is calibrated at M6. Under the df-3 null, df 4 is as well calibrated as df 3 (0.4–0.8 %) and has **higher** power for the narrow +0.1 bump (2 bins: 89/89/58 % vs 73/79/37 %). df 3 wins only for broad bumps (6 bins +0.1: 74/94/71 vs 33/24/14 %). The `ns(df)` spaces are not nested (quantile knots), so "more flexible" is not monotone either. The **substantive E3 conclusion is unaffected**: false-alarm inflation cannot create a negative finding, and df 3, 4 and 5, width weighting and both grids all give "no localised departure". | `08` sheet `C_df_simulation` (72 rows); `C_variant_decisions` | Reword to: "df 3 trades sensitivity to broad departures against narrow ones; df 4 is equally defensible. No departure is flagged with df 3–5, so the conclusion does not depend on this choice." Replace "calibrated (df ≥ 3)" with "df 3–5". Optional: add the power cells under df_true = 5, and a negative-bump scenario (the only flags seen, with df 1–2, were negative runs at 50–70 mg). |
| F4 | minor | D: `shape_contrast()` comment and report ("the clr level is not identifiable") | The level of the difference curve over the 19 bounded bins **is** identified. It is the change in the geometric mean of the 20–300 mg parts relative to the `<20` and `≥300` parts, which is a real compositional change. Excluding it is a legitimate *choice* (testing shape within 20–300 mg), but the stated reason is wrong. | code review of `06_helpers_revision.R` D block | Reword to "the test deliberately ignores a constant shift, which reflects exchange with the <20 / ≥300 parts". Optionally report the intercept contrast separately. |
| F5 | minor | `shape_contrast()` slope p-value | `2·min(P(T* − t ≥ |t|), P(T* − t ≤ −|t|))` takes the **thinner** of the two opposite tails. Under skew it is smaller than the usual `P(|T* − t| ≥ |t|)`, i.e. anti-conservative. The impact is nil here (all Holm p ≥ 0.32). | code review | Use `mean(abs(slope_b - slope) >= abs(slope))`. |
| F6 | minor | D sample definition | Patients without M0 are excluded from the M3/M6/M18 contrasts, although M0 cancels in `Δ = (clr_b − clr_0) − (clr_a − clr_0)`. This is defensible for consistency with E3, but it is an unnecessary loss of patients. | code review | Either state the reason or add a sensitivity analysis without the M0 requirement. |
| F7 | minor | D report text | The common-sample M6 − M3 Wald test (p = 0.024, **Holm 0.071**) is borderline. The report says the common sample "gives the same picture" and cites the raw p only. | `08` sheet `D_shape_summary` | Report the Holm value for the common sample too. |
| F8 | minor | E: LRTs | ML LRTs with 47 patients are anti-conservative. The "any subtype effect" p = 0.036 (ΔT100, ph_VAS) and the MMT8 severity → ΔIG result (LRT p = 0.016 while the bootstrap CI includes 0) show the mismatch between LRT and bootstrap. | `09` sheets `lrt`, `models` | Add a parametric-bootstrap LRT (e.g. `pbkrtest::PBmodcomp`) or Kenward–Roger F tests, or state the limitation explicitly. |
| F9 | minor | E: model E_c output | In `chg ~ … + sev * visit` with treatment contrasts, `b["sev"]` is the severity slope **at M3**, but it is labelled "severity (per unit)", visit "all". | `fit_subtype_models()` | Relabel it, or report the visit-averaged slope. |
| F10 | minor | Report E Evaluation ("ASyS gains slightly less IG than IMNM… which disappears after adjustment") | Restricting the sample alone already takes the CI over 0: −0.082 (−0.181 to 0.005) in the ph_VAS sample, unadjusted. The text attributes this to adjustment. | `09` `models`, dIG, ASyS − IMNM | Attribute it to sample restriction plus adjustment, using rows 1 vs 2 vs 3. |
| F11 | minor | Report KM6 and E construct box | "ph_VAS reflects **mostly** extramuscular activity in DM and ASyS" is an overstatement. ASyS: ρ with muscle VAS 0.70 vs extramuscular 0.87; DM: 0.38 vs 0.51; n = 14–18. The construct argument itself (IMNM ρ_muscle 0.97, median extramuscular VAS 2 mm) is sound. Also, "MMT8 and FI2… ρ ≈ 0.85" should be 0.82. | `09` `construct_ph_vas`, `collinearity_spearman` | Write "substantially driven by extramuscular activity". Correct ρ to 0.82. |
| F12 | minor | Report KM6, E Evaluation, H5 | "The shape of change does not differ between subtypes" reads as evidence of absence. The CIs are wide (tilt IMNM − DM −0.03, −0.33 to 0.28). | `09` `models` | Write "no evidence of a difference (wide CIs)". |
| F13 | minor | B: net/gross 0.20 flag; F6 text "That figure used the net gain as denominator, which is inflated…" | (a) The 0.20 cut-off is arbitrary. At net/gross = r the net share can reach up to (1 + r)/(2r), i.e. 300 % at r = 0.2, so "interpretable" cells can still exceed 100 %. (b) Above 200 mg every cell has net/gross 0.49–0.85 and is flagged "interpretable" by the author's own rule. There the net share (46–75 %) and `share_pos` (35–65 %) differ only modestly. Saying the first report's figure was "inflated" and is "no longer interpreted" is stronger than the data support. | `07` `B_concentration` | State the bound. Say that above 200 mg the original figure was moderately higher, not invalid. |
| F14 | minor | Report: discrepancy-log note ("small in absolute minutes and concentrated in a minority of patients") | This contradicts the revised H3 ("not specifically at high intensities", "not a responder subgroup"). | Rmd §Discrepancy log | Align the wording with H3. |
| F15 | minor | Validation (script 10) | Gaps against request §4: (a) R8 compares the E3 patient-set hash produced by the **same code twice** (tautological). It should compare against IDs/hashes saved by 01–06 (e.g. `qc_objects`). (b) R7 LOO covers only T100/T200 concentration and the linear tilt. It does not cover bin-level `share_pos` at high intensity (k = 4–6, most influence-prone), which is the basis of the "M3 high-intensity" claim (F2), nor the D Wald test. (c) R10b claim 6 uses the **lower** MC-CI bound ≤ 5 %, which is a weak criterion. (d) No check of the numeric claims in the text ("2–4×", "ρ ≈ 0.85", "attenuated"), which are exactly where F1/F2/F11 arise. (e) No convergence/singularity check of the lmer fits (bootstrap reports 0 failures, but singular fits are not counted). (f) No check of the plasmode construction (e.g. that the population residual is 0 under the null). | code review of `10_REV_VALIDATION_BOTH.R`; `10` tables | Add (a), (b) and (d) as a minimum. Change (c) to the point estimate ≤ 5 %. |
| F16 | minor | Rendered HTML | The HTML does not match the current Rmd. E Evaluation in the HTML reads "Patients with higher baseline **activity or weaker muscles** show a larger tilt…", while the Rmd (line 1217) reads "higher baseline ph_VAS… for MMT8 the association points the same way…". "Baseline activity" is also ambiguous (disease activity vs physical activity). | text extraction from HTML vs Rmd | Re-knit after the fixes. |
| F17 | minor | Report §Reproducibility | It still says "Run the scripts in order 01 → 06, then knit", but v2 needs 07–10 (and the clinical file). | Rmd ~line 1570 | Update to 01 → 10, and name the clinical file / env var. |
| F18 | minor | Author's review E4 | "about 115 observations" — the models use 126–127 (126 ph_VAS, 127 MMT8). | `09` `lrt` | Correct it. |

---

## 4. Items checked and found correct

### A — descriptive profile
- `density = value × 10 / width` holds per recording (`stopifnot` < 1e-8; QC D11 2.8e-14). Widths are 10 mg × 10, 20 mg × 9 and 25 mg × 4 (R2).
- The equivalence claim is correct: quantiles are equivariant to a positive constant, and the width is constant within a bin. Normalising per recording or after summarising is therefore identical (`norm_equals_scaled_raw` TRUE for all cells).
- The raw profile reproduces the original 01 table in 92/92 cells.
- The 120–140 mg rise is a width artefact. Raw step +2.9 to +4.9 min/day; normalised step −2.8 to −3.4. The normalised medians are monotone at every visit, and the small raw rise at 300–325 mg at M6 also disappears.
- No interpolation is used. Decision **accept** is justified.

### B — concentration
- `k = ceiling(0.1·n)` over all pairs of the cell (5 / 6 / 4 at M3 / M6 / M18), selected per bin × visit. Ties at the k-th value do not change the sum. If `n₊ < k`, all positives are summed (100 %).
- `share_pos` ∈ [0, 100]. The original C3 is reproduced in 69/69 cells, and there is an independent base-R recomputation (R4, max diff 0).
- Undefined net-share cells (M3 20–40, M18 30–60 mg) are correctly explained as `Σ Δ ≤ 0`. This answers request item B8.
- Bootstrap: duplicated patients are relabelled (`resample_patients`). My spot check of 2000 replicates found the **bias of `share_pos` is small and not systematic** (mean bootstrap − observed −0.36 pp, range −3.4 to +3.0 pp over 75 cells). Percentile CIs are therefore usable as approximate, and the report's caveat is adequate.
- Rendered numbers match the table: 51–83 %, 32–65 %, 32–47 %, 64 (42–77) % at M3 300–325, and T100 37/41/35 %.

### C — E3
- Closure is correct (max |Σclr| 5e-15; replaced composition = 1440).
- Zeros: 4/217 recordings, max 3.2 % in a part.
- The **invariance argument is correct**, including zero replacement:
  - per-bin constants (log 10/w) cancel in the paired difference;
  - per-record constants (clr centre, multiplicative factor `1 − Σrepl/1440`) are constant across parts and lie in the spline column space via the intercept;
  - the replacement value (0.65 × one epoch per valid day) does not depend on bin width, so applying it on the raw-minutes scale is right.
- The numeric check (≤ 3e-16) is by construction, because normalisation is applied after replacement. The argument would fail for zero parts if one normalised *before* replacing; with 4/217 recordings this is immaterial.
- Same basis per visit; native grid rebuilt through `e3_comp_grouped()` equals the original exactly.
- The variant classes are sensible: width-weighted = same estimand; 20-mg grid and to-400 = changed estimand. They are labelled as such and quantile re-binning is correctly rejected.
- The matrix engine reproduces the primary residuals (3e-16). Band half-widths are within ±4 % of the primary.
- Plasmode design: true null = df_true fit plus observed patient covariance, resampled. This is sound for its purpose.
- The rewritten panel A/B explanations and the "absence ≠ equal informativeness" box meet request C.

### D — shape contrasts
- The bootstrap Wald on the three non-intercept `ns(df = 3)` coefficients is a valid test that the smooth difference curve is constant (zero coefficients ⇔ constant fit). The covariance and the bootstrap reference distribution come from the same 2000 replicates.
- Bin-wise max-t band computed correctly. Holm is applied within each sample across the three pairs. M6 − M3 reproduces the original C1 exactly.
- The season/timing table is correct and the structural confounding argument is right: 0 % same season as M0 at M6/M18; M3 vs M6 2 %, M6 vs M18 71 % same season.
- The report correctly states "not significant after Holm" (0.111) and uses descriptive wording for panel A.
- M18 − M3 at 250–300 mg is −0.22 to −0.27, matching the "−0.2 to −0.3" in the text.

### E — subtypes
- Linkage 71/71, subtype codes agree 71/71. Clinical "x"/blank → NA. The clinical file md5 is logged.
- Availability and distributions are correct. ph_VAS is missing mostly in ASyS (14/17 vs 18/18 DM).
- The model is parsimonious (≤ 12 fixed effects, 126–127 obs, 47–48 patients, 12–18 per cell).
- Severity is centred at patient level, and subtype contrasts do not depend on the centring (no subtype × severity term).
- E_e is not fitted for tilt, which is correct. LRT nesting and df are correct (4, 6, 2, 1).
- The cluster bootstrap resamples patients with their covariates (0 failed replicates).
- FI2 is not combined with MMT8 (ρ 0.82).
- No subtype × visit interaction p < 0.05 (min 0.081): correct. The severity-effect numbers in the text match the table.
- The wording "exploratory / hypothesis-generating / no responder subgroup" is appropriate apart from F1, F11 and F12.

### Validation and restrictions
- R1–R6 and R9/R9b are meaningful. Primary observed E1–E4 are reproduced with max diff 0 in all 8 estimand groups, primary E3 bands are re-derived exactly, and results are reproducible from 1 to 4 workers.
- **Section 6:**
  - no silent change of populations: pairs 47/52/39 are unchanged and the E subsets are documented;
  - thresholds, bin definitions and primary model specifications are unchanged;
  - primary results are unchanged;
  - new estimands (B `share_pos`, D contrasts, E models) are declared as new;
  - the original report and code are preserved;
  - the replacement of the F3/F6 figures in the report body is documented in the change log.
- The only restriction-related issue is documentation (F17).
