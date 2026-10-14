# MYO-II / Paper A — Phase 3: independent methodological review

- **Reviewer model identity:** Claude **Fable 5.1** (model ID `claude-fable-5-1`), running as a Claude Code subagent. This is a different model from the author of the reviewed document (Claude Opus 5.5).
- **Date:** 2026-10-13
- **Object of review:** `261013_phase1_audit_plan_ABC.md` (sections A–D)
- **Role:** independent critical biostatistical reviewer. I was not asked to agree and I did not modify any existing file.

---

## 1. Reviewer identity and materials actually read

| Material | How it was used |
|---|---|
| `working/docs/261013_phase1_audit_plan_ABC.md` | Read in full |
| `working/docs/data_dictionary_master_v1.csv` | Read in full |
| `working/docs/repro_r05_baseR_check.txt` | Read in full. Note: the file ends with an R error (`object 'visit_code.x' not found`), so the reproduction script did not finish cleanly after the per-bin tables. |
| `MYOII_thresholds_agent_handoff_prompt.md` | Read in full |
| `MYOII_ACC_OUTCOMES_MASTER_V1.rds` | Analysed independently in base R (`Rscript --vanilla`). Scripts are in the session scratchpad and were not added to the project. No patient-level rows are reproduced here. |
| `A04_FINAL_CUMULATIVE_THRESHOLD_ANALYSIS_v2_3.R` | Lines ~150–215 (seasonality block) |
| `A07_FINAL_EXCLUSIVE_BIN_ANALYSIS_v2_3.R` | Searched only, for the same `Chi Df` pattern |
| `PaperA_FINAL_A04_CUMULATIVE/A04_FINAL_cumulative_analysis.xlsx`, `sessionInfo.txt` | Sheet `season_anchor_test` read from the unzipped XML. lme4 version taken from sessionInfo. |
| `R01_R04/PaperA_R03_PARALLEL_CLUSTER_BOOTSTRAP_BOTH.R` | Read in full |
| `R01_R04/PaperA_R00_helpers.R`, `R01_R04/R00_MYOII_common.R` | Path-handling lines only |
| `MYOII_R05_R06_revised…/R00_MYOII_common.R`, `R05_…R`, `R06_…R` | R00 header, R05 bootstrap block, R06 RNG and standardization lines |

**Not read:** the R01, R02 and R04 bodies, the `.docx` summaries, `RULE_rmd_report_style.md`, and `older_dataset/`.

**Not run:** lme4 and glmmTMB do not load in this environment, so I fitted no mixed models. All my numbers come from complete pairs in base R, with a patient-level paired bootstrap (B = 2000, one draw shared by all bins within a visit) and a max-t band.

---

## 2. Verdict

**REVISE**

**What the audit gets right:**
- The data-integrity audit is solid. I reproduced D3–D5, D8–D10, D13 (pooled), D18 and K3, and V1 and V2 hold up.
- Moving away from sliding windows and argmax summaries is the right direction.

**Why the plan cannot be locked yet:**
1. **The central scientific finding K1 is mis-stated.** On my reanalysis, the high-intensity pattern is not an artefact of SD standardisation. It persists, and in fact becomes stronger, on the relative scale that the plan itself proposes as co-primary.
2. **The absolute-minutes estimand depends on the grid.** It varies with bin width, so the reported 120–140 mg "absolute peak" is an artefact of that bin being twice as wide.
3. **The decision rule would recover the 80–200 mg band by construction.** Its coverage gate and its ≥1 min/day gate both fall off monotonically with intensity. That is exactly the post-hoc band the brief forbids.
4. **The plan misses several things that change interpretation:**
   - wrist-side switches between visits;
   - the time-use (compositional) structure of the bins;
   - the upper end of the grid at 400 mg acting as a censoring edge;
   - day-level GGIR output that would make reliability estimable.

---

## 3. Major objections

### MO1. K1 ("SD standardisation drives the high-intensity peak") is not supported by the data

**Issue.** I recomputed paired M0→Mx changes per bin (M3 n = 47, M6 n = 52, M18 n = 39 pairs).

At **M3**:

| Bin (mg) | Δ pooled-SD | dz = mean Δ / SD(Δ) | Ratio of means Mx/M0 | Geometric ratio (c = 0.1) |
|---|---|---|---|---|
| 100–110 | 0.37 | 0.50 | 1.20 | 1.19 |
| 180–200 | 0.51 | 0.64 | 1.52 | 1.44 |
| 260–280 | 0.59 | 0.48 | 1.94 | 1.66 |
| 300–325 | 0.56 | 0.37 | 2.16 | 1.58 |

At **M6 and M18**, the ratio of means rises monotonically with intensity, to 1.93 and 1.66 at 375–400 mg.

The time share of each bin within 20–400 mg (a compositional view) also rises monotonically with intensity at M3: 0.93 at 20–30 mg, 1.11 at 100–110 mg, 1.79 at 260–280 mg, 1.96 at 375–400 mg.

**Why it matters.** Three things follow:
- **Pooled-SD standardisation does not create the high-intensity response. It flattens it.** The relative change keeps rising across the high bins, while Δ/SD stays roughly flat at about 0.5 from 120 to 400 mg.
- **The paired signal-to-noise ratio (dz) has a broad plateau of 0.4–0.64 from about 90 to 400 mg.** It peaks at 180–200 mg (M3), 140–160 mg (M6) and 160–180 mg (M18), not at 260–280 mg.
- **The absolute scale localises change in the middle of the spectrum** only because that is where most of the time is spent.

So the three scales answer three different questions; none of them is an artefact. The audit's suspension of H1 is therefore justified for the wrong reason. The correct statement: "M3 shows the largest relative increase at high intensities, but that increase is a fraction of a minute per day, and it is concentrated in a minority (see MO6)."

**Fix.**
- Rewrite K1 as "the location of the maximum depends on the estimand; the estimands must be chosen a priori for distinct questions."
- Do not justify the change of scale by calling it an artefact.

### MO2. The absolute estimand must be normalised to bin width, and 400 mg is a grid edge, not an endpoint

**Issue (bin width).** Raw Δ minutes per bin peak at 120–140 mg at every visit. This happens because the bin width doubles from 10 to 20 mg at 120 mg. Per 10 mg, the absolute peak is at 90–100 mg at all three visits:
- M3: 2.37 per 10 mg at 90–100 mg, versus 1.76 at 120–140 mg.
- K1's own example ("+3.5 min/day at 120–140 mg") is affected by this.

**Issue (upper edge).** In the bootstrap, the argmax of the relative (log ratio-of-means) profile sits in the last bins:
- M6: 375–400 mg in 95 % of replicates.
- M18: 375–400 mg in 74 %.
- M3: 300–400 mg in 88 %.

The relative profile is therefore truncated by the grid, and its upper boundary is not identified. The audit (D7) says there is no open >400 mg bin. However, `threshold_long` contains T400 (cumulative minutes above 400 mg, which equals an open bin ≥400 mg). Time below 20 mg (1440 − T20, about 1000–1040 min/day) is also available.

**Why it matters.**
- Any localisation statement on the absolute scale is a property of the grid unless it is expressed per mg.
- Any "interval" found on the relative scale will have an upper limit at 400 mg because of where the data stop, not because of biology.

**Fix.**
- Make the primary absolute estimand Δ minutes per 10 mg (equivalently, a density).
- Add the open ≥400 mg bin (= T400) and the <20 mg bin to the bin set. The bins then form a closed 24-hour composition.
- State explicitly that the upper boundary of any relative-scale region is right-censored by the measurement grid.

### MO3. The a-priori decision rule (B3) is not neutral: two of its four gates truncate the spectrum mechanically

**Issue.** Each gate behaves as follows:

- **Gate 2 (coverage ≥80 % at 0.25 min/10 mg).** Using the "either" coverage definition, this gate fails above about 240–280 mg at every visit. That happens because high-intensity time is low, not because the change is uninformative.
- **Gate 4 (|Δ| ≥ 1 min/day per bin).**
  - The gate is bin-width dependent: a 20 mg bin passes twice as easily as a 10 mg bin.
  - It fails above about 200–220 mg, again because volume falls with intensity.
  - The 1 min/day figure has no anchor (no MCID, no measurement-error estimate).
- **Gate 1 ("lower simultaneous bound >0 in at least one of two scales").** This picks the better of two scales and inflates the family-wise error rate.
- **Gate 3 (sign agreement with M-C and M-E).** This is a weak criterion. M-C (complete pairs) and M-A (LMM) share most of the data, so they will almost always agree in sign.

In my data, the simultaneous band excludes zero from about 90–100 mg up to 400 mg on both scales at every visit. Intersecting that with Gates 2 and 4 yields roughly **90–220 mg**. That is close to the 80–200 mg band the brief explicitly says must not be selected post hoc, and it would be produced by construction.

**Fix.**
- Coverage and magnitude should be **reported as separate descriptive dimensions** (as the brief asks), not used as inclusion gates.
- One primary scale for inference must be fixed in advance (see section 6).
- If an a-priori threshold is wanted:
  - define it per 10 mg;
  - tie it to a measurement-error quantity (for example, the smallest detectable change from day-level data, MO8);
  - present it as a sensitivity grid (0.25 / 0.5 / 1 min per 10 mg), not a single cut-off.

### MO4. The proposed relative-scale models do not deliver a stable estimand

**Issue.**

**log(minutes + c): c defines the estimand.** The geometric ratio at M3 260–280 mg is:
- 1.66 with c = 0.1;
- 1.29 with c = 1;
- versus a ratio of means of 1.94.

When the bin mean is 0.2–0.8 min/day, c = 1 is not a small offset.

**Gamma GLMM with an offset for zeros.** This brings the same problem back. Exact zeros are also not "≤6.5 %" in the relevant sense:
- The pooled figure of 6.5 % is correct for main visits.
- At **M0** the 375–400 mg bin has **14.3 % exact zeros**, against 0–5.8 % at follow-up.
- So zeros sit mainly at baseline, which biases any offset-based ratio.

**The mean–variance relationship differs between bins.** Within each bin, the slope of log|Δ| on log(pair mean) is:
- 0.22 at 100–110 mg (close to additive);
- 0.86 at 260–280 mg (close to constant CV).

No single variance family fits the whole spectrum.

**Fix.** The primary relative estimand should be the **ratio of marginal means** E[Y_Mx]/E[Y_M0]. Estimate it with either of these, both of which accept exact zeros and need no c:
- a log-link GEE with quasi-Poisson working variance and an independence or exchangeable working correlation, with patient-cluster sandwich SE and bootstrap; or
- a Tweedie GLMM (`glmmTMB`, `family = tweedie(link = "log")`).

Then:
- Treat model families as working models and base uncertainty on the cluster bootstrap.
- Keep log(y + c) only as a sensitivity analysis, with c tied to the bin (for example, half the smallest positive value or 1 % of the bin median) and reported across a c-grid.
- For a log link with a random intercept only, the fixed-effect ratio is also the marginal ratio. Say this in the plan, and do not add random slopes on the log scale without discussing marginalisation.

### MO5. Missingness: the specification is not MAR-adequate, and the data already show informative dropout

**Issue.**
1. **M-A is a random-intercept LMM, which assumes compound symmetry.** If the covariance structure is misspecified, MAR-based estimates of later contrasts can be biased. Four unequally spaced visits (M18 is far from M0–M6) make this structure implausible.
2. **The data contradict the plan's working hypothesis for C2** ("non-completers do not differ"):
   - Patients with M3 who later lack M18 (n = 13) had larger M3 gains than completers (n = 34).
   - Gains for non-completers vs completers: 3.2 vs 2.1 min at 90–100 mg; 2.7 vs 1.8 min at 160–180 mg; 0.85 vs 0.61 min at 260–280 mg.
   - On the log scale, their baseline high-intensity time was 0.1–0.3 SD different from completers.
   - Dropout therefore depends on observed trajectory. That is compatible with MAR, but it means the complete-pair M18 estimate (n = 39) is a selected subsample, and the LMM/complete-pair comparison is not a "sensitivity" check but an expected disagreement.
3. **Two small groups have no reason for missingness recorded:**
   - 8 patients have no M0 (11 measurements);
   - 8 patients have only M0.

   In IIM, a missing baseline (too ill, enrolled late) and a missing follow-up (death, ILD, hospitalisation) are both plausibly MNAR.
4. **No enrolment flow is described** (how many patients were eligible or enrolled versus how many had accelerometry), so selection bias into the accelerometry cohort cannot be judged.

**Fix.**
- Primary model: **MMRM-type**, with visit as a factor and an unstructured covariance (or heterogeneous Toeplitz if unstructured does not converge), for example `glmmTMB(... + us(visit + 0 | patient))` or `nlme::gls` with `corSymm` + `varIdent`. Keep the random-intercept model only as a sensitivity analysis.
- Obtain documented reasons for every missing M0 and M18. Report a CONSORT-style flow from the clinical cohort.
- Delta-adjusted (tipping-point) pattern-mixture sensitivity for M18 on the primary scale.
- Report the M3 change by later M18 status as a pre-specified diagnostic.

### MO6. "Which part of the spectrum carries change" is confounded by volume and composition, and the plan has no shape estimand

**Issue.**
- Exclusive bins plus <20 mg plus ≥400 mg sum to about 1440 min. They are compositional parts.
- At M3, total time at 20–400 mg rose by 27.7 min/day (ratio 1.07).
- Per-patient Δ of mean ENMO correlates 0.77 with Δlog(100–110 mg) and 0.65 with Δlog(260–280 mg). Δ intensity gradient correlates 0.67 with Δlog(260–280 mg). Mean IG changed by +0.048.
- So most of the per-bin signal is consistent with one smooth stretch of the intensity distribution. That makes the "unique interval" question possibly ill-posed: no localised band is needed to explain the data.
- Concentration is also uneven:
  - At 260–280 mg and 350–375 mg, the top 10 % of patients contribute 60 % and 69 % of the total gain. At 100–110 mg the figure is 51 %.
  - At high intensity the median Δ is about half the mean or less.

**Why it matters.** Per-bin marginal tests cannot distinguish:
- (a) a proportional volume increase;
- (b) a smooth shift of the distribution (an IG change);
- (c) a genuinely localised change.

Only (c) would support an "informative interval".

**Fix.** Add a pre-specified **shape estimand**:
- per-bin log-ratio change *relative to the volume-matched expectation*, for example clr- or ilr-transformed bin shares (compositional data analysis); or
- a test of the null hypothesis "all bins change by a common factor" against "shape changes".

Treat a localised interval as demonstrated only if the shape estimand shows a contiguous region that departs from a smooth monotone trend in intensity. Present C5 (IG, M-metrics) as the parsimonious competitor model, not as a tertiary question. If a single scalar such as IG or mean ENMO captures the change, report that as the honest answer.

### MO7. Multiplicity and region inference: replace "a bin passes" with confidence sets for the region

**Issue.**
- Δ correlates strongly between adjacent bins (median adjacent r = 0.95) and weakly between distant bins (60–70 vs 300–325 mg: r = 0.01). The effective number of tests is about 5–16, depending on the estimator.
- The pointwise-versus-max-t choice is fine, but the plan's family is unclear: per visit only? Across visits? Across scales?
- My max-t critical values per visit were 2.67–2.77 with B = 2000.
- What a "region" means statistically is not defined.

**Fix.**
- Define the family as 23 + 2 bins **per visit, per primary scale**. Report visits as separate pre-specified estimands, plus one joint band across visits as a sensitivity analysis.
- Infer the region as **inner and outer simultaneous confidence sets** for {intensity : effect > δ}, using the max-t band (the excursion-set logic of Sommerfeld et al.). Report "the region certainly includes X, and plausibly extends up to Y". This answers the user's question directly and replaces argmax and "selection frequency".
- Argmax selection frequency, if reported at all, belongs in a supplement with a warning that it is unstable on the absolute scale (best bin at M3 in only 25 % of replicates) and edge-dominated on the relative scale.

### MO8. Reliability is declared non-estimable, but GGIR day-level output probably makes it estimable

**Issue.** A1(iv) says ICC and SEM cannot be computed because the master has no daily data. But GGIR part 2 writes `part2_daysummary.csv`, which contains per-day time above each `MVPA_E5S_T*` threshold. That is GGIR output already on disk, not a rerun on raw `.bin` files. With 6–7 valid days per record, per-bin day-level ICC, the SEM of the weekly mean (Spearman–Brown) and a smallest detectable change can all be computed.

**Why it matters.** "Most informative for longitudinal change" is fundamentally a ratio of change to within-person noise. Without within-person noise, the group-level precision (bootstrap SE) mixes true between-person heterogeneity of change with measurement noise. This also gives the only non-arbitrary anchor for Gate 4 (MO3).

**Fix.** Ask the user for the GGIR `output_*/results/QC/part2_daysummary*.csv` files and the GGIR config (`config.csv`), and add a reliability module. If the files are unavailable, state the limitation. Do not call reliability impossible in principle.

### MO9. Measurement confounders missing from the audit: wrist side and device

**Issue.**
- Filenames show that **11 patients** were measured on the left wrist at some visits and the right wrist at others.
- **12 main-visit records** have no wrist side in the filename (some carry an undocumented `_upr` suffix).
- 60 of 72 patients wore different devices across visits (27 devices in total).
- With absolute ENMO cut-points, dominant versus non-dominant wrist placement shifts time above thresholds, and plausibly does so more at high intensity. That is precisely where M3's relative response lies.

**Fix.**
- Add a QC check (wrist side per record and change of side within patient).
- Run a sensitivity analysis that excludes records whose side differs from the patient's modal side or is unknown.
- Ask the clinic which wrist is dominant.
- Device changes are mitigated by autocalibration (cal.error.end ≤ 0.006 g). Document this but do not model it.

---

## 4. Minor objections

1. **Season (D19).** Season is better described as *structurally* linked to visit than as partially confounded with it.
   - No M0→M6 pair falls in the same season. The transitions are 12 spring→autumn, 12 autumn→spring and 9 summer→winter.
   - M6 and M18 sit about half a year from baseline, so a harmonic adjustment is identified only through variation in baseline calendar month, and M0 has only 6 winter records.
   - Report baseline-month balance, and interpret M-E at M6/M18 cautiously.
2. **Visit timing (D18).** 6 M3 records are later than 4.3 months and 14 M6 records are earlier than 6.0 months. Keep nominal visit as primary (it is the protocol estimand). M-F should use time since M0 with a natural spline, and the overlapping M3/M6 windows should be discussed.
3. **Valid days (D14).** Pre-specify a primary rule. GGIR's `includedaycrit` and weekend representation (WD/WE weighting in "AD" summaries) should come from the config. Run a sensitivity analysis with ≥5 valid days (this drops 1 + 9 records).
4. **SD standardisation in bootstrap (V3).** If SD effects are kept, use dz (SD of the change), recomputed in every replicate. The pooled SD across all records mixes in between-person variance and is not a responsiveness index. A fixed denominator understates uncertainty.
5. **Bootstrap with small samples.**
   - With 39–52 clusters and skewed bins, percentile intervals for ratios are biased. Use BCa or studentised intervals; for log ratio-of-means, the bootstrap-t on the log scale is cheap.
   - Check coverage of the max-t band by simulation: parametric resampling from the fitted MMRM, with n = 39, about 500 simulated datasets.
   - Record and report the number of failed or singular LMM fits per bin. Sparse bins with near-zero between-patient variance will be singular.
5a. **Decide how the simultaneous band is built: re-fit inside each replicate, or bootstrap only the paired means.** Then keep bootstrap SEs (studentisation) consistent with it.
6. **C1 model.** `minutes ~ visit × bin + (1|patient) + (1|patient:bin)` on raw minutes across bins whose scales differ about 400-fold is not interpretable, and an LRT with 66 interaction df at n = 71 is unreliable. Use the shape estimand (MO6), and bootstrap the difference of the M3 and M6 log-ratio profiles within patients who have both.
7. **C3 priority.** Raise C3 to priority 1. My concentration numbers (MO6) indicate that the high-intensity result is about a minority. Report the median Δ and the proportion of patients with Δ > 0 per bin; at high bins at M3 this proportion is 72–83 %.
8. **D16 (UCTD/SSc).** Lock the inclusion criteria before any model is run. Run the primary analysis with and without that patient.
9. **Coverage definition.** "Either visit" is a statement about pairs; per-visit coverage describes the population. Report both. Coverage describes how prevalent exposure is, not whether a bin is informative.
10. **Repro file.** `repro_r05_baseR_check.txt` ends with an error. Rerun it to completion, or mark the trailing part as not executed.
11. **Location citation for V1.** The audit cites R06 lines 176–177. In both R06 versions present, the RNG call is at lines 81–82 (and in R03 at lines 53–56, which is correct).
12. **Worker count.** Taking `min(myo_ncores, detectCores() - 1)` is fine. With `furrr` and `seed = SEED`, results do not depend on it. Add the 1-worker vs N-worker identity test (already in B4).

---

## 5. Verification of audit claims

| Claim | How I checked | Result |
|---|---|---|
| K3: `season_anchor_test` is empty because of `` comparison$`Chi Df` `` | Read A04 lines 179–200 and A07 line 277 (same code). Unzipped `A04_FINAL_cumulative_analysis.xlsx`: sheet `season_anchor_test` has dimension `A1:J1`, so it holds only the 10-column header. Both early-exit branches return a one-row tibble with a `note`, so 0 rows can only come from a zero-length column. lme4 is 2.0-1 according to `sessionInfo.txt`; in current lme4, `anova.merMod` labels the column `Df`. I could not load lme4 here to inspect the column names directly. | **CONFIRMED** (mechanism inferred from code and output, not from executing lme4) |
| V1: RNG depends on worker count in R03/R06 | Read R03 lines 53–56 and R06 lines 81–82 (`clusterSetRNGStream` + `parLapply`). Ran a base-R `parallel` demo with the same pattern: 2 vs 3 workers gave different replicate values from replicate 5 on; 2 vs 2 workers gave identical values. | **CONFIRMED** (line numbers for R06 are wrong, minor #11) |
| V2: R05 resamples independently per visit × bin | Read R05 `boot_ci` inside `group_modify` by (visit, descriptor) | **CONFIRMED** |
| K2: revised R00 hard-codes the `R:/` path and ignores `MYOII_MASTER_RDS` | Revised R00 line 6 vs older `R00_MYOII_common.R` line 6 (uses `Sys.getenv`) | **CONFIRMED** |
| D3/D4/D5: 63/52/57/45 records; 47/52/39 pairs; 8 patients without M0 | Base R on `record_level` | **CONFIRMED.** Also found: 8 patients with M0 only. |
| D13: zeros ≤6.5 % | Base R | **CONFIRMED pooled** (6.5 % over main visits, 6.8 % over all 251), but **misleading**: M0 at 375–400 mg = 14.3 % |
| D18: timing ranges | Base R from `start_time` | **CONFIRMED** (M3 2.0–6.0, M6 4.3–9.4, M18 16.5–21.8 months) |
| K1 numbers (pooled SD 24.9 → 0.40; 0.585 SD = 0.68 min; 0.41 SD = 3.5 min) | Base R | **Numbers CONFIRMED** |
| K1 interpretation ("high-intensity response is mostly a scale artefact") | Relative, dz and compositional profiles (MO1) | **NOT CONFIRMED** |
| D7: no open >400 mg bin available | `threshold_long` has T400 (= open ≥400 mg bin) | **NOT CONFIRMED as stated:** absent from `bin_long`, but derivable from `threshold_long` |
| A1(iv): reliability not estimable | Depends on GGIR day summaries outside the master | **UNABLE** to verify the files exist. The claim is too strong (MO8). |

---

## 6. Recommended changes to the locked analysis plan

### Estimands (pre-specified, each answering a named question)

| Code | Question | Estimand | Role |
|---|---|---|---|
| **E1** | Where in the spectrum is time gained? | Δ E[min/day per 10 mg], Mx − M0, per bin, including <20 mg and ≥400 mg | Primary absolute |
| **E2** | Where is the relative change largest? | log ratio of marginal means E[Y_Mx]/E[Y_M0] per bin | Primary relative. Upper edge reported as censored at the grid. |
| **E3** | Does the shape of the distribution change beyond volume? | Change in clr-transformed bin shares (24 h composition); test of the "common factor" null | Primary shape. This is the only estimand that can support a "specific interval" claim. |
| **E4** | Is a change at the group level detectable against day-to-day noise? | Per-bin day-level ICC, SEM, SDC (if `part2_daysummary` is available); otherwise dz as a secondary index | Secondary |

Secondary or descriptive: coverage per visit (limit grid), median Δ, proportion of patients with Δ > 0, concentration (top-10 % share), SD effect (dz) for comparison with A04/A07.

### Models

- **E1:** MMRM (visit factor, unstructured covariance), all main-visit records, MAR.
- **E2:** log-link GEE (quasi-Poisson working variance) and/or Tweedie GLMM. Both handle exact zeros.
- **Sensitivities:**
  - complete pairs;
  - random-intercept LMM;
  - log(y + c) over a c-grid;
  - season harmonics;
  - actual time since M0 (spline);
  - ≥5 valid days;
  - consistent wrist side;
  - without the UCTD/SSc patient;
  - delta-adjusted MNAR for M18.
- **Cumulative thresholds:** same estimands E1/E2 as a complementary representation. Add the identity T_x = Σ bins ≥ x as a check.

### Uncertainty

- Patient-cluster bootstrap with B ≥ 2000. One draw is shared across all bins and visits. Duplicated patients get new IDs.
- Studentised (bootstrap-t) or BCa pointwise intervals.
- Max-t simultaneous band per visit per primary estimand.
- Inner and outer confidence sets for {intensity : effect > 0} and for {effect > δ} on a δ-grid.
- Coverage of the band checked by simulation at n = 39.
- `furrr` with `seed = SEED`. Log B, seed, workers and failed or singular fits per bin.

### Decision rule (replacing B3)

- **Localised informative interval.** Claim one only if both hold:
  - the E3 inner confidence set contains a contiguous region of ≥2 adjacent bins whose shape departure from a monotone smooth trend in log-intensity is supported by the simultaneous band;
  - the region keeps its direction under the MMRM, complete-pair and wrist-side sensitivities.
- **Otherwise, report one of two descriptive answers:**
  - "change is a smooth redistribution of intensity, best summarised by [IG / mean ENMO]"; or
  - "cannot be identified uniquely".
- Use the E1/E2 confidence sets to describe where the change lies (for example, "from about 90 mg upwards; the upper limit is not identifiable beyond 400 mg").
- Coverage and magnitude gates (≥80 %, ≥1 min/day) are removed as gates. They are reported alongside the result as sensitivity grids per 10 mg.

### Sliding windows

Drop them as inference; I agree with the plan. Keep at most one supplementary figure showing that window sums are linear combinations of E1, with max-t over all windows jointly if any CI is shown. Region inference (MO7) replaces their purpose.

### Figures

1. E1 per-10 mg profile with pointwise and simultaneous bands.
2. E2 log ratio-of-means profile with the bands and a marked grid edge.
3. E3 clr share-change profile.
4. Per-visit coverage and median Δ / proportion with Δ > 0 panels.
5. Confidence-set bar (inner / outer region) per visit.
6. Flow diagram including reasons for missingness and wrist-side QC.

Do not show argmax-frequency figures in the main text.

---

## 7. Literature (only items I believe are real; DOIs marked "to verify" where I am not certain)

- Rowlands AV, Edwardson CL, Davies MJ, et al. Beyond cut points: accelerometer metrics that capture the physical activity profile. *Med Sci Sports Exerc* 2018;50:1323–1332. DOI 10.1249/MSS.0000000000001561.
- Rowlands AV, Dawkins NP, Maylor B, et al. Enhancing the value of accelerometer-assessed physical activity: meaningful visual comparisons of data-driven translational accelerometer metrics. *Sports Med Open* 2019;5:47. DOI 10.1186/s40798-019-0225-9 (to verify).
- Hildebrand M, van Hees VT, Hansen BH, Ekelund U. Age group comparability of raw accelerometer output from wrist- and hip-worn monitors. *Med Sci Sports Exerc* 2014;46:1816–1824. DOI 10.1249/MSS.0000000000000289.
- Migueles JH, Rowlands AV, Huber F, Sabia S, van Hees VT. GGIR: a research community-driven open source R package for generating physical activity and sleep outcomes from multi-day raw accelerometer data. *J Meas Phys Behav* 2019;2:188–196. DOI 10.1123/jmpb.2018-0063.
- Chastin SFM, Palarea-Albaladejo J, Dontje ML, Skelton DA. Combined effects of time spent in physical activity, sedentary behaviors and sleep on obesity and cardio-metabolic health markers: a novel compositional data analysis approach. *PLoS One* 2015;10:e0139984. DOI 10.1371/journal.pone.0139984.
- Dumuid D, Stanford TE, Martin-Fernández JA, et al. Compositional data analysis for physical activity, sedentary time and sleep research. *Stat Methods Med Res* 2018;27:3726–3738. DOI 10.1177/0962280217710835 (to verify pages).
- Sommerfeld M, Sain S, Schwartzman A. Confidence regions for spatial excursion sets from repeated random field observations, with an application to climate. *J Am Stat Assoc* 2018;113:1327–1340. DOI 10.1080/01621459.2017.1341838 (to verify).
- Westfall PH, Young SS. *Resampling-Based Multiple Testing*. Wiley, 1993 (book; max-T method).
- Field CA, Welsh AH. Bootstrapping clustered data. *J R Stat Soc B* 2007;69:369–390. DOI 10.1111/j.1467-9868.2007.00593.x (to verify).
- Benjamini Y, Yekutieli D. False discovery rate-adjusted multiple confidence intervals for selected parameters. *J Am Stat Assoc* 2005;100:71–81. DOI 10.1198/016214504000001907.
- Brooks ME, Kristensen K, van Benthem KJ, et al. glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. *R Journal* 2017;9:378–400. DOI 10.32614/RJ-2017-066.
- Mallinckrodt CH, et al. Recommendations for the primary analysis of continuous endpoints in longitudinal clinical trials. *Drug Inf J* 2008;42:303–319 (to verify).
- Dominant versus non-dominant wrist effect on ENMO / time above thresholds: a citation is needed (to verify; for example, Dieu et al., *Clin Physiol Funct Imaging* 2017, not verified).

The plan's unverified PMIDs (32920327, 33023619, 26881287) were **not checked** by me and should not be cited until they are resolved.

---

## 8. Re-review of disposition

- **Reviewer:** Claude Fable 5.1
- **Object:** `261013_phase3_disposition_and_locked_plan.md`
- **Checks:** new base-R checks on the master, paired M0→M3 and M0→M6. I fitted no MMRM.

**Verdict: REVISE.** The scope is narrow. MO1–MO3, MO7 and MO9 are adequately resolved: MO9 by the user's statement that the device is always on the non-dominant wrist, and the rest by the changes in the plan. MO8 is acceptably deferred. The new decision rule in §3.4, as written, would produce false "localised intervals", for the reasons in R1 and R2.

**R1. "Linear in log-intensity" is not an adequate smooth null. Major.**

On the data, the mean paired Δlog profile is sigmoid, not linear. At M3 (zeros replaced by 1/12 min/day) it runs:
- 20–30 mg: −0.01
- 90–100 mg: 0.13
- 200–220 mg: 0.45
- 260–400 mg: about 0.65, a plateau

Residuals from a linear fit therefore form a U shape at every zero-replacement value I tried. For example, M3 gives `+++------------++++++++` across the 23 bins, with max |r| about 0.19. So the inner set of residuals >0 will tend to flag the two **grid ends**, 20–40 mg and about 260–400 mg, as "localised intervals". That is curvature of a smooth monotone shift, not localisation.

The linear null has one point in its favour: it is exactly the intensity-gradient (IG) model, so it is still useful as a descriptive "is IG sufficient?" check. It cannot be the null for a localisation claim.

Fix:
- Use a **smooth monotone null**: a natural spline in log-midpoint with df = 3 (or a shape-constrained monotone fit), refitted inside every bootstrap replicate.
- Define "localised" as a **bump**: ≥2 contiguous bins with simultaneous residual >0 (or <0; the test must be two-sided), flanked on both sides by bins where it is not. Exclude regions that touch 20 mg or 400 mg.
- Pre-register a small power simulation, for example a +0.2 log bump over 2–3 bins at n = 39–52. If power is low, outcome (b) is the expected result and should be stated as such.
- Rephrase (a) as "no evidence of a localised departure; consistent with smooth redistribution". It must not read as proof that none exists.

**R2. Zero replacement drives E3 in the sparse tail. Major.**

| Visit | Zero replacement (min/day) | Mean Δclr, 375–400 mg |
|---|---|---|
| M6 | 0.01 | 0.55 |
| M6 | 1/12 (one 5-s epoch per day) | 0.28 |
| M6 | 0.5 | 0.06 |
| M3 | 0.01 | 0.45 |
| M3 | 0.5 | 0.27 |

M6 itself has no zeros in this bin; the swing comes entirely from the 14.3 % baseline (M0) zeros. The ≥400 mg part is stable, at 0.43–0.47 for M6 (T400 zeros are only 0–3.2 %).

Fix:
- **Amalgamate the sparse tail into one part**, ≥300 mg (= T300). That leaves 300–325 mg as the only bin with any baseline zeros.
- Pre-specify multiplicative replacement with a detection limit of one epoch per recording (5 s / N valid days), for example `zCompositions::multRepl`.
- Require the decision to be invariant across a replacement grid {0.01, 1/12, 0.5}, and add that as condition 2 of the rule.

**R3. Centring and closure are well-defined but partly irrelevant. Clarify.**

With an intercept in the trend fit, residuals are invariant to clr centring versus plain Δlog: the maximum difference I found was 4e-16. The 25-part closure and the geometric-mean centre therefore do not affect E3 residuals.

The <20 mg and ≥400 mg parts have no intensity coordinate. Keep them out of the trend fit and report them as separate parts; assigning them an x-position would be arbitrary.

Name the x-coordinate in the plan, either the log of the arithmetic midpoint or the log of the geometric mean of the bin bounds. The choice matters slightly for the 25-mg-wide bins.

**R4. The E2 primary assumes MCAR, which contradicts MO5. Major.**

The ratio of means of all available measurements (an independence GEE without weights) compares M0, which includes the 8 patients with M0 only, against selected follow-up samples. That is valid only under MCAR, and MO5 found trajectory-dependent dropout.

Fix: make one of these the primary, and use the all-available ratio and complete pairs as sensitivities:
- the Tweedie GLMM (log link, MAR), the current S3;
- an IPW-weighted GEE.

E3 and E4 are paired, which is complete-case. State that explicitly, and give at least the M18 caveat for both.

**R5. Rule condition 2 is partly redundant. Minor.**

S1 (complete pairs) is the same thing as the paired E3, so it adds nothing. Replace it with:
- invariance across zero replacements (R2);
- S7 (≥5 valid days);
- S8 (without the UCTD/SSc patient).

Bootstrap requirement: the trend or spline must be refitted in every replicate, and residuals studentised with bootstrap SEs. The residual covariance is singular (residuals are orthogonal to the design), which is fine for max-t but must not be inverted.

**R6. Implementation notes. Acceptable.**
- MMRM on y/sd(y) is scale-equivariant, so it is fine provided the same constant is used for back-scaling inside each replicate. Log the BFGS convergence and positive-definite Hessian per replicate, not only for the original fit (23/23 bins converged).
- Pre-generating the bootstrap draws in the main process is a correct fix for V1.

---

## 9. Round 3 — review of results and report

- **Reviewer:** Claude Fable 5.1
- **Date:** 2026-10-14
- **Materials:**
  - `reports/261013_thresholds_report.Rmd`, read in full;
  - the saved results: `output/RData/{qc,primary,sensitivity,descriptive,validation}_objects.RData`, loaded in base R;
  - `06_VALIDATION_BOTH.R`, plasmode section.
- **Not done:**
  - I did not re-fit any model.
  - I did not inspect the rendered HTML; Rmd line numbers below refer to the source file.

**Verdict: REVISE.**

The primary analyses are sound and I do not ask for them to be re-run. The required changes concern **interpretation and text**, plus **two additional validation runs** (items 5 and 6).

### 9.1 What checks out against the saved objects

- E1 time above 100 mg (T100), with 95 % CI:
  - M3: 20.86 (11.60–31.05)
  - M6: 21.00 (11.51–30.96)
  - M18: 16.76 (6.76–26.48)
- E1 per-bin values.
- E2 ratios: 2.20 at 375–400 mg at M3; 1.40–1.68 at M6 and 1.32–1.49 at M18 above 240 mg.
- E3 tilt: 0.28 / 0.25 / 0.22. The verdict is invariant across all five variants.
- C1 estimates and bands.
- Coverage at 260–280 mg: 55.7 % at M0 vs 74.0 % at M3.
- C5 correlations: ρ with ΔIG 0.68 / 0.69 / 0.71; ρ with Δ mean ENMO 0.43 / 0.53 / 0.62.
- Season-adjusted E1 is within 0.38 bootstrap SE of the primary (my check of S5; the claim on line 915 holds).
- Leave-one-patient-out: maximum shift 0.93 SE.

### 9.2 Required changes

1. **The issue-log note contradicts the withdrawn K1 (line 1029).**
   - The text says the M3 high-intensity claim was "partly an artefact of scaling by the standard deviation". This is the claim the disposition withdrew.
   - Replace with: "The location of the maximum depended on the estimand; the high-intensity relative response is real but small in absolute minutes and concentrated in a minority."

2. **Bin text sits under the cumulative figures (lines 456–576).**
   - After `F1c_E1_cumulative` (line 457), the first note and finding boxes describe *exclusive bins per 10 mg*. After `F2c` (line 519), the first boxes describe *bin* ratios.
   - The bin figures `F1`/`F2` (lines 581–585) have no explanatory text.
   - Move the bin boxes (lines 460–484 and 524–549) under F1/F2.

3. **The text implicitly recommends the forbidden post-hoc band (line 642).**
   - The text reads: "the mid range (~80–200 mg) … suits group-level monitoring".
   - This band is high-coverage and high-precision *by construction*, because that is where most of the time is.
   - The brief explicitly excludes post-hoc selection of 80–200 mg. Delete the line, or replace it with: "coverage and precision are highest where most time is spent; this is a property of the time distribution, not evidence for a preferred band."

4. **"IG is sufficient" is overstated (lines 259–262, 794–795).**
   - ρ ≈ 0.7 means the two measures share only about half of the between-patient ranking. Both are slopes of log time against log intensity, so a high correlation is expected by construction.
   - At group level they also **diverge in time course**:
     - mean ΔIG *increases* M3 → M6 → M18 (0.048 / 0.079 / 0.105);
     - the clr tilt *decreases* (0.277 / 0.247 / 0.218).
     - A plausible cause is that GGIR computes IG on 0–4000 mg in 25-mg bins.
   - Rephrase: "IG is a closely related, cut-point-free summary that tracks the tilt moderately (ρ ≈ 0.7); the two diverge in their time course, so IG is a reasonable but not equivalent summary." Report the divergence in the C5 box.

5. **The "null false-flag" is not a null, and "not just a lack of power" is overstated (lines 239, 1003, 1011–1014).**
   - `sim_one` with `null_no_bump` resamples the *observed* data. Those data carry their own small residuals, for example +0.03 to +0.04 at 100–140 mg at M3.
   - The 6 / 7 / 15 % figures are therefore the bootstrap stability of a flag under the observed data, not a type-I error under a true null.
   - With N_SIM = 100 and only B_SIM = 200 bootstrap replicates, the Monte-Carlo uncertainty is large: the 95 % interval for 15 % is about 9–23 %.
   - Required:
     - (a) rerun a **true-null plasmode**: subtract the observed E3 residual from each patient's paired clr change before resampling;
     - use N_SIM ≥ 500 and B_SIM ≥ 1000;
     - report the type-I rate with a Monte-Carlo CI at each visit.
     - (b) If the rate exceeds 5 % at M18, state that the max-t bands, and hence the E1/E2 **inner sets at M18** (for example E1 120–260 mg), may be anti-conservative at n = 39.
     - (c) Replace "therefore not just a lack of power" with: "a large (+22 %), narrow (2–3 bins) band would very likely have been detected at M3/M6; smaller or broader departures were not tested."

6. **The spline null with df = 3 can absorb broad bumps (lines 239, 1013; scope of the E3 conclusion).**
   - Over the bounded bins, log-intensity spans about 2.4 log units (25–290 mg), and the natural spline has 2 interior knots.
   - A departure about 0.5–0.7 log units wide (for example 70–140 mg, 4–6 bins) will be largely absorbed. The test can therefore only detect **narrow** bands.
   - The linear-null residuals are about 0.16 at 20–30 mg, about −0.08 at 80–90 mg and about 0.18 at 280–300 mg. That is exactly the curvature the spline absorbs.
   - Required:
     - add power scenarios with 4- and 6-bin bumps (for example 70–110 and 90–160 mg) and an amplitude of +0.1;
     - state the resolution explicitly: "E3 excludes narrow (≈20–60 mg) localised departures of ≥ +0.2 log; it cannot distinguish a broad mid-range emphasis from smooth curvature."
   - Also: the `detected` flag uses `grepl` on the lower and upper edge strings. A bump in an adjacent bin, such as 160–180 or 240–260 mg, can count as detected. Require overlap with the injected bins.

7. **The S9 tipping point is a mechanical identity (lines 949–951).**
   - tipping δ = E / p with p = 24/63, so "about 2.6 times the gain in every bin" is 63/24 by construction, not a finding.
   - The tipping is also defined on the point estimate, not on the inference.
   - Report the δ at which the **95 % lower bound** reaches 0. Approximately, δ ≈ lower bound / p:
     - about 1.5 min/day per 10 mg for 90–100 mg, against a point-estimate tipping of 7.2;
     - about 18 min/day for T100, against 44.
   - Interpret these smaller values. Remove "in every bin, about 2.6 times".

8. **The M18 "replacing pottering" interpretation rests on non-significant changes (lines 478–482).**
   - The M18 changes at 20–50 mg (−3.5 to −2.0) have CIs including 0. The E2 band at 20–30 mg is 0.83 upwards, so it includes 1.
   - The 20–50 mg range also includes restless sleep, since the windows are 0–24 h.
   - Delete the "replace slow pottering with purposeful walking" sentence, or move it to a hypothesis explicitly labelled as such.

9. **Activity labels and functional claims (lines 545–546, 640).**
   - Remove "brisk walking or stair climbing" and "rising from a chair". There is no validation of wrist ENMO against activity type in IIM, and these labels contradict the caution box on line 1040.
   - Line 640, "More patients *reach* higher intensities after treatment, not only those who already did":
     - coverage compares different sets of recordings (n = 61 vs 50), so it cannot say "not only those who already did" without a paired transition table;
     - "after treatment" is causal language.
   - Rephrase as "a larger share of follow-up recordings".

10. **The location of the absolute gain is partly a precision effect (lines 242, 474, 668).**
    - At M3, E1 at 40–70 mg (2.2–2.3) is as large as at 90–110 mg (2.4). The difference is precision, not size: the SE is 1.25–1.86 at 40–70 mg against 0.76 at 90–100 mg.
    - The lower edge of the inner set (about 90 mg) and the "core 90–160 mg" therefore reflect precision.
    - Say "certainly positive from about 90 mg; below that, gains of similar size cannot be distinguished from zero" instead of "gained mostly between 70 and 140 mg".

11. **The IMNM subtype observation is a post-hoc selection (lines 254–257, 832–841).**
    - 90–100 mg is one of 23 bins × 3 subtypes × 3 visits.
    - The IMNM mean change is about 0 (mean Δ 0.01–0.05 min/day per 10 mg), but the medians are −1.2 at M6 and +1.4 at M18, indicating heterogeneity rather than "no change".
    - Required:
      - replace the single-bin highlight with a pre-specified summary per subtype (T100, tilt, ΔIG), shown as full profiles;
      - drop "most clinically provocative observation";
      - support or remove the clinical claim that IMNM has "slower or less complete recovery", which needs a reference;
      - fix "13–18 patients" (line 257) to 12–18.

12. **The "early burst" hypothesis needs its confounders named (lines 245–248, 869).**
    - Add that M3 and M6 differ in season contrast: no M0→M6 pair falls in the same season, whereas M0→M3 is about a quarter-year shift.
    - They also overlap in timing.
    - So the shape difference could be seasonal or a timing effect. The existing hedge on significance is fine.

13. **The claim "numbers are read live" is false for many statements (line 226).** These figures are hand-written:
    - "2–3 min/day per 10 mg", "×1.9–2.2", "×1.3–1.7";
    - "55–79 %" (observed 38.5–83 %), "half to two-thirds" (observed 0.42–0.75);
    - "94–100 %", "6–7 %", "2.6 times", "90–160 mg", "~20 min/day".

    Either compute them inline or state that they are hand-checked ranges from the 2026-10-14 run. Correct "55–79 %".

14. **Reference 12 (line 1121).** The title given for PMID 32920327 is unverified and may not be the paper's actual title. The brief forbids invented citations. Remove it until it has been verified on PubMed.

15. **Minor wording.**
    - Line 236–237: "every higher intensity range gains proportionally more" can be read as "every range gains share". In fact, at M3, 20–120 mg *loses* clr share.
    - Line 540: at M18 the simultaneous E2 band *includes* 1 at 90–140 mg; the inner set is 140–260 mg. Fix the sentence.
    - Line 894: "do not depend on these analytic choices". Containment within the primary CI is weak evidence, because the samples overlap. Soften to "are consistent under".
    - Line 230 ("the gain lasts"): add "under MAR; M18 has informative dropout".

### 9.3 Answers to the coordinator's questions

1. **Numbers that are not supported or are overstated:** items 1, 4, 7, 10, 11, 13 and 15. No inline-computed number disagrees with the objects. The problems are in the hand-written ranges and in the interpretation.
2. **Is "no localised band" warranted?** Yes, for **narrow, large** departures at M3 and M6. It is not established for narrow small bumps, for broad bumps, or at M18 until the true-null calibration is run (items 5–6). With the scope stated, the verdict is defensible: "no evidence of a narrow localised band; change is consistent with a smooth tilt."
3. **Can the df = 3 spline absorb bumps?** Yes: departures wider than about 4 bins (item 6).
4. **Clinical interpretations:**
   - "early burst": adequately hedged, but needs the season and timing caveat;
   - IMNM: post-hoc, needs reframing;
   - "IG is sufficient": contradicted by the opposite time course;
   - "pottering": unsupported.
5. **Required before sharing with clinicians or using in a manuscript:** items 1–9 and 14. Items 10–13 and 15 are strongly recommended.

### 9.4 Round 3 closure

**Materials checked:**
- disposition §6;
- `R/05_helpers_e3fast.R` and V7 in `06_VALIDATION_BOTH.R`;
- the revised Rmd (grep of all 15 items, plus the key messages, E3, C5, S9 and power boxes);
- the new `validation_objects$power_summary` (in `scratchpad/wfull/output/RData`).

**Verdict: ACCEPT WITH MINOR CHANGES.**

**(a) V7 engine.** I found no errors.
- **True null.** Subtracting R·mean(D) from every row gives R·mean(D0) = 0 exactly, because R is idempotent. The patient covariance and the smooth tilt are kept.
- **Bump injection.** Adding `size` to clr changes without re-closing the composition is equivalent to multiplying minutes and then re-closing. Re-closure shifts every part by the same constant, and R annihilates constants.
- **Detection.** It now requires `pos %in% hit` (overlap with the injected bins), and the edge exclusion is respected.
- **Seeds.** Seeds do not collide across visits, and the scenarios share random numbers within a visit, which is a sensible design.
- **Resampling only the complete pairs of each visit is acceptable.** E3 is a complete-pair estimand. Patients without that visit contribute nothing in the primary bootstrap apart from a random pair count, so the difference is negligible.
- **Item 6 conceded.** My concern that the df = 3 spline absorbs broad bumps was overstated for +0.2. Six-bin bumps at 70–140 mg are detected in 99.6–100 % of simulations.

**(b) Items 1–15.** All are resolved in the Rmd. I found no residual instance of:
- the "SD artefact" wording, "pottering", the activity labels, or the 80–200 mg recommendation;
- "IG sufficient", the IMNM single-bin claim, "55–79 %", or "2.6 times";
- PMID 32920327.

The time-course divergence between IG and the tilt and the lower-bound tipping point are now reported.

**(c) Remaining minor points.** None requires a re-run.
1. **Broad, small departures can be flagged in the wrong place** (lines 274 and 1082).
   - For "6 bins 70–140, +0.1", any bump is flagged in 71–94 % of simulations, but the flag overlaps the injected bins in only 49–76 %.
   - So broad, small departures are partly absorbed by the spline, and residual flags can appear at the wrong location.
   - Add half a sentence to the scope statement: "small broad departures may be mislocalised".
2. **The calibration covers the ≥2-bin bump rule only** (type-I 0.6–2 %, i.e. conservative). It does not cover the per-bin coverage of the max-t bands. Line 1081 already says so.
   - Repeat this caveat next to the confidence-set table in the executive summary (`exec_sets`) for M18 (n = 39).
   - Optionally, report the single-bin family-wise rate under the true null; the matrix engine makes that cheap.
3. **"Overshoot" wording (line 921).** Keep it labelled as a hypothesis with the same season and timing caveat as Key message 4, or reword it to "early high-intensity emphasis".
