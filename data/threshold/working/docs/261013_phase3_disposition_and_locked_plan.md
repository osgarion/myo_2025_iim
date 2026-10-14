# MYO-II / Paper A — Vypořádání oponentury (fáze 3) a uzamčený analytický plán v1

- **Datum:** 2026-10-13
- **První role:** Claude Opus 5.5
- **Oponent:** Claude Fable 5.1 (`claude-fable-5-1`), viz `261013_phase3_independent_review.md`. Verdikt **REVISE**.
- **Uživatelské vstupy po oponentuře (2026-10-13):**
  - Náramek je **vždy na nedominantním zápěstí**.
  - Subtypy: 0 + 4 = DM, 2 + 3 = IMNM, 5 = ASyS. PM (1) a UCTD/SSc jsou vyřazeni jen z analýz podtypů.

## 1. Vypořádání hlavních námitek

| ID | Námitka | Vlastní ověření (base R) | Rozhodnutí |
|---|---|---|---|
| MO1 | K1 („high-intensity peak je artefakt SD“) neplatí | **Potvrzeno.** Poměr průměrů Mx/M0 roste s intenzitou: M3 1,15 (90–100 mg), 1,52 (180–200), 1,94 (260–280), 2,16 (300–325); M6 1,19 → 1,93 (375–400). | **Přijato.** K1 přeformulováno: poloha maxima závisí na estimandu. Absolutní, relativní a tvarový estimand odpovídají na různé otázky a volí se předem. Původní interpretaci H1 „artefakt“ **stahuji**. |
| MO2 | Absolutní estimand musí být na 10 mg; 400 mg je okraj mřížky | **Potvrzeno.** Δ na 10 mg má peak v 90–100 mg (M3 2,37; M6 2,96), 120–140 mg jen 1,76 / 2,15. T400 medián 0,57 min/den. | **Přijato.** E1 = Δ min/den na 10 mg. Přidány části < 20 mg (1440 − T20) a ≥ 400 mg (= T400) pro kompozici (E3) a relativní škálu. Horní mez je u relativní škály uváděna jako **cenzurovaná mřížkou**. |
| MO3 | Rozhodovací pravidlo B3 mechanicky ořezává spektrum | Logicky správné: brány coverage ≥ 80 % a ≥ 1 min/den klesají s objemem | **Přijato.** Brány zrušeny. Coverage a magnituda se reportují deskriptivně (citlivostní mřížka). Nové pravidlo je v §3.4. |
| MO4 | log(y + c) a Gamma s offsetem nedávají stabilní estimand | **Potvrzeno.** Nuly v 375–400 mg: M0 14,3 %, M3 5,8 %, M6 0 %, M18 4,4 %. | **Přijato.** E2 = poměr marginálních průměrů. Primární odhad je neparametrický (independence-GEE ekvivalent = poměr průměrů všech dostupných měření). Citlivost: Tweedie GLMM (log link, MAR), log(y + c) na mřížce c. |
| MO5 | Random intercept nestačí pro MAR; informativní drop-out | Počty potvrzeny (8 bez M0, 8 jen M0) | **Přijato.** E1 primárně **MMRM** (glmmTMB, `us(visit + 0 \| patient)`, REML). Random intercept je citlivost. Delta-adjusted tipping point pro M18. Diagnostika „změna v M3 podle pozdější dostupnosti M18“. **Důvody chybění chybí v masteru** → otevřený požadavek na uživatele. |
| MO6 | Chybí tvarový estimand; volume vs. shape | Neověřováno nad rámec oponenta | **Přijato.** E3 = změna clr podílů 25dílné 24h kompozice (< 20, 23 binů, ≥ 400). Test nulové hypotézy „společný faktor“ a odchylky od hladkého lineárního trendu v log-intenzitě. IG / mean ENMO je úsporný konkurenční model. |
| MO7 | Inference oblasti: confidence sets místo „bin projde“ | — | **Přijato.** Max-t pás per visit per estimand. Inner / outer confidence sets pro {efekt > δ} na δ-mřížce. Argmax frekvence se nereportuje v hlavním textu. |
| MO8 | Reliabilita je odhadnutelná z GGIR `part2_daysummary` | Soubory nejsou v masteru ani ve složce | **Přijato podmíněně.** E4 = dz (Δ / SD(Δ), přepočítané v každé replikaci) jako sekundární index. Modul reliability (ICC/SEM/SDC) se připraví, **až uživatel dodá `part2_daysummary*.csv` a `config.csv`**. Formulace A1(iv) opravena na „není odhadnutelné z masteru“. |
| MO9 | Strana zápěstí, zařízení | Strany ve jménech souborů: L 188, R 16, neuvedeno 13; 10 pacientů má L i R | **Vyřešeno uživatelem:** vždy nedominantní zápěstí. Citlivostní analýza se nedělá. Nesoulad L/R značek je QC informace (pravděpodobně leváci nebo chyba v názvu souboru); D20 v QC. |

## 2. Vypořádání drobných námitek

| # | Rozhodnutí |
|---|---|
| 1 Sezóna | Přijato. „Strukturálně svázaná s návštěvou.“ Report balance baseline měsíce. Harmonická adjustace je citlivost s opatrnou interpretací M6/M18. |
| 2 Časování | Přijato. Primárně nominální návštěva. Citlivost: LMM s `ns(months_from_m0, df = 3)`, predikované rozdíly v 3 / 6 / 18 měsících. |
| 3 Validní dny | Přijato. Primárně všechny záznamy (GGIR kritéria už aplikována). Citlivost ≥ 5 dnů. Config GGIR vyžádán. |
| 4 SD v bootstrapu | Přijato. dz se přepočítává v každé replikaci. Pooled-SD efekt jen pro srovnání s A04/A07 (supplement). |
| 5 / 5a Bootstrap | Přijato částečně. E2 percentile na log škále (transformation-respecting). E1 percentile. Max-t pás ze studentizace bootstrap SE. MMRM se **refituje v každé replikaci**. Simulační kontrola pokrytí pásu (n = 39) je zařazena jako volitelný validační skript (07) s menším počtem simulací. Selhání a singularity se logují per bin. |
| 6 C1 model | Přijato. Nahrazeno srovnáním clr profilů M3 vs M6 u pacientů s oběma návštěvami (bootstrap rozdílu). |
| 7 C3 priorita | Přijato. **Priorita 1**: medián Δ, % pacientů s Δ > 0, podíl top-10 % pacientů na celkovém přírůstku. |
| 8 UCTD/SSc | Přijato. Hlavní analýza všech 71 pacientů (dle QC eligibility). Citlivost bez pacienta UCTD/SSc. |
| 9 Coverage | Přijato. Per-visit coverage (populační) primárně, „either“ (párová) sekundárně. |
| 10 Repro soubor | Přijato. Chyba byla v posledním bloku (R06 merge). Opravený běh reprodukoval všech 9 R06 peaků; doplněno do `repro_r05_baseR_check.txt`. |
| 11 Čísla řádků | Přijato. V1 v R06 je na ř. 81–82 souboru R06 (moje číslování bylo z konkatenovaného výpisu). |
| 12 Workery | Shoda. Navíc jsou draws generovány předem v hlavním procesu: test 1 vs 3 workery dává identické výsledky (ověřeno). |

## 3. Uzamčený analytický plán v1

### 3.1 Populace
- QC-eligible měření M0–M18: 217 měření, 71 pacientů.
- M30 / M48 / M60 jen deskriptivně.
- Analýza podtypů (deskriptivní, supplement): DM / IMNM / ASyS, bez PM a UCTD/SSc.

### 3.2 Estimandy

| Kód | Otázka | Estimand | Škála | Role |
|---|---|---|---|---|
| E1 | Kde ve spektru přibývá čas? | Δ E[min/den na 10 mg] Mx − M0, per bin (23 binů) | absolutní | primární absolutní |
| E2 | Kde je relativní změna největší? | log(E[Y_Mx] / E[Y_M0]), per část (23 binů + ≥ 400) | relativní | primární relativní; horní okraj cenzurovaný mřížkou |
| E3 | Mění se tvar rozložení nad rámec objemu? | Δ clr podílů 25dílné 24h kompozice (párově); centrovaný kontrast „bin vs. průměr aktivních částí“; reziduum od lineárního trendu v log-intenzitě | kompoziční | **jediný estimand, který může podpořit tvrzení o lokalizovaném intervalu** |
| E4 | Je změna velká proti variabilitě změny? | dz = mean Δ / SD(Δ), párově | bezrozměrná | sekundární |
| CUM | Komplementární kumulativní reprezentace | E1 (Δ min/den nad T) a E2 pro 24 prahů | absolutní / relativní | komplementární |

**Deskriptivně:**
- coverage per visit na mřížce limitů 0 / 0,1 / 0,25 / 0,5 / 1 min/den/10 mg;
- medián Δ;
- % pacientů s Δ > 0;
- podíl top-10 % pacientů na přírůstku;
- pooled-SD efekt pro srovnání s A04/A07.

### 3.3 Modely a nejistota

**Primární modely:**
- **E1:** MMRM per bin. `y/sd(y) ~ visit + us(visit + 0 | patient)`, glmmTMB, REML, optimizer BFGS, zpětné škálování. Fallback: random intercept (logováno).
- **E2:** poměr průměrů všech dostupných měření (independence GEE).
- **E3 / E4:** párové (M0 + Mx).

**Citlivosti:**

| Kód | Citlivost |
|---|---|
| S1 | complete pairs (E1, E2) |
| S2 | random-intercept LMM (E1) |
| S3 | Tweedie GLMM (E2) |
| S4 | log(y + c), c ∈ {0,01; 0,1; 1} (E2) |
| S5 | sezóna (faktor) a harmonické (E1) + opravený anchor test |
| S6 | čas od M0, `ns(df = 3)` (E1) |
| S7 | ≥ 5 validních dnů (E1, E2) |
| S8 | bez UCTD/SSc (E1, E2) |
| S9 | delta-adjusted tipping point pro M18 (E1) |
| S10 | změna v M3 podle pozdější dostupnosti M18 (diagnostika) |

**Nejistota:**
- Cluster bootstrap přes pacienty, B = 2000, seed 20261013.
- Draws se generují předem, jeden draw pro všechny biny, návštěvy a obě reprezentace. Duplicitní pacienti dostávají nová ID.
- Pointwise percentile CI. Max-t simultánní pás per visit per estimand (rodina = biny jedné návštěvy).
- Inner / outer confidence sets pro {efekt > δ}:
  - E1: δ ∈ {0; 0,25; 0,5; 1} min/den/10 mg;
  - E2: δ ∈ {0; log 1,1; log 1,25}.

### 3.4 Rozhodovací pravidlo (předem)

**Lokalizovaný informativní interval** je deklarován jen tehdy, když platí obě podmínky:
1. Pro E3 (reziduum od lineárního trendu clr změn v log-intenzitě) obsahuje inner confidence set (max-t pás > 0) **souvislou oblast ≥ 2 sousedních binů**.
2. Oblast si zachová směr v S1 a S8.

**Jinak** je závěr jeden z dvou:
- (a) „změna je hladké přerozdělení intenzit (posun / sklon profilu), shrnutelný IG / mean ENMO“, pokud je lineární trend E3 nenulový a rezidua nesignifikantní;
- (b) „nelze jednoznačně identifikovat“.

**E1 / E2 confidence sets** popisují, *kde* změna leží, například „od ~90 mg výše; horní mez za 400 mg neidentifikovatelná“.

### 3.5 Sliding windows
Vypuštěny jako inference. Nahrazuje je inference oblasti (§3.3).

### 3.6 Otevřené požadavky na uživatele
1. Důvody chybění M0 a M18; CONSORT-like flow z klinické kohorty.
2. GGIR `part2_daysummary*.csv` + `config.csv` (modul reliability E4b; validní dny).
3. Kontrola značek L/R ve jménech souborů u 10 pacientů (QC informace).

## 4. Skripty (fáze 4)

| Skript | Rozsah | Obsah |
|---|---|---|
| `01_QC_data_dictionary_BOTH.R` | BOTH | slovník, QC D1–D20, flow, deskriptivní profily, coverage |
| `02_PRIMARY_bootstrap_BOTH.R` | BOTH | E1 (MMRM), E2, E3, E4, CUM; jeden bootstrap; pásy, confidence sets, rozhodovací pravidlo |
| `03_SENSITIVITY_BOTH.R` | BOTH | S1–S10 |
| `04_DESCRIPTIVE_heterogeneity_BINS.R` | BINS | C3 (koncentrace), C1 (M3 vs M6 tvar), podtypy |
| `05_FIGURES_TABLES_BOTH.R` | BOTH | publikační figury a tabulky, figure legends |
| `06_VALIDATION_BOTH.R` | BOTH | jednotkové / reprodukční testy, issue log |
| `reports/261013_thresholds_report.Rmd` | — | report dle `RULE_rmd_report_style.md` |

---

## 5. Re-review (Fable 5.1, kolo 2: REVISE) — vypořádání a plán v2 (uzamčen)

| ID | Námitka | Rozhodnutí → plán v2 |
|---|---|---|
| R1 | Lineární trend v log-intenzitě je špatná nulová hypotéza (rezidua ve tvaru U → falešné „intervaly“ na okrajích). Ověřeno i v mém pozorovaném běhu (M3 rezidua + na 20–40 a 260–400 mg). | **Přijato.** Nulová hypotéza = natural spline `ns(df = 3)` v log geometrického průměru mezí binu, refit v každé replikaci. „Lokalizovaná odchylka“ = **bump**: ≥ 2 sousední biny se simultánním pásem reziduí > 0 (nebo < 0, oboustranně), z obou stran ohraničené biny bez příznaku, nedotýká se prvního ani posledního binu. Lineární trend se reportuje jen jako deskriptivní kontrola „stačí IG?“. Varianta (a): „bez důkazu lokalizované odchylky; konzistentní s hladkým přerozdělením“. **Simulace síly** (bump +0,2 log přes 2–3 biny) v `06_VALIDATION`. |
| R2 | Náhrada nul řídí E3 v řídkém konci | **Přijato.** E3 používá 21dílnou kompozici: < 20 mg, 19 binů 20–300 mg, ≥ 300 mg (= T300). Multiplikativní náhrada, detekční limit = jedna 5s epocha na záznam (5 s / N validních dnů), faktor 0,65. Rozhodnutí musí být invariantní na mřížce {epocha; 0,01; 0,5 min/den} (`03_SENSITIVITY`). |
| R3 | Centrování / uzávěr neovlivňují rezidua | **Přijato.** x-souřadnice = log(√(dolní × horní mez)). Části < 20 a ≥ 300 mg jsou mimo trend, reportují se zvlášť. |
| R4 | E2 primárně (all-available) předpokládá MCAR | **Přijato s jinou implementací.** E2 primárně = poměr marginálních průměrů z **téhož unstructured MMRM jako E1**: log((β₀ + βₓ)/β₀). Unstructured MMRM (GLS) dává konzistentní marginální průměry při MAR i bez normality reziduí. Je výpočetně dostupný v bootstrapu (Tweedie GLMM ~2 s/fit a u binu 20–30 nekonvergoval). Tweedie (S3) a all-available (E2_avail) jsou citlivosti. E3 a E4 jsou **complete-case** (výslovně uvedeno, s výhradou pro M18). Pro ≥ 400 mg se E2 bere z MMRM fitu T400. |
| R5 | Podmínka 2 pravidla je redundantní | **Přijato.** Podmínka 2 = invariance napříč náhradou nul + S7 (≥ 5 dnů) + S8 (bez UCTD/SSc). Rezidua jsou studentizována bootstrap SE, kovariance se neinvertuje. |
| R6 | Implementace | **Přijato.** Škálovací konstanta se počítá v každé replikaci a zpětně aplikuje. Metoda fitu (MMRM_us / LMM_ri / failed) se loguje pro každou replikaci (`fit_log`). |

**Plán v2 je uzamčen 2026-10-13.** Odlišnosti od v1 jsou jen ty uvedené v této tabulce. §3.2–3.4 platí s těmito změnami.

**Pozorovaný (ne-bootstrapový) kontrolní běh:** E1 a E2 konvergují jako MMRM_us u všech 23 binů a 24 prahů. Poměr podle MAR je u vysokých binů vyšší než všechna dostupná měření, např. M3 260–280 mg: 1,91 vs 1,73. To odpovídá MO5 (baseline zahrnuje pacienty jen s M0).

---

## 6. Round 3 (Fable 5.1, review of results and report: REVISE) — disposition (2026-10-14)

The primary analyses were not re-run (no reviewer request). Changes made:

| # | Item | Resolution |
|---|---|---|
| 1 | Issue-log note contradicted the withdrawn K1 | Rewritten: "location of the maximum depended on the scale; the high-intensity relative response is real, small in minutes, concentrated in a minority". |
| 2 | Bin boxes under the cumulative figures | Cause: chunk search `fig_f1` matched `fig_f1c`. The section was rewritten, and the boxes now sit under F1/F2. |
| 3 | Implicit recommendation of the 80–200 mg band | Removed. Replaced with "coverage/precision are highest where time is spent — a property of the distribution, not a preferred band". |
| 4 | "IG is sufficient" | Reworded to "closely related, not equivalent". Added the divergent group time course (ΔIG ↑, clr tilt ↓) and the GGIR IG range (0–4000 mg). |
| 5 | The "null" plasmode was not a null; Monte-Carlo noise | **V7 rewritten** with a matrix E3 engine (`R/05_helpers_e3fast.R`, V7a = identity with `est_E3()`, max difference 2.8e-16). True null = the observed mean residual subtracted per patient. N_SIM 500 × B_SIM 1000, binomial MC 95% CIs. **Result:** type-I M3 2.0% (1–4), M6 0.6% (0–2), M18 0.8% (0–2). The rule is not anti-conservative. |
| 6 | Spline null and broad bumps; `grepl` detection | Added scenarios: 2/3/4/6 bins × +0.1/+0.2. Detection now requires overlap with the injected bins. **Result:** +0.2 detected in 87–100% for 2–6 bins; +0.1 detected in 35–98% (lowest at M18 and for 6 bins). The conclusion is scoped to "≥ ~+22% excluded; ~+10% not excluded". |
| 7 | Tipping point E/p is an identity | Added δ at which the lower 95% bound reaches 0 (90–100 mg ≈ 1.5 min/day/10 mg; T100 ≈ 17.7 min/day). The identity is stated as such. |
| 8 | "Pottering" | Removed. Added a caution that 20–50 mg in 0–24 h includes restless sleep. |
| 9 | Activity labels, causal coverage wording | Removed. Coverage is now described as "a larger share of follow-up recordings", with the note that the recordings are not paired. |
| 10 | The core of the absolute gain is partly precision | Stated explicitly (SE 40–70 mg vs 90–100 mg). |
| 11 | IMNM single-bin selection | Replaced with pre-specified per-subtype summaries (ΔT100 mean/median/IQR/% gain, clr tilt, ΔIG; `04`) and full profiles (F7). The "most provocative" wording and the unsupported clinical claim were removed. n is now given as 12–18. |
| 12 | "Early burst" confounders | Added the timing overlap and season contrast. |
| 13 | "Numbers are read live" | Ranges are now computed inline (`c3_range`, `pow_range`, `null_txt`, `tip_lo`). The note distinguishes computed numbers from checked wording. "55–79%" was replaced by the computed range. |
| 14 | Reference 12 (unverified PMID title) | Removed. |
| 15 | Minor wording | Applied: "consistent under" instead of "do not depend on"; the M18 band at E2 now uses the inline set; the "gain lasts" claim carries the MAR caveat; the low-bin share loss is stated. |

**Status:** all required items (1–9, 14) and the recommended items (10–13, 15) are addressed. The user re-runs `04`, `06`, `05` and the report on Windows.

**Round 3 closure (Fable 5.1): ACCEPT WITH MINOR CHANGES.** All three minor points are applied in the Rmd:
1. a mislocalisation caveat for small broad departures (key message 2 and the power section);
2. a calibration caveat (bump rule ≠ per-bin band coverage) under the executive-summary confidence-set table, highlighting M18;
3. "overshoot" reworded to "early high-intensity emphasis", with the timing and season caveat.

The plan v2 analysis is final, pending the user's Windows re-run of `04`, `06`, `05` and the report.
