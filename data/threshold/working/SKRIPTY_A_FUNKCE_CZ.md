# MYO-II / Paper A — přehled skriptů a funkcí

- **Autor:** Veronika Balounová (připraveno s Claude Opus 5.5)
- **Datum:** 2026-10-20
- **Stav kódu:** větev `thresholds`, commit `ccc46d5` (původní pipeline 01–06 + revize A07 07–10)

Dokument doplňuje `README.md`. README popisuje požadavky, umístění dat a spuštění. Tady je **co který skript dělá, na čem závisí, co vytváří** a **katalog všech funkcí**. Podrobnosti metod jsou v hlavičkách skriptů (PURPOSE / INPUTS / ANALYSES / OUTPUTS) a v komentářích nad funkcemi.

---

## 1. Architektura

```
R/00_config.R ─ cesty, balíčky, master .rds + md5, konstanty, bezpečné ukládání
      │  (každý skript ho načte jako první; ten pak načte R/01–R/06)
      ▼
01 QC ──► 02 PRIMARY ──► 03 SENSITIVITY ──┐
   │          │                            ├──► 05 FIGURES ──► report 261013 (v1)
   │          ├──► 06 VALIDATION ──────────┤
   └──► 04 DESCRIPTIVE ───────────────────┘
                 │
Revize A07:      ▼
   07 REV A/B  (qc, descriptive)      ─┐
   08 REV C/D  (primary, descriptive) ─┼──► 10 REV VALIDATION ──► reporty 261016_v2 (EN) a v2_cz (CZ)
   09 REV E    (+ klinický .xlsx)     ─┘
```

- **Výpočet a report jsou oddělené.** Skripty ukládají výsledky do `output/RData/*_objects.RData`. Reporty jen načítají hotové objekty a nic nepočítají.
- **Původní skripty 01–06 revize nezměnila.** Revize A07 je v nových souborech 07–10 a `R/06_helpers_revision.R` a jen čte původní objekty.
- **Výstupy bez dat pacientů.** Každé `RData` se ukládá přes `save_public()`. Ta odstraní řádky na úrovni pacienta, i ty skryté v prostředích ggplot objektů. Pokud nějaké zůstanou, odmítne soubor zapsat.
- **Shodné výsledky bez ohledu na počet jader.** Výběry pacientů pro bootstrap se generují předem v hlavním procesu (seed `20261013`). Počet workerů proto výsledek nemění.

---

## 2. Skripty

Pořadí spuštění: `run_all.R` (01 → 10 + tři reporty). Doby běhu odpovídají Windows s 31 workery a B = 2000.

| Skript | Rozsah | Účel | Vstupy (kromě masteru) | Hlavní výstupy | Doba |
|---|---|---|---|---|---|
| `01_QC_data_dictionary_BOTH.R` | biny + kumulativní | Datový slovník, kontroly integrity D1–D19 (PASS/WARN/FAIL), flow kohorty, deskriptivní profily a pokrytí. **Při FAIL se zastaví.** | — | `qc_objects.RData`, `*_01_QC_data_dictionary.xlsx`, figury `01_descriptive_profile_*` | < 1 min |
| `02_PRIMARY_bootstrap_BOTH.R` | obojí | Primární analýza (uzamčený plán v2): E1 (absolutní změna, MMRM), E2 (poměr marginálních průměrů), E3 (clr změna tvaru vs hladký null), E4 (dz). Jeden pacientský cluster bootstrap pro vše, bodové a simultánní (max-t) intervaly, inner/outer množiny, rozhodnutí o lokalizaci. | `qc_objects` | `primary_objects.RData`, `primary_boot_replicates.rds`, `*_02_primary_results.xlsx` | 1,5–2 h |
| `03_SENSITIVITY_BOTH.R` | obojí | Citlivostní analýzy S1–S10 (kompletní páry, RI-LMM, Tweedie, log(y+c), sezóna, čas od M0, ≥ 5 validních dní, bez UCTD/SSc, tipping point M18, dropout M3) a invariance rozhodnutí E3. | `primary_objects` | `sensitivity_objects.RData`, `*_03_sensitivity_results.xlsx` | 10–15 min |
| `04_DESCRIPTIVE_heterogeneity_BINS.R` | biny | C3 (stabilní posun vs menšina pacientů), C1 (tvar změny M6 vs M3), C5 (tilt vs ΔIG / mean ENMO), deskriptivně podtypy a M30–M60. | — | `descriptive_objects.RData`, `*_04_descriptive_heterogeneity.xlsx` | ~2 min |
| `05_FIGURES_TABLES_BOTH.R` | obojí | Publikační figury F1–F7 (ČB tisk) a legendy. Nic nefituje. | `qc`, `primary`, (`sensitivity`, `descriptive`) | `figure_objects.RData`, figury `05_F*`, `*_05_figure_legends.md` | < 1 min |
| `06_VALIDATION_BOTH.R` | obojí | Validace V1–V8: identity bin/kumulativní, reprodukovatelnost (1 vs N workerů), konvergence modelů, leave-one-out, kalibrace a síla pravidla E3 (plasmode), issue log. | `primary`, (`sensitivity`) | `validation_objects.RData`, `*_06_validation.xlsx` | 5–10 min |
| `07_REV_A_B_profile_concentration_BINS.R` | biny | **A07-A:** deskriptivní profil normalizovaný na 10 mg šířky binu. **A07-B:** koncentrace kladných zisků v top 10 % pacientů, net/gross, původní estimand s příznakem nestability, bootstrap CI. | `qc`, `descriptive` | `rev_AB_objects.RData`, `*_07_REV_A_B_*.xlsx`, figury `07_A1`, `A2`, `B1`, `B2` | ~13 min |
| `08_REV_C_D_E3_shape_BOTH.R` | obojí | **A07-C:** konstrukce a invariance E3, varianty hladkého nullu (df 1–5, váhy šířkou) a mřížky, plasmode volby df. **A07-D:** párové srovnání tvaru změny M3/M6/M18 (biny, sklon, globální test), Holm, společný vzorek, načasování a sezóna. | `primary`, `descriptive` | `rev_CD_objects.RData`, `*_08_REV_C_D_*.xlsx`, figury `08_C1`, `C2`, `C3`, `D1` | ~1 min |
| `09_REV_E_subtypes_clinical.R` | souhrny pacientů | **A07-E (exploratorní):** podtypy DM / IMNM / ASyS s adjustací na výchozí ph_VAS (primárně) nebo MMT8. Modely E_a–E_e, LRT, bootstrap kontrastů. | `klinicka_data_12102026.xlsx` | `rev_E_objects.RData`, `*_09_REV_E_*.xlsx`, figura `09_E1_F7b` | ~17 min |
| `10_REV_VALIDATION_BOTH.R` | obojí | Validace revize R1–R10: nezávislé přepočty v base R, reprodukce původních čísel, bootstrap pásy, velikosti vzorků, LOO, shoda figur, tabulek a textu. | vše výše + klinický `.xlsx` | `rev_validation_objects.RData`, `*_10_REV_validation.xlsx` | ~15 min |
| `run_all.R` | — | Najde `working/`, spustí 01 → 10 v čistých prostředích a vyrenderuje tři reporty. | — | vše | ~3 h |

Označení v názvu: **BINS** = exkluzivní biny, **CUMULATIVE** = kumulativní prahy, **BOTH** = obojí.

**Reporty** (`reports/`) jen načítají `output/RData/*_objects.RData`:

| Report | Načítá | Jazyk |
|---|---|---|
| `261013_thresholds_report.Rmd` | objekty 01–06 | EN |
| `261016_thresholds_report_v2.Rmd` | objekty 01–10 | EN |
| `261016_thresholds_report_v2_cz.Rmd` | objekty 01–10 + vrstva překladu popisků figur a tabulek (`CZ`, `cz_plot()`, `cz_fix()`) | CZ |

**Přepínače přes proměnné prostředí** (nastavují se před spuštěním):
- `MYOII_MASTER_RDS`, `MYOII_CLINICAL_XLSX`, `MYOII_WORKING_DIR` — cesty;
- `MYOII_B` — počet bootstrap replikací (výchozí 2000);
- `MYOII_NSIM`, `MYOII_BSIM` — velikost plasmode simulací (výchozí 500 / 1000).

---

## 3. Katalog funkcí

Sloupec „Volá“ uvádí skripty mimo vlastní soubor. Pomlčka znamená interní pomocnou funkci. Lokální funkce skriptů jsou v §3.8.

### 3.1 `R/00_config.R` — konfigurace a bezpečné ukládání

Načítá se jako první. Najde `working/`, zkontroluje balíčky, najde master a spočítá jeho md5, nastaví konstanty a výstupní složky a načte `R/01–R/06`.

| Funkce | Co dělá | Volá |
|---|---|---|
| `locate_working_dir()` | Najde složku `working/`: `MYOII_WORKING_DIR` → aktivní dokument RStudia → pracovní adresář. | — |
| `check_packages(pkgs)` | Ověří balíčky. Při chybějícím se zastaví a vypíše příkaz `install.packages(...)`, sám nic neinstaluje. | 09 |
| `locate_master()` | Najde master `.rds`: `MYOII_MASTER_RDS` → `working/data/` → nadřazená složka. Vrací cestu, z ní se počítá `INPUT_MD5`. | — |
| `date_prefix()` | Datum jako `YYMMDD` pro názvy výstupů. | 05, R01 |
| `holds_patient_rows(v)` | Pozná, zda objekt (i vnořený) obsahuje tabulku se sloupcem `patient_id`. | — |
| `walk_envs(x, on_env)` | Rekurzivně projde všechna prostředí dosažitelná z objektu (ggplot, quosures, closures). | — |
| `scrub_patient_envs(x)` | Na soukromé kopii objektu smaže z prostředí vazby s řádky pacientů. | — |
| `check_no_patient_envs(x)` | Vrátí názvy vazeb s řádky pacientů, které v objektu zůstaly. | — |
| `save_public(obj_name, file)` | Uloží objekt do `RData` bez řádků pacientů. Pokud nějaké zůstanou, zápis odmítne. | 01–10 |

Konstanty: `MAIN_VISITS` (M0, M3, M6, M18), `FOLLOWUP_VISITS`, `SUBTYPE_LABELS` (0 + 4 = DM, 2 + 3 = IMNM, 5 = ASyS), `COVERAGE_LIMITS`, `B_BOOT`, `SEED`, `N_WORKERS`, `PATHS`, `PKG_ANALYSIS`, `PKG_REPORT`.

### 3.2 `R/01_helpers_data.R` — data, QC, výstupy

| Funkce | Co dělá | Volá |
|---|---|---|
| `load_master(path)` | Načte master a ověří komponenty (`record_level`, `threshold_long`, `bin_long`, `qc`, `metadata`). | 01–04, 06–10 |
| `record_info(acc)` | Kanonické kovariáty záznamu: návštěva, způsobilost, měsíc a sezóna, počet validních dní, mean ENMO, podtyp, pohlaví, měsíce od M0. Jediný zdroj pro QC i modely. | 01, 08, 09, 10 |
| `qc_row(...)` | Jeden řádek QC tabulky (id, kontrola, stav, detail). | — |
| `qc_master(acc)` | Uzamčené kontroly integrity D1–D19 (mřížka, monotonicita, identity binů, nuly, validní dny, časování, značky zápěstí…) → PASS/WARN/FAIL. | 01, 06, 10 |
| `prepare_spectrum(acc, kind, visits)` | Dlouhá tabulka spektra. `bins` = min/den v exkluzivním binu, `cumulative` = min/den nad prahem T. Obsahuje šířku binu, intenzitu a hustotu na 10 mg. | 01, 04, 07, R03 |
| `make_pairs(d)` | Páry follow-up × M0 na pacienta a deskriptor, `delta = value − baseline`. | 04, 07, 10 |
| `theme_myo()` | Projektové ggplot téma (`theme_bw` + nastavení dle AGENTS.md). | 01, 05, 07–09 |
| `scale_visit_linetype()` | Typy čar návštěv (M0 čerchovaná, M3 tečkovaná, M6 přerušovaná, M18 plná). | 01, 05, 07 |
| `x_axis_mg(lim)` | Osa x v mg (20–400, kroky 25). | 01, 05, 07, 08 |
| `save_figure(p, name)` | Uloží figuru jako PNG (300 dpi) a TIFF (600 dpi, LZW) s prefixem data. | 01, 05, 07–09 |
| `save_tables(lst, name)` | Uloží seznam tabulek do jednoho `.xlsx` (list = sheet). | 01–04, 06–10 |
| `write_run_meta(script, extra)` | Log běhu: vstup a md5, verze masteru, B, seed, workery, `sessionInfo()`. | 01–04, 06–10 |
| `acc_meta_string()` | Popis verze masteru z `acc$metadata`. | 01, 02 |

### 3.3 `R/02_helpers_bootstrap.R` — pacientský cluster bootstrap

| Funkce | Co dělá | Volá |
|---|---|---|
| `make_draws(patients, B, seed)` | Předem vygeneruje matici B × n indexů pacientů (izolovaný seed; nezávisí na počtu workerů). | — |
| `resample_patients(d, patients, idx)` | Sestaví jednu replikaci. Každá kopie pacienta dostane vlastní ID (`<id>__b<k>`), takže náhodné efekty nejsou sdílené. | — |
| `run_cluster_bootstrap(d, stat_fun, B, seed, workers, …, checkpoint_dir, tag, chunk_size)` | Paralelní bootstrap (`furrr`) po blocích s checkpointy. Přerušený běh naváže. Chyby replikací zachytí. | 02, 03, 04, 06, 07, 09, 10 |
| `digest_key(x)` | Krátký stabilní hash (base R) pro názvy checkpointů a otisky množin pacientů. | 10 |
| `summarise_bootstrap(observed, boot, keys, band_by, level)` | Percentilové CI, bootstrap SE a simultánní max-t pás v rámci skupiny `band_by`. | 02, 03, 04, 07, 09 |

### 3.4 `R/03_helpers_estimands.R` — estimandy E1–E4 a rozhodovací pravidlo

Konstanta `ZERO_FRAC = 0,65` (multiplikativní náhrada nul).

| Funkce | Co dělá | Volá |
|---|---|---|
| `build_analysis_long(acc, visits)` | Společná analytická tabulka s bloky `bins`, `cumulative` a `composition` (biny + ≥ 400 mg) a kovariátami. | 02–04, 06–10, R04, R06 |
| `descriptor_levels(dl, block)` | Pořadí deskriptorů (binů / prahů) daného bloku. | 02, 04, 07 |
| `fit_mmrm_visit(x, response, structure)` | MMRM (glmmTMB, `us(visit + 0 \| patient)`, BFGS) na y/sd(y). Při selhání fallback na random intercept, metoda se loguje. | — |
| `fit_tweedie_visit(x, response)` | Tweedie GLMM (log link, random intercept): exp(koeficient) = poměr průměrů. | — |
| `per_descriptor(d, fitter, …)` | Spustí fitter pro každý deskriptor. Selhání vrátí jako NA s `method = "failed"`. | — |
| `est_E1_E2(dl, block, structure)` | E1 = absolutní změna (min/den na 10 mg u binů; min/den nad T u prahů) a E2 = log poměru marginálních průměrů z téhož MMRM. | 03, 06 |
| `est_E2_tweedie(dl, block)` | Citlivost S3: E2 z Tweedie GLMM. | 03, R04 |
| `est_E2_avail(dl, block)` | Citlivost: poměr průměrů všech dostupných měření (MCAR). | — |
| `e3_composition(dl, zero, frac)` | 21dílná 24h kompozice (< 20, 19 binů 20–300, ≥ 300), náhrada nul (epocha nebo pevná hodnota), clr. | 02, 04, 06, 08, 10, R05, R06 |
| `est_E3(dl, zero, frac)` | E3: průměrná párová clr změna, rezidua od hladkého nullu (natural spline df = 3 v log intenzity) a lineární sklon (tilt). | 03, 06, R05 |
| `est_E4(dl)` | E4: dz = průměr / SD párové změny (biny). | — |
| `estimate_all(dl, structure)` | Celá statistika (E1–E4 všech bloků) pro pozorovaná data i každou replikaci. | 02, 06, 10 |
| `confidence_sets(summ, deltas)` | Inner množina (dolní simultánní mez > δ) a outer množina (horní mez > δ). | 02 |
| `find_bumps(flag, labels, min_len)` | „Hrb“ = ≥ 2 sousední označené biny, ohraničené neoznačenými, mimo okraje mřížky. | R05 |
| `decide_localisation(e3_resid, e3_slope, bin_levels)` | Uzamčené rozhodovací pravidlo (plán v2 §3.4): lokalizovaný hrb vs hladká redistribuce. | 02, 03 |

### 3.5 `R/04_helpers_sensitivity.R` — citlivostní modely S1–S10

| Funkce | Co dělá | Volá |
|---|---|---|
| `wald(est, se)` | Waldovy 95% intervaly. | — |
| `sens_complete_pairs(dl, block)` | S1: E1/E2 jen z kompletních párů. | 03, 06 |
| `sens_log_c(dl, block, cs)` | S4: LMM na log(min + c) pro c ∈ {0,01; 0,1; 1} → poměr geometrických průměrů. | 03 |
| `sens_season(dl, block, adjust)` | S5: MMRM s adjustací na sezónu (faktor nebo roční harmoniky). | 03 |
| `season_anchor_test(dl, block, descriptor)` | Opravený test sezóny z A04/A07 (chyba K3): společný 2-df LRT harmonik, ML. | 03, 06 |
| `sens_time_spline(dl, block, at)` | S6: čas od M0 jako natural spline (df = 3), kontrasty ve 3/6/18 měsících. | 03 |
| `sens_tipping_m18(dl, block, deltas)` | S9: tipping point M18 (imputace podmíněným průměrem z MMRM, posun −δ). | 03 |
| `sens_dropout_m3(dl, block)` | S10: změna v M3 podle pozdější dostupnosti M18. | 03 |

### 3.6 `R/05_helpers_e3fast.R` — maticový výpočet E3 (bootstrap a plasmode)

Rychlá maticová verze pravidla E3. Používá se pro kalibraci a sílu (06 V7) a pro revizi (08, 10).

| Funkce | Co dělá | Volá |
|---|---|---|
| `e3_pair_matrix(comp, visit)` | Matice párových clr změn (pacienti × 19 binů); atribut `x` = log intenzity. | 06, 08, 10 |
| `e3_projector(x)` | Reziduální operátor R = I − X(X'X)⁻¹X' pro spline null df = 3. | 06 |
| `bump_index(flag, min_len)` | Indexy binů tvořících hrb (stejná definice jako `find_bumps()`). | — |
| `e3_band_decision(D, R, B, seed, level)` | Pozorovaná rezidua, max-t pás z vážených bootstrap průměrů a rozhodnutí o kladných a záporných hrbech. | 10, R06 |
| `simulate_e3(D, R, hit, size, n_sim, b_sim, seed)` | Plasmode: výběr z matice s nulovým populačním reziduem (true null), volitelně přidaný hrb, aplikace pravidla → type-I / síla. | 06 |

### 3.7 `R/06_helpers_revision.R` — funkce revize A07 (A–E)

Konstanty: `NET_GROSS_STABLE = 0,20`, `CLIN_NAME`, `SUBTYPE_MODELS` (E_a–E_e).

| Funkce | Co dělá | Volá |
|---|---|---|
| `describe_profile_norm(d)` | **A:** kvantily denního profilu v min/den na 10 mg šířky binu. | 07 |
| `top_k(n, q)` | k = ceiling(q · n), alespoň 1. | — |
| `concentration_metrics(x, q)` | **B:** pro jednu buňku bin × návštěva: % zisků, ztrát a nul, průměr, medián, net/gross, `share_net_orig` (původní estimand), `share_pos` (podíl top k na kladných ziscích) a `share_pos_ref` (matematické minimum). | 07, 10 |
| `e3_comp_grouped(dl, edges, zero, frac)` | **C:** kompozice E3 na jiné mřížce složené z nativních hran GGIR (jen slučování, bez interpolace). | 08 |
| `e3_projector_w(x, df, w)` | **C:** reziduální operátor (váženého) spline nullu s libovolným df (df = 1 je přímka). | 08, 10 |
| `e3_variant_decision(D, R, B, seed, level)` | **C:** rezidua, max-t pás a rozhodnutí pro jednu variantu × návštěvu. | 08 |
| `simulate_e3_df(D, x, df_true, df_test, hit, size, n_sim, b_sim, seed)` | **C:** plasmode volby df, včetně chybně specifikovaného nullu (true df_true, test df_test). | 08 |
| `shape_contrast(Da, Db, B, seed, keep, level)` | **D:** párové srovnání tvaru změny mezi dvěma návštěvami: rozdíly po binech s max-t pásem, rozdíl sklonu a globální bootstrap Wald test 3 ne-interceptových spline koeficientů. | 08 |
| `locate_clinical()` | **E:** najde klinický `.xlsx`: `MYOII_CLINICAL_XLSX` → `working/data/` → nadřazená složka. | 09 |
| `load_clinical(path)` | **E:** načte výchozí (M0) klinické kovariáty; „x“ a prázdné buňky → NA; sjednotí ID. | 09, 10 |
| `patient_change_summaries(acc, dl)` | **E:** souhrny změny na pacienta a návštěvu: ΔT100, clr tilt, ΔIG (+ výchozí hodnoty). | 09, 10 |
| `fit_subtype_models(d, models)` | **E:** random-intercept LMM změny (ML) podle `SUBTYPE_MODELS`; průměrné změny podtypů, kontrasty IMNM−DM, ASyS−DM, ASyS−IMNM a efekt závažnosti. | 09, 10 |
| `subtype_lrt(d)` | **E:** LRT (ML) pro subtyp × návštěva, globální subtyp a závažnost × návštěva. | 09 |

### 3.8 Lokální funkce ve skriptech

| Skript | Funkce | Co dělá |
|---|---|---|
| 01 | `describe_df()` | Řádek datového slovníku pro každou proměnnou (třída, chybění, unikáty, min/max). |
| 01 | `describe_profile()`, `describe_coverage()` | Kvantily profilu po návštěvách; % záznamů nad limity pokrytí. |
| 01 | `plot_profile()` | Figura deskriptivního profilu (medián + P10–P90). |
| 02 | `describe_set()` | Textový popis souvislých úseků confidence set (např. „70–140“). |
| 03 | `subset_dl()` | Podmnožina analytické tabulky se zachováním atributu `levels`. |
| 03 | `e3_variants_stat()` | E3 pro varianty invariance (náhrada nul, S7, S8); statistika pro bootstrap. |
| 04 | `top_share()` | Původní C3: podíl top 10 % pacientů na čistém zisku. |
| 04 | `c1_stat()` | C1: rozdíl clr změn (M6 − M0) − (M3 − M0) po binech. |
| 05 | `ttl()`, `caption_n()` | Titulek figury; popisek s počty párů a B. |
| 05 | `plot_change_profile()` | Společná šablona figur změny po intenzitách (F1, F1c, F2, F2c). |
| 06 | `vrow()`, `col_sum()`, `bins_of()`, `mc_ci()` | Řádek validační tabulky; součet sloupce; indexy binů v rozsahu mg; Monte Carlo CI (binomické). |
| 06 | `stat_fast()`, `one_loo()` | Statistika E3 pro test 1 vs N workerů; leave-one-patient-out pro E1. |
| 07 | `conc_stat()` | Statistika B pro bootstrap (`concentration_metrics` po bin × návštěva). |
| 08 | `lr_mat()` | Matice párových změn po transformaci (kontrola invariance E3). |
| 08 | `bins_of()`, `mc_ci()` | Viz 06. |
| 09 | `q3()` | Medián [IQR] (n) pro tabulku rozložení. |
| 09 | `model_data()` | Datová sada pro jeden výsledek × ukazatel závažnosti (centrování, škálování na 10 mm VAS / 5 bodů MMT8). |
| 09 | `model_set()` | Všechny modely E_a–E_e pro všechny výsledky a ukazatele závažnosti. |
| 09 | `stat_E()` | Statistika E pro pacientský bootstrap. |
| 10 | `vrow()`, `pf()` | Řádek validace; převod podmínky na PASS / WARN / FAIL. |
| 10 | `num_cols()`, `nonfinite()` | Kontrola nekonečných a NaN hodnot v tabulkách. |
| 10 | `stat_b()` | `share_pos` pro test reprodukovatelnosti bootstrapu. |
| 10 | `mdat()` | Data modelu E pro leave-one-out (R7). |
| 10 | `hash_ids()` | Otisk množiny pacientů (R8; jen popisný). |

---

## 4. Konvence kódu

- Hlavička každého souboru: PURPOSE / INPUTS / ANALYSES / OUTPUTS (viz `AGENTS.md` projektu).
- Sekce `# *** … ----` jsou kompatibilní s outline v RStudiu.
- Výstupy mají prefix data `YYMMDD_` a číslo skriptu: `261016_08_…`.
- „MVPA“ v názvech GGIR proměnných je historické označení, nikoli validovaná kategorie intenzity.
