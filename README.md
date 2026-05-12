# myo_2025_iim

## Revize analýzy — 2026-05-07

### Kontext

Kolegyňa Lucia Vernerová upravila vstupní data: ze souboru jsou vyjmuti pacienti, kteří byli příliš dlouho léčeni nebo nemocní (příliš dlho chorí/liečení). Revidované listy mají v názvu příponu `_revize`.

**Dotčené zdrojové soubory** (`R:\MYOZITIDY_VÝZKUM\_IMMET-GRANT\IMMET štatistika\Data\Raw\Data\`):

- `IMMET_terapie a klinická data_01102025.xlsx` — nový list `*_revize`
- `Myokiny sérum_pro stat.xlsx` — nový list `*_revize`

Tyto listy je třeba použít místo původních ve všech nových analýzách.

### Analýza 1 — Hlavní tabulka lmer kovariantních modelů M0 vs M6

- **Skript:** `scripts/Script_myo_2025_iim_02.R` (nový, odvozený od `Script_myo_2025_iim_01.R`)
- **Metoda:** lmer pro všechny kombinace závislých proměnných (`var_dep_01`) a nezávislých proměnných (`var_indep_01`: mstn, fst, fstl3, akta), pouze časové body M0 a M6
- **Varianty konfounderu:** bez konfounderu (`_without`), věk (`_age`), denní dávka GK na kg TH (`_gk`), BCM (`_bcm`), kreatinin (`_creatine`)
- **Výstupní sloupce:** beta, standardizovaná beta (`refit`/`posthoc` fallback), CI, p-hodnota, partial omega², podmíněné R²
- **Vzorová tabulka:** `R:\MYOZITIDY_VÝZKUM\_IMMET-GRANT\IMMET štatistika\Output\251017\m0_m6\260424_covariates_0_m6_misc_01.xlsx`
- **Výstup skriptu:** `output/tables/YYMMDD_covariates_revize_01.xlsx`

### Analýza 2 — Post-hoc subgrouping (vizualizace)

Pro každý signifikantní kovariantní model (p < 0,06) se výsledky zobrazují zvlášť podle čtyř stratifikačních proměnných:

| Proměnná | Popis |
| --- | --- |
| `pohlavi` | pohlaví |
| `podtyp_nemoci_zjednoduseny` | poddiagnóza |
| `jo_1` | přítomnost anti-Jo1 protilátky |
| `anti_hmgcr` | přítomnost anti-HMGCR protilátky |

Subgrouping se provádí zvlášť pro každý kovariantní model (věk, GK dávka, BCM, kreatinin). Vzorové výstupy (z původního datasetu):
`R:\MYOZITIDY_VÝZKUM\_IMMET-GRANT\IMMET štatistika\Output\251017\m0_m6\{age,bcm,gk,creatine}\{sex,subdiag,jo_1,anti_hgmcr}.pdf`

Nové výstupy: `output/figures/YYMMDD_revize_{covariate}_{stratification}.pdf`

### Analýza 3 — PLS-DA M0, podskupiny MDA kategorií

- **Skript:** sekce 6 v `scripts/Script_myo_2025_iim_02.R`
- **Metoda:** `mixOmics::plsda()`, ncomp = 2, scale = TRUE
- **Data:** `d03_rev` filtrovaný na M0 → 57 pacientů, všichni mají `mda_categories`
- **Výstupy:** `output/figures/cluster analyses/YYMMDD_revize_m0_pls_mda_categories_{group,loadings}.tiff`
- **Vzorové výstupy z původního datasetu:** složky `251113` a `251120` na R:\ disku

#### Výběr prediktorů

Jako prediktory vstupují všechny klinické a myokinové proměnné z `var_dep_01` a `var_indep_01`, s výjimkou:

| Proměnná | Důvod vyřazení |
| --- | --- |
| `muscle_disease_activity` | zdrojová kontinuální proměnná pro `mda_categories` — data leakage, triviálně by separovala skupiny |
| `odpoved_na_terapii_m0_vs_m6` | outcome léčby, v M0 neznámý |

Log-transformace před PLS-DA (pravostranně zešikmené markery): `ast`, `alt`, `ck`, `crp`, `haq`, `ld`, `mitax`, `myoact`, `myoglobin` → `log(x + 1)`.

Přejmenování os v loadings grafu: načítá se z listu `260507_MDA categories` v `IMMET_terapie a klinická data_01102025.xlsx`. Pokud list neexistuje, skript pokračuje s původními názvy sloupců.

#### Ošetření chybějících hodnot

Bez ošetření by `complete.cases()` ponechal pouze **22 ze 57 pacientů**. Přehled NA v klíčových prediktorech:

| Proměnná | NA (n) | NA (%) |
| --- | --- | --- |
| `sf36_mcs`, `sf36_pcs` | 16 | 28 % |
| `sf36_gh` | 15 | 26 % |
| `sf36_vt`, `sf36_re`, `sf36_mh` | 13 | 23 % |
| `sf36_pf`, `sf36_rp`, `sf36_bp`, `sf36_sf` | 12 | 21 % |
| `anti_hmgcr` | 10 | 18 % |
| `fi_2`, `borg10` | 8 | 14 % |
| `ld` | 6 | 11 % |
| `fstl3` | 5 | 9 % |

Různí pacienti mají NA v různých proměnných → jejich kombinace vyřadí 35 pacientů.

**Postup:**

1. **Práh 40 %** — proměnné s více než 40 % NA jsou vyřazeny (pojistka pro strukturální absenci; v aktuálních datech žádná proměnná tento práh nepřekračuje, max. je 28 % u SF-36).
2. **Mice imputace** — `mice::mice(method = "pmm", m = 1, maxit = 10, seed = 42)` doplní zbývající NA. Finální dataset pro PLS-DA má **57 kompletních řádků**.

### Report — `reports/260511_iim_revize_02_report.Rmd`

- **Výstup:** `reports/260511_iim_revize_02_report.html` (generován automaticky na konci `Script_myo_2025_iim_02.R` přes `rmarkdown::render()`)
- **Obsah:**
  1. **Demografická tabulka** — popisná statistika revidovaného souboru (původní data, bez mice imputace), stratifikovaná dle M0 a M6 (`gtsummary::tbl_summary()`)
  2. **Kovariantní modely M0 vs M6** — přehled počtu signifikantních asociací (p < 0,05 / p < 0,06) pro každý kovariant + hlavní tabulka (β, std β, p) s víceúrovňovými záhlavími
  3. **Post-hoc vizualizace** — scatter ploty pro každý signifikantní kovariantní model (p < 0,06), stratifikované dle pohlaví, poddiagnózy, anti-Jo1 a anti-HMGCR
  4. **PLS-DA (M0, mda_categories)** — skupinový graf a loadings graf z modelu `pls_m0_rev_mda`
- **Zdrojová data:** RData soubor `reports/YYMMDD_covariates_revize_01.RData` (obsahuje `res_covar_02_rev_tab`, `sub_age`, `sub_gk`, `sub_bcm`, `sub_creatine`, `pls_m0_rev_mda`, `d04_sel2_rev`)
- **Inspirace:** struktura reportu vychází ze šablony `reports/260504_iim_phase_angle_bootstrap_report.Rmd`

---

## Poznámky ke konkrétním výstupům (před revizí dat)

Tato sekce dokumentuje výstupy vzniklé **před** revizí vstupních dat (tj. z původního, nerevidovaného souboru). Pro revidované výstupy viz sekci "Revize analýzy — 2026-05-07" výše.

### `output/tables/260106_covariates_0_m6_misc_01.xlsx` / `output/tables/260424_covariates_0_m6_misc_01.xlsx`

- **Zdrojové skripty:** `misc/script_covariates_m0_m6_01.R`, `scripts/Script_myo_2025_iim_01.R`
- **Obsah:** Výsledky lineárních smíšených modelů (lmer, M0 vs. M6) pro všechny kombinace závislých a nezávislých proměnných, ve čtyřech variantách dle konfounderu: bez konfounderu (`_without`), s věkem (`_age`), s denní dávkou GK (`_gk`), s BCM (`_bcm`) a s kreatininem (`_creatine`). Každá varianta obsahuje partial omega², podmíněné R², nestandardizovanou betu, standardizovanou betu (metoda `refit`), CI a p-hodnotu pro `var_indep_value`.
- **Použití v článku:** Data z této tabulky byla použita pro článek v rámci **major revision** (verze ze dne 2026-01-06, aktualizovaná verze s přidanými standardizovanými betami ze dne 2026-04-24).

## Změny provedené 2026-04-24 / 2026-04-27

### `misc/script_covariates_m0_m6_01.R` — přidány standardizované bety

- Do funkcí `model_min_2` a `model_covar_2` byl přidán výpočet standardizované bety pomocí `parameters::standardize_parameters(mod1$fit, method = "refit")` s fallbackem na `method = "posthoc"` pro skupiny, kde refit selže (malý N, nulová variance).
- Nový sloupec `std_beta_covar` byl přidán do návratové tibble i do všech pěti `rename_with` bloků (`_without`, `_age`, `_gk`, `_bcm`, `_creatine`).
- Paralelizace přepsána z `multidplyr` na `furrr::future_map` + `future::plan(multisession)` — původní `multidplyr::new_cluster()` selhávalo timeoutem při startu 31 worker procesů.

### `.Rprofile` — trvalá nastavení session

Přidána volba `options()` platná pro každou R session i pro furrr workery:

- `renv.config.auto.snapshot = FALSE` — zastaví automatické skenování závislostí při startu, které způsobovalo timeout workerů.
- `myo_ncores = 31L` — nastaví počet jader (dříve bylo ořezáno na 4 výchozím `min(4L, ...)`).
- `myo_pacman_update = FALSE` a `myo_pacman_install = FALSE` — zakáže automatické updaty balíčků při každém `source()` v `FUN_01.R` (updaty přes `pacman` způsobily corrupted `cli` package).

### `.renvignore` — nový soubor

Vytvořen `.renvignore` vylučující `output/`, `data/`, `scripts/back_up/` a `misc/` ze skenování závislostí renv, což dále zkracuje inicializaci nových R session.

## Ostatní skripty

### `scripts/lmer_srovnani_2_promennych_interakce_bootstrap_bmi.R`

- **Účel skriptu:** Bootstrapová analýza lineárního smíšeného modelu (`lmer`) pro GGIR proměnnou `ig_gradient_enmo`, zaměřená na vztah k času vyšetření, věku, pohlaví a mezi-pacientským složkám `bcm` a `bmi`.
- **Model:** Hlavní model má tvar `ig_gradient_enmo ~ exam + age + sex + bcm_between + bmi_between + (1 | boot_id)`. Skript fituje také standardizovanou variantu modelu, kde jsou spojité prediktory a volitelně i odpověď převedeny na z-skóry.
- **Bootstrap:** V každém bootstrap cyklu se náhodně vybere 50 pacientů, defaultně s návratem. Pro každý cyklus se znovu vytváří within-between složky `bcm` a `bmi`, aby bootstrap respektoval opakovaná měření a případné opakované vybrání stejného pacienta.
- **Paralelizace:** Bootstrap běží přes `future::multisession` a `furrr::future_map_dfr`, počet workerů je nastaven podle dostupných jader.
- **Výstupy:** Skript ukládá full-sample model, souhrny bootstrapu, všechny bootstrap koeficienty, status úspěšnosti běhů, případné chyby, Excel report, CSV/RDS soubory a grafy hustot a forest ploty pro raw i standardizovaný model.
