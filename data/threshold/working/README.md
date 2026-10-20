# MYO-II / Paper A — analýza intenzitního spektra ENMO (thresholds)

**Autor:** Veronika Balounová. Kód připravil Claude Opus 5.5, nezávislou oponenturu provedl Claude Fable 5.1.

**Datum:** 2026-10-13

**Otázka:** Které části spektra intenzit ENMO (20–400 mg) nesou reprodukovatelnou longitudinální změnu běžné pohybové aktivity u pacientů s IIM (M0 → M3 / M6 / M18)?

**Analytický plán:** uzamčený plán v2 je v `docs/261013_phase3_disposition_and_locked_plan.md`.

---

## 1. Obsah složky

```
working/
├── README.md                     ← tento soubor
├── SKRIPTY_A_FUNKCE_CZ.md        ← přehled skriptů (účel, vstupy, výstupy, doba) a katalog všech funkcí
├── run_all.R                     ← spustí 01 → 10 a vyrenderuje všechny tři reporty
├── 01_QC_data_dictionary_BOTH.R  ← datový slovník, QC (PASS/WARN/FAIL), flow, deskripce
├── 02_PRIMARY_bootstrap_BOTH.R   ← E1–E4 + cluster bootstrap (B = 2000), pásy, rozhodnutí
├── 03_SENSITIVITY_BOTH.R         ← S1–S10 + invariance rozhodovacího pravidla
├── 04_DESCRIPTIVE_heterogeneity_BINS.R ← C1, C3, C5, podtypy, M30–M60
├── 05_FIGURES_TABLES_BOTH.R      ← publikační figury + legendy
├── 06_VALIDATION_BOTH.R          ← testy, diagnostika, simulace síly, issue log
├── 07_REV_A_B_profile_concentration_BINS.R ← revize A07: A profil na 10 mg, B koncentrace zisků
├── 08_REV_C_D_E3_shape_BOTH.R    ← revize A07: C varianty a invariance E3, D tvar změny M3/M6/M18
├── 09_REV_E_subtypes_clinical.R  ← revize A07: E podtypy + klinická závažnost (exploratorní)
├── 10_REV_VALIDATION_BOTH.R      ← validace revize R1–R10
├── R/                            ← sdílené funkce (pouze source, nespouštět samostatně)
│   ├── 00_config.R               ← cesty, balíčky, konstanty (B, seed, workery, subtypy), save_public()
│   ├── 01_helpers_data.R         ← načtení masteru, QC, analytické tabulky, výstupy
│   ├── 02_helpers_bootstrap.R    ← pacientský cluster bootstrap (checkpointy), max-t pásy
│   ├── 03_helpers_estimands.R    ← estimandy E1–E4, rozhodovací pravidlo
│   ├── 04_helpers_sensitivity.R  ← modely citlivostních analýz
│   ├── 05_helpers_e3fast.R       ← maticový výpočet E3 (bootstrap pás, plasmode type-I / síla)
│   └── 06_helpers_revision.R     ← funkce revize A07 (A–E)
├── reports/
│   ├── 261013_thresholds_report.Rmd      ← report v1 (jen načítá výsledky, nic nepočítá)
│   ├── 261016_thresholds_report_v2.Rmd   ← revidovaný report (EN)
│   ├── 261016_thresholds_report_v2_cz.Rmd ← revidovaný report (CZ)
│   └── fonts_header.tex                  ← fonty pro PDF výstup
├── docs/                         ← audit trail (audit, oponentura, vypořádání, slovník)
└── output/                       ← generováno skripty (RData, tables, figures, logs)
```

Označení v názvu skriptu: **BINS** = exkluzivní biny, **CUMULATIVE** = kumulativní prahy, **BOTH** = obojí.

---

## 2. Požadavky

- **R ≥ 4.4.** Používá se základní operátor `%||%`. Původní výstupy vznikly v R 4.5.2 (Windows).
- **Balíčky pro analýzu:**
  `dplyr tidyr purrr tibble stringr ggplot2 scales patchwork lme4 glmmTMB future furrr rio writexl`
- **Balíčky pro report:**
  `rmarkdown knitr kableExtra sessioninfo pacman`. Report pomocí `pacman` načítá také `tidyverse conflicted renv janitor gridExtra`.
- **PDF výstup reportu (volitelný):** XeLaTeX (např. TinyTeX).
- Chybějící balíček skript **neinstaluje sám**. Zastaví se a vypíše přesný příkaz `install.packages(...)`.

Instalace všeho najednou:

```r
install.packages(c("dplyr","tidyr","purrr","tibble","stringr","ggplot2","scales","patchwork",
                   "lme4","glmmTMB","future","furrr","rio","writexl",
                   "rmarkdown","knitr","kableExtra","sessioninfo","pacman","tidyverse",
                   "conflicted","janitor","gridExtra"))
```

> **renv:** Uvnitř projektu `myo_2025_iim` se aktivuje projektové renv přes `.Rprofile`. Pokud v něm některý balíček (např. `glmmTMB`) chybí, doinstalujte ho do renv (`renv::install("glmmTMB")`). Na jiném stroji, kde je složka `working/` samostatně, stačí běžná knihovna.

---

## 3. Kde má být master `.rds`

Skripty hledají `MYOII_ACC_OUTCOMES_MASTER_V1.rds` v tomto pořadí:

1. proměnná prostředí `MYOII_MASTER_RDS`, např. `Sys.setenv(MYOII_MASTER_RDS = "R:/.../MYOII_ACC_OUTCOMES_MASTER_V1.rds")`;
2. `working/data/MYOII_ACC_OUTCOMES_MASTER_V1.rds`;
3. nadřazená složka `working/` (v projektu je to `data/threshold/`).

Cesta a **md5** použitého souboru se ukládají do všech logů i do reportu. Master, nad kterým byl kód připraven, má md5 `a4821db9803a46c31eb4dc6639add5b9`. Vytvořil ho P06 V3 dne 2026-10-12 14:08. Jiný md5 znamená jiný vstup.

Složka `working/` se najde sama, pokud skript otevřete v RStudiu nebo nastavíte pracovní adresář do `working/`. Jinak nastavte `Sys.setenv(MYOII_WORKING_DIR = "<cesta>/working")`.

---

## 4. Spuštění

**Doporučený postup v RStudiu:**

1. Otevřete `run_all.R` a dejte **Source**.
2. Nebo spusťte jednotlivé skripty v pořadí **01 → 02 → 03 → 04 → 05 → 06**.
3. Potom knitujte `reports/261013_thresholds_report.Rmd`.

**Rychlý testovací běh** (ověří, že vše funguje, za pár minut):

```r
Sys.setenv(MYOII_B = 30, MYOII_NSIM = 5, MYOII_BSIM = 30)
source("run_all.R")
```

Výsledky testovacího běhu nejsou interpretovatelné. Před ostrým během smažte `output/RData/boot_chunks/`, případně ho spusťte v nové R session bez `MYOII_*`.

**Ostrý běh:** bez nastavení `MYOII_*` (B = 2000, N_SIM = 500, B_SIM = 1000).

| Skript | Odhad doby (31 logických jader → 30 workerů) |
|---|---|
| 01 | < 1 min |
| 02 | **~1,5–2 h** (Windows, 31 workerů: ~5,5 min na blok 100 replikací, 20 bloků) |
| 03 | ~10–15 min (Tweedie a sezóna jednovláknově; invariance-bootstrap rychlý) |
| 04 | ~2 min |
| 05 | < 1 min |
| 06 | ~5–10 min (leave-one-out 71 × 23 MMRM + kalibrace a síla pravidla E3: 8 scénářů × 3 návštěvy × 500 × 1000, maticově) |

- **Počet workerů** = `min(detectCores(logical = TRUE) − 1, getOption("myo_ncores"))`, fallback 1.
- **Reprodukovatelnost:** výběry pacientů pro bootstrap se generují předem v hlavním procesu (seed 20261013). Výsledky jsou proto **identické při libovolném počtu workerů** (ověřuje test V4).
- **Checkpointy:** bootstrap se ukládá po blocích do `output/RData/boot_chunks/`. Přerušený běh při opětovném spuštění naváže. Klíč checkpointu obsahuje B, seed, seznam pacientů a označení analýzy, takže se nemíchají běhy s jiným nastavením.

---

## 5. Výstupy

| Kam | Co |
|---|---|
| `output/RData/*_objects.RData` | Objekty pro report (qc, primary, sensitivity, descriptive, validation, figure) |
| `output/RData/primary_boot_replicates.rds` | Všechny bootstrap replikace primární analýzy |
| `output/tables/YYMMDD_0N_*.xlsx` | Tabulky každého kroku |
| `output/figures/YYMMDD_05_F*.png/.tiff` | Publikační figury (TIFF 600 dpi, LZW) |
| `output/YYMMDD_05_figure_legends.md` | Legendy figur |
| `output/logs/*_sessionInfo.txt` | Vstup + md5, B, seed, workery, `sessionInfo()` |
| `reports/261013_thresholds_report.html` | Report |

Výstupy **neobsahují řádky ani identifikátory pacientů** (jen agregované výsledky), protože se commitují do veřejného repozitáře. Master `.rds` se do gitu nedává.

---

## 6. Stav ověření (k 2026-10-13)

| Krok | Stav |
|---|---|
| 01 QC | **Spuštěno nad masterem** (R 4.6.1, Linux): 0 FAIL, 4 WARN |
| 02–06 | **Technicky otestováno** s malým B (16) a malou simulací. Kód proběhl bez chyb. **Čísla z testů nejsou výsledky.** |
| Ostrý běh (B = 2000) | **NOT EXECUTED / NOT VERIFIED.** Spouští uživatel. |
| Report | Struktura připravena dle `RULE_rmd_report_style.md`. Text závěrů se generuje z objektů. |
| Nezávislý oponent | Claude Fable 5.1: dvě kola, obě REVISE, všechny připomínky vypořádány (docs/). Po ostrém běhu doporučeno třetí kolo nad výsledky. |

---

## 7. Otevřené požadavky (mimo data v masteru)

1. **Důvody chybění** M0 (8 pacientů) a M18; CONSORT-like flow z klinické kohorty.
2. **GGIR `part2_daysummary*.csv` + `config.csv`** pro modul reliability (ICC / SEM / SDC) a pro dokumentaci kritérií validních dnů.
3. **Kontrola značek L/R** ve jménech souborů u 10 pacientů. Náramek je podle uživatele vždy na nedominantním zápěstí.

---

## 8. Kódování a konvence

- Návštěvy: M0 (reference), M3, M6, M18; M30 / M48 / M60 jen deskriptivně.
- Čáry návštěv ve figurách: **M0** čerchovaná, **M3** tečkovaná, **M6** přerušovaná, **M18** plná (vhodné pro ČB tisk).
- Subtypy (jen pro deskriptivní analýzu podtypů): 0 + 4 = DM, 2 + 3 = IMNM, 5 = ASyS. PM (1) a UCTD/SSc jsou vyřazeni z analýz podtypů, v hlavních analýzách zůstávají.
- „MVPA“ v názvech GGIR proměnných je historické označení, nikoli validovaná kategorie intenzity.

---

## 9. Revize A07 (2026-10-16)

Podklad: `../MYO_A07_narhy_na_zmeny.txt`. Původní skripty 01–06 a report `261013_*` zůstávají beze změny. Revize je ve verzovaných souborech:

| Soubor | Obsah |
|---|---|
| `R/06_helpers_revision.R` | funkce revize (A–E) |
| `07_REV_A_B_profile_concentration_BINS.R` | A: profil na 10 mg; B: koncentrace zisků (top 10 % kladných zisků, net/gross, původní estimand) |
| `08_REV_C_D_E3_shape_BOTH.R` | C: invariance a varianty E3, plasmode volby df; D: párové srovnání tvaru M3/M6/M18 |
| `09_REV_E_subtypes_clinical.R` | E: podtypy s adjustací na ph_VAS / MMT8 (exploratorní) |
| `10_REV_VALIDATION_BOTH.R` | validace R1–R10 |
| `reports/261016_thresholds_report_v2.Rmd` | revidovaný report (anglicky) |
| `reports/261016_thresholds_report_v2_cz.Rmd` | revidovaný report v češtině (stejné výsledky; překlad popisků figur a tabulek v reportu, desetinná čárka) |
| `docs/261016_A07_*.md` | kritické review, validace + change log + shrnutí, nezávislá oponentura |

**Klinický soubor** `klinicka_data_12102026.xlsx` (jen pro 09) se hledá v pořadí: `Sys.setenv(MYOII_CLINICAL_XLSX = "...")` → `working/data/` → nadřazená složka (`data/threshold/`). Do výstupů jdou jen agregované hodnoty. Navíc je potřeba balíček `readxl`.
