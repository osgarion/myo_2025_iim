# Claude Code Instructions — myo_2025_iim

This file is automatically loaded by Claude Code at project startup.
All coding style, analysis philosophy, and project conventions are defined in [AGENTS.md](AGENTS.md).

<!-- Import AGENTS.md content -->
@AGENTS.md

## Communication language

Respond in Czech unless the user writes in English.

## Project context (quick reference)

- **Project:** myo_2025_iim — vliv fyzické aktivity (aktigrafie / GGIR) a klinických proměnných na výsledky pacientů s idiopatickými zánětlivými myopatiemi (IIM), sledování M0 → M6
- **Hlavní výsledky (outcomes):** SF-36 PCS/MCS, funkce svalů (MMT, HAQ, FI), fázový úhel (phase angle)
- **Hlavní prediktory:** GGIR aktigrafické proměnné (ig_gradient_enmo, m5_enmo, p9931_enmo, mvpa_t100_time_min, acc_day_spt_wei), klinické kovariáty (věk, pohlaví, BCM, BMI, denní dávka GK, kreatinin)
- **Hlavní metoda:** lmer (lme4) + bootstrap resampling (furrr / future::plan(multisession)), standardizované bety via `parameters::standardize_parameters(method = "refit")`
- **Paralelizace:** `getOption("myo_ncores")` (nastaveno v `.Rprofile` na 31L)

### Klíčové skripty

- `scripts/Script_myo_2025_iim_01.R` — hlavní analytický skript (lmer modely M0 vs. M6)
- `scripts/lmer_srovnani_2_promennych_interakce_bootstrap_bmi.R` — bootstrapová analýza lmer pro GGIR proměnné
- `misc/script_covariates_m0_m6_01.R` — lmer modely pro všechny kombinace závislých / nezávislých proměnných (5 variant konfounderu)
- `scripts/Script_myo_2025_iim_02_ggir.R` — zpracování GGIR dat
- `scripts/Script_myo_2025_iim_tis_01.R` — analýzy TIS

### Klíčové funkce

- `scripts/functions/FUN_01.R` — hlavní balík funkcí projektu (načítání dat, lmer wrappery, aj.)
- `scripts/functions/OBJ_01.R` — sdílené objekty (labels, palety, konfigurace)
- `scripts/functions/FUN_ceil_flo_01.R` — funkce pro analýzu ceiling/floor efektů

### Klíčové reporty

- `reports/260124_ggir_clinics_01.Rmd` — GGIR vs. klinické proměnné (hlavní přehledový report)
- `reports/260504_iim_phase_angle_bootstrap_report.Rmd` — nejnovější report (fázový úhel, bootstrap lmer); **použij jako šablonu pro nové `.Rmd` reporty**

### Výstupní složky

- `output/figures/` — obrázky (drafty); `output/figures/final/` — finální verze
- `output/tables/` — tabulky (Excel / CSV)
- Témata výstupů mají vlastní podsložky: `output/tis/`, `output/sf36/`, `output/ggir/`, `output/ceil_flo/`, atd.

## Modelovací konvence

- Opakovaná měření: M0 (baseline) a M6 (follow-up), proměnná `exam`
- Between/within dekompozice pro BCM a BMI: `bcm_between`, `bcm_within`, `bmi_between`, `bmi_within`
- Bootstrap: náhodný výběr 50 pacientů s návratem; between/within složky se přepočítávají v každém cyklu
- Standardizovaná beta: `parameters::standardize_parameters(mod, method = "refit")`, fallback na `"posthoc"` při selhání refitu (malý N, nulová variance)
- Výsledkové tabulky obsahují: estimate, std_beta, CI, p-hodnotu, partial omega², podmíněné R²
