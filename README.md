# myo_2025_iim

## Poznámky ke konkrétním výstupům

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
