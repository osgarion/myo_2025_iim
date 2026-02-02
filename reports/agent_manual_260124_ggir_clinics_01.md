# MANUÁL PRO R AGENTA: generování reportu `260124_ggir_clinics_01` (rCCA + DIABLO + MOFA)

> Cíl: vytvořit **HTML report** použitelný jako podklad pro **manuskript** i **grantovou zprávu**.  
> Report bude v **češtině** a zaměřený pouze na **M0 a M6**.  
> Subtypy jsou číselné (**1, 2, 3…**), bez jmenných labelů.

---

## 0) Kontext, role a očekávání kvality

### Role agenta
- **R programátor**: používá primárně `tidyverse`, `tidymodels`, `purrr`, `rio` (+ běžně `here`, `conflicted`, `glue`).
- **Vědec a lékař (myozitidy/revmatologie)**: píše interpretace výsledků a klinický význam.

### Nejdůležitější princip
- Analýzu **nepřepočítávat v reportu**, pokud to není nezbytné.  
  Report má primárně:
  - načíst uložené objekty z `.RData`,
  - vložit existující obrázky z `output/figures/final`,
  - načíst existující tabulky (nebo je jednorázově vygenerovat a uložit do `reports/` pro upřesnění textu).

---

## 1A) POVINNÝ PRVNÍ KROK: vytvoření sekce v hlavním `.R` skriptu pro export objektů do `.RData`

### Zásadní pravidlo
**Ještě před zahájením tvorby reportu musí agent upravit hlavní analytický `260124_ggir_clinics_01.R` skript**
(ten, který se nachází ve složce `scripts/` a je přiložen k tomuto zadání).

Cílem tohoto kroku je:
- uložit **všechny objekty**, které mohou být potřeba pro tvorbu reportu,
- do **jednoho `.RData` souboru**,
- uloženého do složky **`reports/`**,
- který bude následně **načítán v `.Rmd` reportu**.

Bez tohoto kroku **nesmí agent pokračovat** k tvorbě reportu.

---

### Název sekce v `.R` skriptu
Agent **vytvoří novou samostatnou sekci** v hlavním `.R` skriptu s názvem:


# **********************************************************************
# 7. REPORT EXPORT: objects for 260124_ggir_clinics_01 ----
# **********************************************************************


## 1) Vstupy a kde je hledat

### 1.1 Hlavní analytický skript (.R)
- Umístění: složka **`scripts/`** (uživatelem potvrzeno).
- Agent musí:
  1) **vyhledat** příslušný `.R` skript ve `scripts/`,
  2) spustit jej / zkontrolovat jeho běh,
  3) **na konec** vložit blok pro uložení všech „report-ready“ objektů do `.RData` v `reports/`.

> Pozn.: Kód v reportu má být **vzorově podobný** tomu ve skriptu (styl, názvosloví, volání funkcí).

### 1.2 Podpůrné soubory s funkcemi/objekty (povinné)
Agent musí najít a zdrojovat (source) i tyto soubory:
- `scripts/functions/OBJ_01.R`
- `scripts/functions/FUN_01.R`

Důvod: obsahují funkce/objekty používané skriptem, které nemusí být definované přímo v hlavním `.R`.

### 1.3 Styl reportu (.Rmd šablona)
Agent musí z přiloženého `.Rmd` vytáhnout a **přenést** do nového reportu:
- YAML hlavičku,
- `<style>…</style>` (CSS),
- logiku „Session info“ sekce.

Report se bude jmenovat:
- `reports/260124_ggir_clinics_01.Rmd`
a kompiluje se do HTML.

---

## 2) Změna hlavního skriptu: uložení objektů do `.RData`

### 2.1 Co přesně udělat
1) Na začátek (pokud tam není) zajistit:
   - `source("scripts/functions/OBJ_01.R")`
   - `source("scripts/functions/FUN_01.R")`
2) Na konec hlavního `.R` skriptu vložit blok:
   - vytvoří `reports/` (pokud neexistuje),
   - uloží všechny důležité objekty do jednoho `.RData`,
   - pojmenuje soubor fixně: `reports/260124_ggir_clinics_01.RData`

### 2.2 Co má být v `.RData` (minimální sada)
Ulož minimálně (podle reálných názvů objektů ve skriptu):
- Data / join tabulky pro M0/M6:
  - `d12_join`, `d13_join` (nebo ekvivalenty pro aktivitu/čas),
  - jakékoli tabulky, které report používá pro popisy.
- Modely:
  - rCCA: `rcc_fit_activ`, `rcc_fit_time`
  - DIABLO: `diablo_fit_activ`, `diablo_fit_time` (+ `perf_*` pokud existuje)
  - MOFA: `mofa_trained` (a vše potřebné pro `plot_dimred()`; např. metadata/subtype/exam_order)
- Výstupní tabulky pro interpretaci:
  - loadings, korelace, výběry feature setů, signatury, p-hodnoty apod.
- Pomocné objekty, které report očekává:
  - vektory sloupců, nastavení modelů (keepX/design), mapy názvů apod.

> Pokud je environment velký, preferuj „white-list“ výběr objektů (jen potřebné pro report).

### 2.3 Doporučený blok k vložení na konec skriptu
Agent upraví seznam objektů podle skutečných názvů ve skriptu:

```r
# ============================================================
# REPORT EXPORT: save all objects needed for HTML report
# ============================================================
dir.create("reports", showWarnings = FALSE, recursive = TRUE)

report_rdata <- file.path("reports", "260124_ggir_clinics_01.RData")

# White-list (doplnit podle reálných objektů):
objs_to_save <- c(
  # --- data ---
  "d12_join", "d13_join", "d12_mofa",
  # --- rCCA ---
  "rcc_fit_activ", "rcc_fit_time",
  # --- DIABLO ---
  "diablo_fit_activ", "diablo_fit_time",
  "perf_diablo_activ", "perf_diablo_time",
  # --- MOFA ---
  "mofa_trained", "mofa_obj_activ", "mofa_obj_time",
  # --- tables/dfs (pokud existují) ---
  ls(pattern = "^(cor_df_|load_|tab_|tbl_|res_|out_)"),
  # --- figures/plots (pokud existují jako objekty) ---
  ls(pattern = "^(fig_|plt_)")
) |> unlist() |> unique()

objs_to_save <- objs_to_save[objs_to_save %in% ls()]

save(list = objs_to_save, file = report_rdata)
message("Saved report objects to: ", normalizePath(report_rdata))
```

---

## 3) Nový report `reports/260124_ggir_clinics_01.Rmd`

### 3.1 Povinné vlastnosti
- Jazyk: **čeština**
- Formát: **HTML** (PDF není potřeba)
- Struktura: podobná šabloně `.Rmd` (YAML + CSS + Session info)
- Data: načítat přes `load("reports/260124_ggir_clinics_01.RData")`
- Zahrnout pouze **M0 a M6**
- Subtypy: používat pouze **čísla** (např. „Subtyp 1“, „Subtyp 2“…), žádná jména.

### 3.2 Setup chunk (na začátku)
- `knitr::opts_chunk$set(echo = TRUE/FALSE dle potřeby, message = FALSE, warning = FALSE)`
- načíst knihovny (minimálně `tidyverse`, `purrr`, `rio`, a balíčky pro grafy dle skriptu)
- `load("reports/260124_ggir_clinics_01.RData")`

### 3.3 Vkládání obrázků
- Primární zdroj obrázků: `output/figures/final`
- Agent má:
  - vylistovat soubory v `output/figures/final`,
  - přiřadit je k sekcím rCCA/DIABLO/MOFA podle názvů/prefixů,
  - vkládat je přes:

```r
knitr::include_graphics("output/figures/final/NAZEV_SOUBORU.png")
```

> Cíl: co nejvíc využít již hotové obrázky; nevyrábět znovu, pokud existují.

### 3.4 Tabulky
- Pokud skript exportuje tabulky (CSV/XLSX) do `output/tables/final`, report je načte a použije pro:
  - přesný popis (např. top loadings, top korelace, signatura feature setů).
- Pokud tabulka není exportovaná, agent ji může:
  - jednorázově z objektů v `.RData` sestavit,
  - uložit do `reports/` (dočasně i finálně),
  - a vložit do reportu jako `knitr::kable()` / `gt` / `DT`.

---

## 4) Povinný obsah reportu podle metod

Každá metoda musí mít:
1) **Vysvětlení pro neodborníky** (biologové a lékaři), česky, s příkladem z klinické praxe  
2) **Výstupy** (obrázky + relevantní tabulky)  
3) **Interpretaci pod každým výstupem**  
4) **Význam pro myozitidy** (klinický/biologický kontext)

### 4.1 rCCA (odděleně: AKTIVITA a ČAS)
- Report musí držet paralelní strukturu:
  - **AKTIVITA outputs**
  - **ČAS outputs**
- U každé části:
  - canonical scores/variates,
  - heatmap korelací X vs Y,
  - loadings (které proměnné táhnou osu),
  - popis a klinická interpretace pro myozitidy.

### 4.2 DIABLO
- Vysvětlit integrativní multi-blok přístup.
- Zahrnout:
  - separace subtypů v latentním prostoru,
  - selected features (pokud existuje cim/heatmap),
  - stručně: co je „signatura“ a jak se dá interpretovat u myozitid.

### 4.3 MOFA
- Vysvětlit latentní faktory (unsupervised).
- Zahrnout dim-reduction (UMAP) minimálně:
  - barvení podle `subtype`
  - barvení podle `exam_order` (M0 vs M6)
- Interpretovat faktory:
  - co naznačují o biologii/klinice,
  - jaké proměnné se typicky vážou k faktorům (pokud je v datech dostupné).

---

## 5) Interpretace výstupů (povinný styl textu)

Pod každým grafem/tabulkou:
- „Co vidíme“ (popis)
- „Jak to číst“ (navigace pro čtenáře)
- „Co to znamená pro myozitidy“ (klinický význam, hypotézy)
- „Co je nejisté“ (omezení, potřeba validace)
- „Jaký další krok“ (validace, doplňkové analýzy/experimenty)

---

## 6) Finální sekce: revmatologické zhodnocení (detailní)

Na konci reportu musí být samostatná kapitola:
- syntéza napříč rCCA/DIABLO/MOFA,
- co se shoduje vs. co je specifické pro metodu,
- implikace pro myozitidy (fenotypy/subtypy, klinická stratifikace),
- doporučení pro manuskript (co dát do main vs. supplement),
- doporučení pro grantovou zprávu (srozumitelné shrnutí dopadu).

---

## 7) Reprodukovatelnost

V závěru reportu:
- sekce „Session info“ ve stylu šablony:
  - `devtools::session_info()` nebo ekvivalent,
  - doplnit verze balíčků klíčových pro analýzu.

---

## 8) Co agent musí odevzdat (artefakty)

1) Upravený `.R` skript ve `scripts/` s přidaným exportem:
   - `reports/260124_ggir_clinics_01.RData`
2) Nový report:
   - `reports/260124_ggir_clinics_01.Rmd`
   - render do HTML

---

## 9) Kontrolní seznam před odevzdáním

- [ ] `.RData` existuje v `reports/` a jde načíst bez chyb
- [ ] report se renderuje do HTML bez přepočtu celé analýzy
- [ ] všechny obrázky se načítají z `output/figures/final` (kde je to možné)
- [ ] každá metoda má vysvětlení + klinický příklad + interpretace výstupů
- [ ] závěrečné revmatologické zhodnocení je detailní a česky
- [ ] report obsahuje session info

Konec manuálu.
