# Ceiling/Floor Clinical Report: ig_gradient pro klinickou praxi

Datum: 2026-03-04  
Autor: AI analytik (biostatistika + klinicka interpretace)  
Hlavni zdroje:  
- `output/ceil_flo/tables/260303_ig_gradient_hurdle_conversion_01.xlsx`  
- `output/ceil_flo/tables/260304_disease_activity_models_01.xlsx`  
- `output/ceil_flo/tables/260304_phvas_ig_gradient_ceiling_models_01.xlsx`  
- `output/ceil_flo/tables/260304_res_disease_model3_posthoc_all_01.xlsx`  
- `output/ceil_flo/tables/260304_res_disease_model3_posthoc_cei_flo_01.xlsx`

## 1. Shrnuti pro klinika (1 strana)

Cil: zhodnotit, zda lze `ig_gradient` pouzit jako praktickou nahradu klinickych skal (`mmt8`, `mmt10`, `fi2`, `fi3`, `haq`) zejmena u pacientu se stropem/dnem.

Hlavni zaver podle skaly:

| Skala | Verdikt | Strucne vysvetleni |
|---|---|---|
| mmt8 | Castecne | Nejlepsi signal mimo strop, ale boundary detekce jen stredni a BA-like trend chyby je vyrazny. |
| mmt10 | Castecne | Podobne jako mmt8, navic casovy efekt (M0 vs M6) se muze menit. |
| fi2 | Omezene | Mimo strop ma vztah, ale ceiling cast je velmi mala a klasifikace stropu selhava. |
| fi3 | Omezene | Mimo strop signal pritomen, ale chyba prevodu je vysoka a ceiling cast slaba. |
| haq | Spise ne | Slaba individualni predikce, vyrazny proportional bias, v modelech PhVAS nepridava stabilni informaci. |

Prakticke doporuceni:
- `ig_gradient` je vhodny jako doplnkovy digitalni marker trendu.
- Jako samostatna nahrada klinickych skal zatim vhodny neni.
- Nejvetsi potencial ma pro `mmt8` a `mmt10`, zejmena u pacientu, kteri se blizi stropu.

## 2. Klinicky problem

Strop/dno efekt snizuje citlivost skal:
- u MMT/FI pri stropech nelze dobre zachytit dalsi zlepseni,
- u HAQ pri floor (0) nelze dobre zachytit dalsi zlepseni.

`ig_gradient` je prubezna, objektivni metrika aktivity a muze pridat informaci tam, kde klinicka skala saturuje.

## 3. Data a kohorta

### 3.1 Vstupni data

- 69 pacientu
- 1435 radku v long formatu (pacient x navsteva x outcome)
- navstevy: `M0`, `M3`, `M6`, `M18`, `M30`, `M48`
- pouzite skaly: `mmt8`, `mmt10`, `fi2`, `fi3`, `haq`
- prediktor: `ig_gradient`

Definice hranic:
- `mmt8`: strop 80
- `mmt10`: strop 100
- `fi2`: strop 100
- `fi3`: strop 100
- `haq`: floor 0 (rozsah 0-3)

### 3.2 Velikost analyzovanych kohort

All cohort (PhVAS modely):

| Outcome | n mereni | n pacientu |
|---|---:|---:|
| mmt8 | 178 | 60 |
| mmt10 | 169 | 59 |
| fi2 | 155 | 57 |
| fi3 | 166 | 59 |
| haq | 179 | 60 |

No ceiling/floor cohort:

| Outcome | n mereni | n pacientu |
|---|---:|---:|
| mmt8 | 131 | 53 |
| mmt10 | 130 | 54 |
| fi2 | 150 | 55 |
| fi3 | 148 | 55 |
| haq | 136 | 52 |

Ceiling/floor cohort (podskupina pro Model 3 post-hoc):

| Outcome | n mereni | n pacientu |
|---|---:|---:|
| mmt8 | 47 | 23 |
| mmt10 | 41 | 21 |
| fi2 | 3 | 1 |
| fi3 | 13 | 6 |
| haq | 27 | 16 |

Poznamka: `fi2` a `fi3` maji v ceiling/floor kohorte velmi maly vzorek, proto cast analyz nebyla stabilne odhadnutelna.

## 4. Analyticky plan (klinicky srozumitelne)

- EDA:
  - distribuce a trendy v case,
  - Kendall korelace `ig_gradient` vs klinicke skaly,
  - specialni vizualizace pro pacienty se stropem/floor.
- Conversion model (`ig_gradient -> klinicka skala`):
  - hurdle pristup (boundary pravdepodobnost + non-boundary hodnota),
  - grouped CV po pacientech (OOF predikce bez leakage).
- Disease activity modely (anchor `ph_vas`):
  - Model A: `ph_vas ~ clin + exam + covariates + (1|patient)`
  - Model B: `ph_vas ~ clin + ig_gradient + exam + covariates + (1|patient)`
  - Model C: `ph_vas ~ ig_gradient + exam + covariates + (1|patient)`
  - porovnani `compare_performance`, plus post-hoc mezi navstevami pro Model 3.

## 5. Vysledky podle klinicke otazky

### 5.1 Jak velky je problem stropu/dna?

Boundary prevalence v OOF conversion analyze:

| Outcome | Typ hranice | Prevalence | n boundary / n total |
|---|---|---:|---:|
| mmt8 | ceiling | 28.0% | 58 / 207 |
| mmt10 | ceiling | 25.3% | 50 / 198 |
| fi2 | ceiling | 2.7% | 5 / 182 |
| fi3 | ceiling | 11.2% | 22 / 197 |
| haq | floor | 23.8% | 50 / 210 |

Interpretace:
- klinicky nejvetsi saturace je u `mmt8`, `mmt10` a `haq`.
- `fi2` ma velmi malo ceiling pripadu, coz omezuje spolehlivou analyzu stropove casti.

### 5.2 Ma `ig_gradient` informaci mimo strop/dno?

Kendall tau (non-boundary cast, smer `ig_gradient`):

| Outcome | n | Kendall tau | p-hodnota | Delta P90-P10 (95% CI) |
|---|---:|---:|---:|---:|
| mmt8 | 149 | 0.229 | <0.001 | +5.93 (2.88; 9.30) |
| mmt10 | 148 | 0.263 | <0.001 | +9.51 (5.73; 13.46) |
| fi2 | 177 | 0.345 | <0.001 | +30.26 (20.10; 41.55) |
| fi3 | 175 | 0.270 | <0.001 | +25.11 (16.90; 34.52) |
| haq | 160 | -0.154 | 0.0049 | -0.393 (-0.639; -0.100) |

Interpretace:
- mimo boundary je vztah statisticky prukazny u vsech outcome,
- smer je klinicky konzistentni (vyssi ig_gradient = lepsi MMT/FI, nizsi HAQ).

Doporucene grafy:
- `../../output/ceil_flo/figures/260303_eda_kendall_mmt8_01.png`
- `../../output/ceil_flo/figures/260303_eda_kendall_mmt10_01.png`
- `../../output/ceil_flo/figures/260303_eda_kendall_fi2_01.png`
- `../../output/ceil_flo/figures/260303_eda_kendall_fi3_01.png`
- `../../output/ceil_flo/figures/260303_eda_kendall_haq_01.png`

### 5.3 Dokaze `ig_gradient` predikovat klinickou skalu?

OOF global metriky conversion modelu:

| Outcome | MAE | RMSE | R2 |
|---|---:|---:|---:|
| mmt8 | 5.63 | 7.22 | 0.198 |
| mmt10 | 7.50 | 9.45 | 0.217 |
| fi2 | 19.27 | 23.67 | 0.183 |
| fi3 | 21.64 | 26.28 | 0.157 |
| haq | 0.562 | 0.679 | 0.063 |

Boundary klasifikace (AUC):

| Outcome | AUC | PR-AUC |
|---|---:|---:|
| mmt8 | 0.668 | 0.429 |
| mmt10 | 0.666 | 0.399 |
| fi2 | 0.227 | 0.022 |
| fi3 | 0.662 | 0.240 |
| haq | 0.617 | 0.333 |

BA-like trend chyby (`diff ~ mean`) ukazuje proporcionalni bias u vsech outcome:
- trend_slope cca `0.84-1.19`, vsechny p << 0.001.

Prakticky dopad:
- model prinasi zlepseni oproti baseline (prumer), ale presnost nestaci na plnou substituci klinicke skaly na individualni urovni.

Decision sheet (usability):
- `mmt8`, `mmt10`: `partially_usable`
- `fi2`, `fi3`, `haq`: `limited`

### 5.4 Co se deje u pacientu se stropem/floor?

#### A) Disease activity modely (PhVAS): all cohort vs no ceiling/floor

Nejlepsi model podle AIC vah:
- All cohort: u vsech outcome vychazi nejlepe Model 2 (`clin + ig_gradient`), u HAQ je vyhoda velmi tesna.
- No ceiling/floor: Model 2 je nejlepsi pro `mmt8`, `mmt10`, `fi2`, `fi3`; pro `haq` je nejlepsi Model 1 (`clin`).

Fixni efekt `ig_gradient` v Modelu 2:
- signifikantni pro `mmt8`, `mmt10`, `fi2`, `fi3` (all i no-cei/flo),
- nesignifikantni pro `haq`.

Vizualizace modeloveho vykonu:
- `../../output/ceil_flo/figures/260304_disease_activity_compare_performance_all_cohort_01.png`
- `../../output/ceil_flo/figures/260304_disease_activity_compare_performance_01.png`

Vizualizace vztahu PhVAS ~ ig_gradient:
- `../../output/ceil_flo/figures/260304_disease_activity_phvas_ig_gradient_exam_all_cohort_01.png`
- `../../output/ceil_flo/figures/260304_disease_activity_phvas_ig_gradient_exam_no_cei_flo_01.png`

#### B) Post-hoc Model 3 (`ph_vas ~ ig_gradient + ...`)

All cohort - signifikantni rozdily mezi navstevami:
- vsechny outcome: silny rozdil `M0 -> M3` (pokles odhadu),
- navic `M6 -> M18` u `fi2`, `fi3`, `haq`.

Ceiling/floor cohort - signifikantni rozdily:
- `mmt8`: `M0-M3` (estimate -16.0; 95% CI -22.9 až -9.12; p=0.000005)
- `mmt10`: `M0-M3` (estimate -8.85; 95% CI -15.9 až -1.82; p=0.0136)
- `haq`: `M0-M3` a `M3-M6` signifikantni
- `fi2`, `fi3`: model nestabilni / nedostupny (maly vzorek)

#### C) Strict ceiling-only modely (`ph_vas ~ ig_gradient + covariates`)

- `mmt8`: beta -10.9, p=0.48
- `mmt10`: beta -13.6, p=0.25
- `fi2`, `fi3`: nedostatecna data
- `haq`: zde se model stropu nepouziva (HAQ ma floor)

Klinicka interpretace:
- v ciste ceiling-only casti neni samostatny signal `ig_gradient` dostatecne robustni,
- prakticky je vhodnejsi kombinovat `ig_gradient` s klinickou skalou a casem, ne ho pouzivat izolovane.

## 6. Rozhodovaci tabulka pro praxi

| Skala | Mira stropu/dna | Sila vazby na ig_gradient (mimo boundary) | Vykon prevodu | Vazba na PhVAS | Doporuceni |
|---|---|---|---|---|---|
| mmt8 | Vysoka (28%) | Stredni (tau 0.229) | Castecne pouzitelny | Model 2 lepsi, ig signifikantni | Nenahradit, ale silny doplnek hlavne pri stropech |
| mmt10 | Vysoka (25%) | Stredni (tau 0.263) | Castecne pouzitelny | Model 2 lepsi, ig signifikantni | Nenahradit, doplnit o ig_gradient pro lepsi sledovani trendu |
| fi2 | Nizka (2.7%) | Stredne silna (tau 0.345) | Omezene (slaba boundary cast) | Model 2 lepsi, ig signifikantni | Pouzitelne jen jako doplnek, ceiling cast nelze spolehlive hodnotit |
| fi3 | Nizka/stredni (11.2%) | Stredni (tau 0.270) | Omezene | Model 2 lepsi, ig signifikantni | Pouzitelne jen jako doplnek |
| haq | Vysoka floor (23.8%) | Slaba (tau -0.154) | Omezene | ig v Modelu 2 nesignifikantni | Pro nahradu nevhodne; drzet klinickou skalu |

## 7. Co delat v ordinaci (implementace)

Doporuceny 3-krokovy postup:

1. Pouzit klinickou skalu jako primarni rozhodovaci metriku.
2. Pokud je pacient na strope/floor (nebo blizko), doplnit interpretaci o trend `ig_gradient` mezi navstevami.
3. Pokud klinicka skala je stabilni, ale `ig_gradient` se vyznamne meni, overit to v kontextu `ph_vas`, funkcniho stavu a klinickeho obrazu.

Prakticke pravidlo:
- `ig_gradient` berte jako "early signal" zmeny aktivity,
- definitivni klinicke rozhodnuti neopirejte jen o `ig_gradient`.

## 8. Omezeni a jistota zaveru

- Pozdni navstevy (zejmena `M48`) maji maly pocet mereni.
- Ceiling/floor podskupiny jsou male, u `fi2/fi3` vede tento problem az k neodhadnutelnym modelum.
- BA-like analyza ukazuje proporcionalni bias u vsech skal, coz omezuje individualni substituci.
- Zavery jsou robustni pro tvrzeni "doplnek ano", ale opatrne pro tvrzeni "nahrada ano".

## Hlavni klinicky zaver

V teto kohorte `ig_gradient` zatim **neni vhodny jako plna nahrada** `mmt8/mmt10/fi2/fi3/haq`.  
Nejvyssi prakticky prinos ma jako **doplnek** k MMT skalam pri stropovych hodnotach a pro longitudinalni monitoring aktivity.
