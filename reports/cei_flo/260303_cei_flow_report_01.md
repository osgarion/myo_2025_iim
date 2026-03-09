---
output:
  pdf_document: default
  html_document: default
editor_options:
  markdown:
    wrap: 72
---

# Ceiling/Floor report: ig_gradient jako nahrada klinickych skal

Datum: 2026-03-03  
Vstupni soubor: `output/ceil_flo/tables/260303_ig_gradient_hurdle_conversion_01.xlsx`

## 1) Co bylo analyzovano

Cil byl zhodnotit, zda lze `ig_gradient` prakticky pouzit jako nahradu
klinickych skal:

- `mmt8`
- `mmt10`
- `fi2`
- `fi3`
- `haq`

Pouzit byl hurdle pristup (2 casti):

1. Model pravdepodobnosti stropu/dna (`ceiling`/`floor`) pres logistickou regresi.
2. Model hodnot mimo strop/dno (`non-boundary`) pres regresi se spline.

Validace je out-of-fold (grouped CV podle `patient_id`), takze predikce
nejsou kontaminovane stejnym pacientem mezi train/test foldy.

## 2) Kvalita a rozsah dat

- 69 pacientu, 287 navstev.
- Validni dvojice (`ig_gradient + outcome`):
  - `mmt8`: 207
  - `mmt10`: 198
  - `fi2`: 182
  - `fi3`: 197
  - `haq`: 210
- Vek (`age_m0_cov`) chybi minimalne (1 zaznam v kazdem outcome bloku).
- `exam_time` nechybi.
- Outliery `ig_gradient` byly reportovany, neodstranovany (cca 4-5 % podle outcome).
- Hrance byly nastavene klinicky korektne (`mmt8` max 80, `mmt10/fi2/fi3` max 100, `haq` 0-3).

## 3) Jak cist metriky (prakticky)

- `MAE`, `RMSE`: velikost chyby v jednotkach dane skaly.
- `R2`: kolik variability skaly model vysvetli.
- `AUC`, `PR-AUC`: jak dobre model rozpoznava strop/dno.
- BA-like (`bias`, `LoA`, `trend_slope`):
  - `bias` blizko 0 je dobre,
  - `trend_slope` blizko 0 je dobre (bez proporcionalni chyby),
  - vysoke absolutni `trend_slope` znamena systematicke pod/nadhodnocovani podle uroven skaly.

## 4) Hlavni vysledky: OOF predikce (global)

| Outcome | N | MAE | RMSE | R2 | Spearman | Kendall |
|---|---:|---:|---:|---:|---:|---:|
| mmt8  | 207 | 5.63 | 7.22 | 0.20 | 0.43 | 0.31 |
| mmt10 | 198 | 7.50 | 9.45 | 0.22 | 0.46 | 0.33 |
| fi2   | 182 | 19.27 | 23.67 | 0.18 | 0.46 | 0.32 |
| fi3   | 197 | 21.64 | 26.28 | 0.16 | 0.42 | 0.29 |
| haq   | 210 | 0.56 | 0.68 | 0.06 | 0.23 | 0.16 |

Interpretace:

- Nejlepsi globalne vychazi `mmt8` a `mmt10`, ale R2 je stale jen nizke az nizko-stredni.
- `fi2`/`fi3` maji velkou absolutni chybu.
- `haq` ma slabou individualni vysvetlovaci schopnost (`R2 = 0.06`).

## 5) Boundary cast (strop/dno)

| Outcome | Typ | Prevalence | AUC | PR-AUC |
|---|---|---:|---:|---:|
| mmt8  | ceiling | 0.280 | 0.668 | 0.429 |
| mmt10 | ceiling | 0.253 | 0.666 | 0.399 |
| fi2   | ceiling | 0.027 | 0.227 | 0.022 |
| fi3   | ceiling | 0.112 | 0.662 | 0.240 |
| haq   | floor   | 0.238 | 0.617 | 0.333 |

Interpretace:

- Diskriminace stropu/dna je celkove slaba.
- `fi2` je boundary casti prakticky nevyuzitelne (velmi nizka prevalence + velmi nizky AUC).
- U `mmt8/mmt10` je boundary signal pritomny, ale nestaci na spolehlive rozhodovani.

## 6) BA-like chyba prevodu (y - y_hat)

Globalne je `bias` maly, ale u vsech outcome je vyznamny trend chyby:

- `mmt8`: trend_slope 0.917, p < 1e-26
- `mmt10`: trend_slope 0.847, p < 1e-22
- `fi2`: trend_slope 0.841, p < 1e-19
- `fi3`: trend_slope 0.872, p < 1e-21
- `haq`: trend_slope 1.192, p < 1e-35

Co to znamena prakticky:

- Model ma proporcionalni chybu (stahovani ke stredu).
- Pri nizsich hodnotach casto nadhodnocuje, pri vyssich podhodnocuje.
- To je zasadni limit pro plnou substituci klinickeho mereni na urovni jednotlivce.

## 7) Klicovy test mimo strop/dno: vztah ig_gradient -> outcome

Asociace v non-boundary casti (`assoc_non_ceiling`, smer `ig_gradient`):

- `mmt8`: tau 0.229, rho 0.328, delta P90-P10 = +5.93 (95% CI 2.88; 9.30)
- `mmt10`: tau 0.263, rho 0.388, delta P90-P10 = +9.51 (95% CI 5.73; 13.46)
- `fi2`: tau 0.345, rho 0.502, delta P90-P10 = +30.26 (95% CI 20.10; 41.55)
- `fi3`: tau 0.270, rho 0.401, delta P90-P10 = +25.11 (95% CI 16.90; 34.52)
- `haq`: tau -0.154, rho -0.224, delta P90-P10 = -0.393 (95% CI -0.639; -0.100)

Interpretace:

- Mimo strop/dno je vztah statisticky vyznamny u vsech 5 skal.
- Smer efektu je klinicky konzistentni:
  - vyssi (mene negativni) `ig_gradient` -> vyssi svalove/vykonove skaly,
  - vyssi `ig_gradient` -> nizsi `haq` (lepsi funkcni stav).

## 8) M0 vs M6 a interakce

Model s interakci (`ig_gradient * exam_order`) ukazuje:

- Signifikantni interakce pouze u `mmt10` (p = 0.048).
- U ostatnich outcome interakce nevysla signifikantne.

Prakticky dopad:

- Efekt `ig_gradient` muze byt u `mmt10` odlisny mezi M0 a M6.
- U ostatnich skal neni dost dukazu, ze by se efekt mezi M0/M6 vyrazne menil.

## 9) Smer promenne (`ig_gradient` vs `ig_grad_pos`)

`ig_direction_choice`:

- `ig_gradient`: 5/5 outcome konzistentni smer.
- Inverze (`ig_grad_pos = -ig_gradient`) konzistenci nezlepsila.

Doporuceni z analyzy: ponechat puvodni smer `ig_gradient`.

## 10) Benchmark: prinos proti jednoduche baseline

Model byl porovnan s baseline (globalni prumer + prumer podle exam):

- Zlepseni MAE i RMSE je kladne u vsech outcome (global i non-boundary).
- Nejvetsi zisky globalne:
  - `fi2`: MAE +2.13, RMSE +2.72
  - `fi3`: MAE +1.51, RMSE +2.64
  - `mmt10`: MAE +1.40, RMSE +1.40
  - `mmt8`: MAE +0.94, RMSE +0.97
- `haq`: zisk proti baseline je maly (MAE +0.036, RMSE +0.038).

Interpretace:

- `ig_gradient` prinasi realny signal nad trivialni predikci prumerem.
- To ale samo o sobe neznamena, ze je presnost dostatecna pro klinickou nahradu skaly.

## 11) Informativeness QC pro ig_gradient

- Celkem: mean -2.771, SD 0.185, IQR 0.208, rozsah -3.594 az -2.207.
- Podle navstev je variabilita podobna, lehce rozdilna podle exam.
- Korelace s dalsimi aktivnimi metrikami jsou casto stredni az vysoke
  (napr. s vyssimi percentily ENMO a nekterymi MVPA metrikami).

Zaver QC: informacni obsah `ig_gradient` je **stredni** (`moderate`).

## 12) Sensitivity analyza: mixed modely

V non-boundary casti byly smery efektu konzistentni i v mixed modelech
(`(1|patient_id)`), tedy zavery nejsou dane jen ignorovanim opakovanych
mereni.

## 13) Klinicky zaver: lze nahradit skaly pouze ig_gradient?

### Prakticka odpoved

- **Plna nahrada klinickych skal samotnym `ig_gradient` se aktualne nedoporucuje.**

### Podle skal

- `mmt8`, `mmt10`:
  - nejblize praktickemu vyuziti,
  - vhodne jako doplnkovy digitalni marker trendu,
  - ne jako samostatna nahrada klinickeho hodnoceni.
- `fi2`, `fi3`:
  - mimo strop je vztah pritomen, ale absolutni chyba predikce je vysoka,
  - boundary cast je slaba,
  - nahrada neni vhodna.
- `haq`:
  - globalni vysvetleni je slabe a trend BA-like je vyrazny,
  - nahrada neni vhodna.

### Jednovetne shrnuti

`ig_gradient` je uzitecny doplnkovy biomarker pro orientacni monitoring
(starsi signal u `mmt8/mmt10`), ale v teto kohorte nema parametry pro
bezpecnou a plnou substituci klinickych skal (`mmt8`, `mmt10`, `fi2`,
`fi3`, `haq`).
