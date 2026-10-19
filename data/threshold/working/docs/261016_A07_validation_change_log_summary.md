# MYO-II / Paper A — Revize A07: validace, change log a vědecké shrnutí

- **Datum:** 2026-10-16
- **Skripty:** `07_REV_A_B_profile_concentration_BINS.R`, `08_REV_C_D_E3_shape_BOTH.R`, `09_REV_E_subtypes_clinical.R`, `10_REV_VALIDATION_BOTH.R`, helper `R/06_helpers_revision.R`
- **Report:** `reports/261016_thresholds_report_v2.Rmd` / `.html`. Původní `261013_thresholds_report.*` zůstává beze změny.
- **Běh:** Windows R 4.6.1, 31 workerů, B = 2000, seed 20261013, plasmode 500 × 1000. Master md5 `a4821db9803a46c31eb4dc6639add5b9`. Klinický soubor `klinicka_data_12102026.xlsx` (md5 je v `output/logs/261016_09_REV_E_subtypes_clinical_sessionInfo.txt`).
- **Rozhodnutí a jejich zdůvodnění:** `261016_A07_critical_review.md`. Oponentura: `261016_A07_independent_review.md`, vypořádání v §5.

---

## 1. Validační zpráva (skript 10)

Výsledek: **0 FAIL, 1 WARN** (17 kontrol, po rozšíření podle oponentury).

| ID | Kontrola | Výsledek |
|---|---|---|
| R1 | Exkluzivní biny rekonstruované z kumulativních prahů (D8–D12) | PASS, max. chyba < 1e-8 |
| R2 | Šířky binů 10 / 20 / 25 mg, faktor 10/šířka, normalizovaný kvantil = surový × faktor | PASS (23 binů: 10 × 10 mg, 9 × 20 mg, 4 × 25 mg) |
| R3 | Nemožné / nekonečné hodnoty v nových výstupech | PASS. Jediná NA jsou dokumentované nedefinované buňky původního estimandu B (5 buněk). |
| R4 | Top-10 % podíly přepočítané nezávisle (base R přímo z `bin_long`) | PASS, 69 buněk, rozdíl ≤ 1e-14 |
| R4b | Původní C3 reprodukováno | PASS |
| R5 | Uzavření clr (Σclr = 0), součet 1440 min, nuly | PASS. Nula se vyskytla jen ve 4 z 217 záznamů, max. 3,2 % v jedné části. |
| R5b | Invariance reziduí E3 vůči normalizaci na 10 mg a centrování clr; stejná báze splinu | PASS (rozdíl < 3e-16) |
| R6 | Simultánní pásy E3 přepočítané z 2000 uložených replikací | PASS (rozdíl 0) |
| R6b | Reprodukovatelnost (stejný seed; 1 vs. N workerů) | PASS |
| R6c | Maticový bootstrap (páry) vs. primární (všichni pacienti) | PASS. Odhady jsou totožné, šíře pásu se liší o ≤ 4 %. |
| R7 | Leave-one-patient-out | **WARN.** V binech ≥ 240 mg posune vyřazení jednoho pacienta podíl top 10 % až o 10,9 p.b. (limit 10). Vysoké hodnoty koncentrace při M3 proto závisí na jednotlivcích, což report uvádí. Koncentrace T100/T200 se posune max. o 6,1 p.b., rozdíl sklonu o 0,51 SE, kontrasty podtypů o 0,52 SE. Kontrast IMNM−DM pro ΔT100 nikdy nezmění znaménko. |
| R7b | Fity lmer modelů podtypů (vzorek s ph_VAS) | PASS. Žádný singulární fit ani konvergenční hlášení, ICC 0,56–0,73. |
| R8 | Velikosti vzorků (srovnání s původním výstupem 01) | PASS. Páry M3 = 47, M6 = 52, M18 = 39; 71 pacientů (hash je jen popisný). Klinika: napojeno 71/71, kódy podtypů souhlasí 71/71. |
| R9 | Primární pozorované odhady E1–E4 (všechny bloky) | PASS, rozdíl 0 |
| R9b | Původní surový profil, C3, C1 | PASS (92/92, 69, rozdíl 0) |
| R10 | Data figur = tabulky | PASS (A1, B1, C1, D1, E1) |
| R10b | 9 tvrzení v reportu ověřených automaticky | PASS, 9/9 |

**Upozornění, která nejsou chybou:**

- Percentilové bootstrapové intervaly pro statistiky typu „top-k“ jsou přibližné, protože duplikace největších zisků v replikaci zvyšuje koncentraci. V reportu je to uvedeno.
- Maticový bootstrap (C, D) převzorkuje kompletní páry dané návštěvy, primární bootstrap všech 71 pacientů. Pro E3 se šíře pásů liší o ≤ 4 % (R6c).

---

## 2. Change log (co se změnilo)

| Položka | Typ | Původně | Revize |
|---|---|---|---|
| Popisný profil BINS (`01_descriptive_profile_bins`) | figura nahrazena v reportu v2 | surové min/den v binu, lokální vzestup 120–140 mg | min/den na 10 mg (`07_A1`) + srovnání (`07_A2`). Profil je při všech návštěvách monotónně klesající. |
| F6 heterogenita | nový estimand, nová figura (`07_B1`) | „podíl top 10 % na celkovém zisku“ = podíl na *čistém* zisku; nad 200 mg 46–75 %; nízké biny NA nebo až 1898 % | podíl top 10 % na *kladných* ziscích: 32–65 % (minimum při stejných ziscích 13–30 %), s 95% CI; index net/gross; panely nízké a vysoké intenzity. Původní estimand je zachován (`07_B2`) a nestabilní buňky jsou označené. |
| Klíčové sdělení 5, H3 | text | „high-intensity zisk nese menšina silných responderů“ | „směrem široká, velikostí nerovnoměrná změna, v podobné míře v celém spektru; vyšší koncentrace při M3 nad ~240 mg není formálně podložená“ |
| F3 (E3) | figura překreslena (`08_C1`), text | šedá křivka, stručný popis | viditelná hladká křivka pro každou návštěvu, výklad panelů A/B a kritéria. **Čísla beze změny.** |
| Robustnost E3 | nové analýzy | invariance jen vůči nulám, ≥ 5 dnům a UCTD | + df 1–5, vážený fit, mřížka 20 mg, biny do 400 mg, plasmode pro volbu df (`08_C2`, `08_C3`). Verdikt „bez lokální odchylky“ platí pro df 3, 4 i 5, vážený fit i obě mřížky. df 3 zůstává jako předem zvolené, jeho optimalita prokázaná není. |
| Tvar mezi návštěvami | nové analýzy (`08_D1`) | jen C1 (M6 − M3 po binech) | párové testy sklonu a hladké křivky (3 df) pro M3/M6/M18 s Holmovou korekcí; společný vzorek; tabulka sezóny a časování. Nic nepřežije Holmovu korekci (nejnižší p = 0,037 → 0,11). |
| Klíčové sdělení 4, H4 | text | „formal M3-vs-M6 comparison not significant“ | odkaz na nové testy, H4 zůstává hypotézou |
| Podtypy | nové analýzy (`09_E1`, F7b) | jen deskripce (F7, tabulka) | modely změny ΔT100, sklonu a ΔIG s adjustací na ph_VAS (hlavní) a MMT8 (citlivost); konstrukt, kolinearita, dostupnost. F7 beze změny. |
| Klíčové sdělení 6, H5 | text | „no subtype conclusion is possible“ | IMNM má menší objemový zisk (ΔT100 −20 min/den po adjustaci na ph_VAS). Ve vzorku s MMT8 se odhad adjustací téměř nemění (−18,8 → −17,5), jen se rozšíří CI. MMT8 s ΔT100 nesouvisí, takže svalová slabost rozdíl nevysvětluje. Rozdíl v tvaru změny mezi podtypy prokázán není (široké CI). |
| Limitace, další otázky, vědecké shrnutí | text | — | doplněno |
| Ukládání RData (`save_public()` v `R/00_config.R`, skripty 01–10) | ochrana dat | ggplot objekty v RData nesly prostředí skriptu (`plot_env`, quosury `aes()`, interní prostředí vrstev) včetně tabulek po pacientech; týkalo se `qc_objects` (a27b505) a `rev_AB/CD/E` | před uložením se objekt zkopíruje a z prostředí se odstraní vazby s řádky `patient_id`; zápis selže, pokud nějaké zůstanou. Existující RData přeuloženy (figury ověřeny, ostatní obsah beze změny), commit a27b505 přepsán. |

**Beze změny:** všechny primární estimandy E1–E4, jejich odhady, intervaly, pásy a confidence sets; populace (71 pacientů, páry 47 / 52 / 39); mřížka prahů a binů; modely; citlivostní analýzy S1–S10; validace V1–V7. Ověřeno v R9.

---

## 3. Vědecké shrnutí

**Co revize potvrzuje**

1. Denní profil intenzity ENMO (na 10 mg) s intenzitou monotónně klesá. Lokální vrchol 120–140 mg byl artefakt šířky binu v popisné figuře. Primární E1 byl normalizovaný už dříve.
2. Zlepšení se projevuje jako plynulé **naklopení** celé 24h kompozice k vyšším intenzitám, bez lokálního pásma. Závěr je robustní vůči reprezentaci šířky binů (algebraicky invariantní), vůči mřížce i vůči každé kalibrované specifikaci hladkého nullu.
3. Směr změny je široký: nad ~70 mg získává 51–83 % pacientů. **Velikost** zisku je ale nerovnoměrná, v podobné míře ve středních i vysokých intenzitách. Top 10 % pacientů nese zhruba třetinu až dvě třetiny všech kladných zisků.

**Co zůstává nejisté**

- Zda se tvar změny liší mezi časnou (M3) a pozdější (M6/M18) kontrolou. Vizuálně ano, statisticky ne. Návštěva, sezóna a účast jsou navíc strukturálně svázané.
- Malé lokální odchylky (≈ +10 %), zejména při M18.
- Rozdíly podtypů. Signál „IMNM menší objem“ je exploratorní a odhad s adjustací na MMT8 je nepřesný. Léčba ani trvání nemoci nejsou modelované.
- M18 stojí na předpokladu MAR.

**Na co navázat v klinicky orientované studii**

- Ukotvit změnu spektra změnou MMT8, FI2, HAQ a PRO (klinicky významná změna).
- Předem specifikovaně otestovat, zda IMNM získává menší objem aktivity než DM při stejné výchozí síle.
- Reliabilita (ICC/SEM/SDC) z denních dat GGIR. `d11_part2daysummary` v `data/260123_backup.RData` je třeba ověřit.
- Sbírat důvody vynechaných návštěv a plánovat kontroly tak, aby se oddělila sezóna od času od M0.

---

## 4. Reprodukce

```r
# v RStudiu ve složce working/ (nebo source("run_all.R"), který obsahuje i 07–10)
for (s in c("07_REV_A_B_profile_concentration_BINS.R", "08_REV_C_D_E3_shape_BOTH.R",
            "09_REV_E_subtypes_clinical.R", "10_REV_VALIDATION_BOTH.R")) source(s, local = new.env())
rmarkdown::render("reports/261016_thresholds_report_v2.Rmd")
```

Doba běhu (31 workerů): 07 ≈ 13 min, 08 ≈ 1 min, 09 ≈ 17 min, 10 ≈ 15 min. Bootstrap se ukládá do checkpointů v `output/RData/boot_chunks/` (tagy `rev_B`, `rev_E`).

---

## 5. Vypořádání nezávislé oponentury

Oponent: Claude Fable 5.1 (`261016_A07_independent_review.md`). Verdikt **MINOR REVISION**. Kritické nálezy žádné, výpočty jsou správné a primární výsledky se reprodukují. Všechny připomínky jsou přijaté.

| ID | Závažnost | Připomínka (zkráceně) | Vypořádání |
|---|---|---|---|
| F1 | major | „Oslabení“ IMNM−DM po adjustaci na MMT8 je srovnání mezi různými vzorky. Ve stejném vzorku −18,8 → −17,5, MMT8 není asociované (p = 0,71). | **Přijato.** Report (KM6, sekce E, H5, „How to use“, change log) a tento dokument přepsány: MMT8 rozdíl nevysvětluje, jen snižuje přesnost. |
| F2 | major | Reference „rovnoměrného sdílení“ je matematické minimum. „Nad referencí všude“ platí triviálně. „M3 nejvíc koncentrované“ není testováno. Nejvyšší podíly jsou ve 20–40 mg. Poměr 2,0–5,1×, ne 2–4×. Vliv šikmosti výchozích hodnot. | **Přijato.** Reference v reportu označena jako minimum. Interpretace změněna na „podobná nerovnoměrnost v celém spektru“; rozdíl M3 vs. ostatní označen jako netestovaný; doplněny nízké biny, poměr 2,0–5,1× a poznámka o výchozích hodnotách. |
| F3 | major | Zdůvodnění df 3 je kruhové (true null = fit df 3). Při df-5 null nejsou df 3 ani df 4 kalibrované (M6). df 4 má větší sílu pro úzké hrby. | **Přijato.** df 3 je v reportu zdůvodněno jako předem zvolené s rozumnou rovnováhou, ne jako optimální. Závěr opřen o shodu df 3–5. Upravena i limitace. |
| F4 | minor | Konstantní posun v D je identifikovaný. | **Přijato.** Komentář v helperu, report a kritické review opraveny. |
| F5 | minor | p-hodnota sklonu z kratšího z obou chvostů je mírně antikonzervativní. | **Přijato.** Kód změněn na `mean(abs(T* − t) ≥ |t|)`, skript 08 znovu spuštěn. Závěry beze změny. |
| F6 | minor | D vyžaduje M0, ačkoli se M0 v rozdílu odečte. | **Přijato jako dokumentace.** Důvod (konzistence s E3 a C1) uveden v reportu. |
| F7 | minor | Holmovo p ve společném vzorku (0,071) neuvedeno. | **Přijato.** Doplněno. |
| F8 | minor | LRT s ML na ~47 pacientech jsou antikonzervativní. | **Přijato jako limitace** (text v reportu). Parametrický bootstrap LRT se nedělal, závěry stojí na bootstrapových intervalech. |
| F9 | minor | Koeficient závažnosti v E_c je efekt v M3. | **Přijato jako dokumentace** (report, komentář v kódu). Popisek ve výstupu je beze změny kvůli konzistenci s checkpointy bootstrapu. |
| F10 | minor | Rozdíl ASyS−IMNM pro ΔIG mizí už omezením vzorku. | **Přijato.** Text opraven. |
| F11 | minor | „ph_VAS mostly extramuscular“ je nadsázka; ρ MMT8–FI2 = 0,82. | **Přijato.** „substantially driven“; obě ρ uvedena; ρ se počítá v textu. |
| F12 | minor | „Tvar se mezi podtypy neliší“ zní jako důkaz absence. | **Přijato.** „no evidence of a difference, wide intervals“. |
| F13 | minor | Hranice 0,20 je arbitrární; nad 200 mg je net/gross 0,49–0,85, takže původní číslo nebylo „nafouknuté“. | **Přijato.** Formulace zmírněna; hranice označena jako pragmatická. |
| F14 | minor | Text issue logu odporuje revidované H3. | **Přijato.** Opraveno. |
| F15 | minor | Mezery validace: hash je tautologický; LOO nepokrývá vysoké biny; chybí kontrola fitů lmer. | **Přijato.** R8 přejmenováno (srovnání počtů s původním výstupem 01, hash jen popisný). Doplněno LOO pro biny ≥ 240 mg (R7) a kontrola singulárních fitů a konvergence (R7b). |
| F16 | minor | HTML neodpovídá Rmd. | **Přijato.** Report znovu vyrenderován po všech opravách. |
| F17 | minor | Reprodukce uvádí 01 → 06. | **Přijato.** 01 → 10 + klinický soubor. |
| F18 | minor | 115 vs. 126–127 pozorování. | **Přijato.** Opraveno v kritickém review. |
