# MYO-II / Paper A — Fáze 1–2: (A) audit, (B) plán, (C) další otázky

- **Datum:** 2026-10-13
- **Autor (první role):** Claude Opus 5.5 (Claude Code). Pracoval v roli biostatistika a revmatologa.
- **Stav:** návrh před nezávislou oponenturou (fáze 3). Finální analytický plán **není uzamčen**.

## 0. Použité vstupy

| Soubor | Použití |
|---|---|
| `data/threshold/MYOII_thresholds_agent_handoff_prompt.md` | Zadání |
| `data/threshold/RULE_rmd_report_style.md` | Normativní styl `.Rmd` (fáze 4) |
| `data/threshold/MYOII_ACC_OUTCOMES_MASTER_V1.rds` | Kanonická data (P06 V3, vytvořeno 2026-10-12 14:08) |
| `kumulativni_cas_nad_prahy/A04_FINAL_CUMULATIVE_THRESHOLD_ANALYSIS_v2_3.R` + výstupy | Audit |
| `biny/A07_FINAL_EXCLUSIVE_BIN_ANALYSIS_v2_3.R` + výstupy | Audit (prohledáno cíleně, nečteno celé) |
| `R01_R04/PaperA_R00_helpers.R`, `PaperA_R01…R04_*.R` + výstupy | Audit |
| `R01_R04/MYOII_R05_R06_revised_figures_and_legends/{R00,R05,R06}` + výstupy | Audit (revidovaná verze) |
| `R01_R04/R00_MYOII_common.R`, `R05_*.R`, `R06_*.R` | Porovnání s revidovanou verzí |
| `MYOII_PaperA_A05-A07_session_methods_results.md` | Kontext (starší dataset, překonáno) |
| `A04_A07_PaperA_statisticke_zavery_a_plan_robustnosti.docx` | Kontext |
| `R01_R04/PaperA_R06_SLIDING_WINDOWS_BINS/MYOII_PaperA_integrated_review_R01_R06.docx` | Kontext |
| `README*.txt`, `R05_R06_PUBLICATION_FIGURE_LEGENDS.txt` | Kontext |

Složka `older_dataset/` je jen historický kontext. Pracuje nad necleaned daty a nebyla auditována do detailu.

**Běhové prostředí:**
- Kontroly probíhaly v R 4.6.1 (WSL/Linux).
- Projektová `renv` knihovna pro Linux **není nainstalovaná**. Mimo renv má `rlang` nekompatibilní binárku, takže tidyverse ani lme4 nejdou načíst.
- Všechny kontroly níže proto proběhly v **base R** přímo nad masterem.
- LMM výsledky A04/A07/R04 **nebyly reprodukovány**: NOT EXECUTED.

---

## A. Audit dat a skriptů

### A1. Estimandy — co lze s daty hodnotit

| Dimenze | Definice | Hodnotitelné? |
|---|---|---|
| (i) Skupinová longitudinální změna | Průměrná změna min/den v intervalu intenzity Mx − M0 v populaci pacientů s IIM | **Ano.** Primární estimand. |
| (ii) Přesnost odhadu | SE / šířka CI estimandu (i) | **Ano.** Pacientský cluster bootstrap a modelové SE. |
| (iii) Zastoupení expozice | Podíl pacientů s expozicí nad popsaným limitem | **Ano**, deskriptivně a s citlivostí na limit. |
| (iv) Individuální reliabilita / responsiveness | Test-retest ICC, SEM, MCID | **Ne.** Master nemá denní data ani opakované měření v rámci návštěvy, takže ICC/SEM nelze spočítat. Lze jen distribuci individuálních Δ (heterogenitu), ne reliabilitu. |
| (v) Prediktivní / klinická informační hodnota | Vztah ke klinickému stavu | **Ne** (mimo scope). Budoucí otázka. |

**Navržený primární estimand:** populační průměrná změna času v exkluzivním binu b mezi M0 a Mx, Mx ∈ {M3, M6, M18}.
- Je vyjádřená na **absolutní škále** (min/den) a na **relativní škále** (poměr průměrů Mx/M0).
- Odhaduje se za předpokladu MAR z LMM se všemi dostupnými návštěvami.
- SD-standardizovaný efekt je **sekundární** harmonizační ukazatel (viz A4-K1).

### A2. Datové QC kontroly (base R, master V1)

| # | Kontrola | Výsledek | Stav |
|---|---|---|---|
| D1 | `record_level`: unikátnost patient × visit | 251 záznamů, 0 duplicit. `qc$duplicates` je ve všech zdrojích bez duplicit. | PASS |
| D2 | QC eligibility | 251/251 `qc_analysis_eligible = TRUE`, 0 hard-fail, 0 calibration review (cal.error.end max 0,0059 g < 0,01 g) | PASS |
| D3 | Hlavní vzorek | 217 měření / 71 pacientů: M0 63, M3 52, M6 57, M18 45. Shoda s docx. | PASS |
| D4 | Páry s M0 | M3 47, M6 52, M18 39. Shoda s integrated review. | PASS |
| D5 | Pacienti bez M0 | 8 pacientů (11 měření) nemá M0. Do párových analýz nevstupují, do LMM ano (přispívají k odhadu přes náhodný intercept). | WARN — popsat v flow diagramu |
| D6 | Mřížka prahů | 24 prahů 20…400 mg dle zadání, 24 řádků/záznam, 0 NA | PASS |
| D7 | Biny | 23 ohraničených binů 20–400 mg (šířky 10/20/25 mg). **Otevřený bin >400 mg v masteru není.** Komentáře v R00/R03 o „top open bin“ jsou neaktuální. | PASS (WARN: zastaralé komentáře) |
| D8 | Monotonicita `minutes_above` | 0 porušení | PASS |
| D9 | Identita `bin = T_lo − T_hi` | max odchylka 0 | PASS |
| D10 | Σ biny = T20 − T400 | max odchylka 0 | PASS |
| D11 | `minutes_per_10mg = minutes_bin × 10 / šířka` | max odchylka 0 | PASS |
| D12 | Nenegativita binů | 0 záporných | PASS |
| D13 | Přesné nuly v binech | ≤ 0,5 % do 260 mg, max 6,5 % (375–400 mg). Malé nenulové hodnoty převažují: ≤1 min/den má 48 % (240–260), 65 % (260–280) a 96 % (375–400). | INFO — hurdle model **není** podložen, problém je škála, ne nuly |
| D14 | Validní dny | 3 d: 1, 4 d: 9, 5 d: 14, 6 d: 177, 7 d: 50 | WARN — chybí doložené pravidlo minimálního počtu dnů; citlivost ≥ 5 d |
| D15 | `acc_date` vs `start_time` | U 64/251 záznamů je `acc_date` o 1 den dřív než lokální datum `start_time` (start 00:15 Europe/Prague, převod přes UTC). `acc_month_num` i `season` jsou správně (0 nesouladů). | WARN — kosmetické, sezónu neovlivňuje |
| D16 | Subtyp | Kódy 0–5 a 1× „UCTD/SSc?“ | **WARN — vyžaduje rozhodnutí:** legenda kódů a zda pacient s „UCTD/SSc?“ patří do IIM populace |
| D17 | Klinické napojení | 3 záznamy bez kliniky (1× M30, 2× M18). Na akcelerometrické analýzy bez kovariát nemá vliv. | PASS (INFO) |
| D18 | Rozestupy návštěv od M0 | M3: medián 3,2 m (2,0–6,0), M6: 6,6 m (4,3–9,4), M18: 19,5 m (16,5–21,8) | WARN — M3 a M6 se časově překrývají (6,0 vs 4,3 m). Citlivost na skutečný čas. |
| D19 | Sezóna × návštěva | M0: zima 6/63 (10 %), M6: zima 16/57 (28 %) | WARN — částečné confounding návštěva × sezóna, ale každá buňka ≥ 6 |
| D20 | Výběrový drop-out (hrubě) | Medián baseline ENMO u pacientů s M18 20,7 mg (n 39), bez M18 20,5 mg (n 24) | INFO — hrubě bez rozdílu; formální analýza viz C2 |

### A3. Reprodukce dosavadních výsledků

| Výstup | Reprodukováno? | Poznámka |
|---|---|---|
| R05 párová změna (pooled SD) per bin | **Ano**, base R | Např. M18 100–110: 0,4334 (xlsx 0,43345). Peaky M3 260–280 0,585, M6 140–160 0,532, M18 160–180 0,479. |
| R06 peaky oken 60/100/140 mg (Δ SD + coverage) | **Ano, všech 9 přesně** | Viz tabulka v integrated review |
| R01/R02/R03 | Ne (bez tidyverse) | Logika R01 je ekvivalentní ověřené R05. |
| A04/A07 LMM, R04 | **NOT EXECUTED** | Vyžaduje lme4 |
| `season_anchor_test` (A04 i A07) | Prázdný list příčinně vysvětlen | Viz A4-K3 |

### A4. Nálezy v kódu

Závažnost: **K** = kritická pro interpretaci / reprodukovatelnost, **V** = vysoká, **S** = střední, **N** = nízká.

| ID | Závažnost | Místo | Nález | Navržená oprava |
|---|---|---|---|---|
| K1 | **K (vědecká)** | R01/R03/R05/R06 (`sd_pool` = SD přes všechny záznamy), A04/A07 (z-skóre přes M0–M18) | **SD-standardizace řídí polohu „peaku“.** Pooled SD binu klesá z 24,9 min (20–30) na 0,40 min (375–400). Peak M3 v 260–280 mg (0,585 SD) odpovídá **+0,68 min/den**, kdežto 120–140 mg (0,41 SD) odpovídá +3,5 min/den. Vysoké SD-efekty ve vysokých binech jsou tedy změny pod 1 minutu. Nadhodnocení „vysokointenzitní odpovědi M3“ je převážně artefakt škály. | Primárně reportovat absolutní Δ min/den a relativní změnu (poměr průměrů). SD-efekt jen sekundárně a vždy vedle absolutní škály. |
| K2 | **K (reprodukovatelnost)** | Revidovaný `MYOII_R05_R06_revised…/R00_MYOII_common.R:6` | `INPUT_RDS` natvrdo `R:/…`; proměnná `MYOII_MASTER_RDS` se **ignoruje**, ač ji README i chybová hláška doporučují. Starší `R01_R04/R00_MYOII_common.R` ji respektuje (regrese). Soubor má CRLF + UTF-8, ostatní ASCII. | V nových skriptech jedna konfigurace cesty (`here`/projektová relativní + env override). |
| K3 | **V** | A04 ř. 194–196 (a obdobně A07) | `comparison$\`Chi Df\`` neexistuje v současném výstupu `anova.merMod` (sloupec je `Df`). `as.numeric(NULL[2])` vrací délku 0, takže `tibble()` vytvoří **0 řádků**. Proto je `season_anchor_test` prázdný (ověřeno: list obsahuje jen hlavičku 10 sloupců). | `comparison$Df[2]`. Test sezónnosti zatím **neexistuje**. |
| V1 | V | R03 ř. 53–56, R06 ř. 176–177 | `clusterSetRNGStream` + `parLapply`: výsledky závisí na **počtu workerů** (jiné rozdělení replikací do streamů). Na jiném PC s jiným počtem jader dají jiné CI. | Seed per replikace (L'Ecuyer streams indexované replikací), nebo `furrr` s `furrr_options(seed = SEED)`. Ten je reprodukovatelný nezávisle na počtu workerů. |
| V2 | V | R05 ř. 26–41 | „Patient bootstrap“ R05 převzorkovává **nezávisle pro každý visit × bin**. Pointwise CI to nezkreslí, ale nejde o cluster bootstrap přes spektrum. Nelze z něj odvodit simultánní pás, stabilitu tvaru ani selekční četnost peaku. Legenda R05 tvrdí „patient-level bootstrap resampling“. | Jeden draw pacientů pro všechny biny × návštěvy (jak dělá R06). Doplnit simultánní (max-t) pásy. |
| V3 | V | R03/R05/R06 | Standardizační SD je v bootstrapu fixní z originálu. Handoff žádá přepočet dle pravidla. Navíc SD počítá ze **všech** záznamů, i nepárových. | Rozhodnout předem. Doporučení: absolutní škála primárně, takže problém zmizí. Pro SD-efekt přepočet v replikaci. |
| V4 | V | R06 | Okna lze tvořit jen tam, kde `start + šířka` padne na hranu nepravidelné mřížky: 14 / 11 / 9 oken. Pokrytí spektra je nerovnoměrné (např. šířka 60 nemá okna začínající v 300 mg). Okna se překrývají, nepřináší nový estimand a selekce maxima vnáší winner's curse. | Sliding windows **omezit na supplementary**, nebo vypustit. |
| S1 | S | A04/A07 vs R01/R05 | Různá analytická populace a škála: A04/A07 = LMM, všech 217 měření, z-skóre (mean/SD M0–M18). R0x = complete pairs, Δ/SD všech záznamů. Čísla nejsou přímo srovnatelná (peaky A07 M3 0,565 vs R05 0,585). | V reportu explicitně oddělit, nebo sjednotit škálu. |
| S2 | S | A04/A07/R04 | Wald CI z REML LMM bez ověření reziduí. Při silné heteroskedasticitě (rozptyl klesá s intenzitou, ale model je per bin) jsou CI potenciálně nepřesné. Rezidua nejsou nikde uložena. | Diagnostika reziduí + cluster bootstrap CI pro LMM. |
| S3 | S | R05 vs R06 | „SE“ v R05 je analytické `sd/√n`, v R06 bootstrapové. Ve stejném typu figury jsou dvě definice. | Jednotně bootstrap SE. |
| S4 | S | R02/R05/R06 | Coverage = „> 0,25 min/den/10 mg v M0 **nebo** follow-up“ (denominátor = páry). U kumulativních prahů je bez normalizace na šířku. Definice „either“ nadhodnocuje zastoupení oproti „both“ nebo per-visit. | Hlavní coverage per visit (podíl měření nad limitem), citlivost na limity 0 / 0,1 / 0,25 / 0,5 / 1. |
| S5 | S | Celkově | Multiplicita: 23 binů × 3 návštěvy, pointwise CI. Žádná simultánní korekce. | Simultánní pásy (max-t bootstrap) jako primární nejistota tvaru. |
| N1 | N | R01 ř. 11 | `pct_nochange = mean(delta == 0)` na float minutách | Tolerance |
| N2 | N | R00 (oba) | Zastaralý komentář o otevřeném binu > 400 mg | Odstranit |
| N3 | N | Grafy | `theme_minimal`, ne projektové `theme_bw`. Linetype M0 = `dotdash` (zadání: „vlastní odlišná čára“, OK). | Sjednotit dle AGENTS.md. |

**Kritické body pro bránu fáze 1:**
- K1 je vědecká a řeší se volbou estimandu v plánu B, ne opravou kódu.
- K2 a K3 jsou kódové chyby, které se opraví v nové implementaci.
- Dosavadní **interpretace** „M3 = odpověď ve vysokých intenzitách“ (H1 v integrated review) **neprošla** kontrolou K1. Dokud nebude přepočtena na absolutní / relativní škále, je pozastavena.

---

## B. Doporučený analytický plán (návrh k oponentuře)

### B1. Hlavní otázka a hierarchie

**Otázka:** Které části spektra ENMO (20–400 mg) nesou reprodukovatelnou longitudinální změnu běžné pohybové aktivity napříč pacienty a návštěvami?

| Úroveň | Strategie | Role | Zdůvodnění |
|---|---|---|---|
| Primární | **Full-spectrum exclusive bins** (23 nativních binů) | Inference | Nepřekrývající se expozice, žádný apriorní výběr pásma, přímá lokalizace |
| Komplementární | Kumulativní prahy (24) | Deskripce + srovnání s literaturou | Překryv expozic, ale srovnatelné s klasickými cut-pointy |
| Supplementary / vypustit | Sliding windows | Jen citlivost na hrubost mřížky | V4: nepravidelná mřížka, žádný nový estimand, riziko post-hoc selekce |
| Nedoporučeno nyní | Nová extrakce mřížky po 10 mg | — | Vyžaduje rerun GGIR nad raw `.bin` soubory (není v masteru). Přidaná hodnota je nejistá: biny 10–25 mg už jsou nad 200 mg řídké, jemnější mřížka sníží expozici na bin. Navrženo jako budoucí krok jen při zjištění ostrých hranic. |

### B2. Modely (per bin, každý estimand explicitně)

| Model | Estimand | Použití |
|---|---|---|
| **M-A: LMM raw minut** `minutes ~ visit + (1 \| patient)`, všechna měření M0–M18 | Δ min/den (MAR) | **Primární absolutní škála** |
| **M-B: Gamma GLMM s log-linkem** (`glmmTMB`), nebo LMM na `log(minutes + c)` s c předem pevným a citlivostí na c | Poměr průměrů (resp. geometrických průměrů) Mx/M0 | **Primární relativní škála**. Řeší pokles rozptylu s intenzitou (K1). Přesné nuly jsou ≤ 6,5 % (D13), proto hurdle model nedává smysl. U Gamma modelu se nuly ošetří malým offsetem a citlivostí. |
| M-C: complete-pair Δ (R01/R05 logika) | Párová průměrná Δ | Citlivost na missingness |
| M-D: SD-standardizovaný efekt | Δ / SD(M0) | Sekundární, jen pro srovnání s A04/A07 |
| M-E: + sezóna (faktor) / harmonické sin–cos | Δ adjustovaná | Citlivost; opravený anchor test (K3) |
| M-F: skutečný čas od M0 (měsíce) místo nominální návštěvy | Trajektorie | Citlivost na D18 |

**Nejistota:**
- Cluster bootstrap přes pacienty, B = 2000. Jeden draw platí pro celé spektrum a všechny návštěvy; duplicitní pacient dostane nové ID pro LMM.
- Pointwise percentile CI a **simultánní max-t pás** přes 23 binů pro každou návštěvu.
- Wald CI se ukládá pro srovnání.

**Reprodukovatelnost:** `furrr` + `future::plan(multisession)`, `furrr_options(seed = SEED)` (řeší V1).
- Počet workerů = `parallel::detectCores(logical = TRUE) - 1` (fallback 1), omezený `getOption("myo_ncores")`.
- **Konflikt:** handoff předepisuje „cores − 1“, projektový AGENTS.md `myo_ncores`. Navrhuju použít minimum z obou.
- Ukládá se B, seed, počet workerů a selhání fitů.

### B3. Dimenze informační hodnoty (předem definované, bez kompozitního skóre)

1. **Magnitude:** Δ min/den a poměr průměrů s CI.
2. **Precision:** bootstrap SE a šířka CI, a to na absolutní i relativní škále.
3. **Coverage:** podíl měření s expozicí > limit, per visit. Limity 0 / 0,1 / 0,25 / 0,5 / 1 min/den/10 mg; 0,25 je jen ilustrační.
4. **Robustness:** shoda směru a tvaru napříč M-A…M-F, complete-pair vs LMM a bootstrapová stabilita *souvislé oblasti* s dolní mezí CI > 0. Nehodnotí se jen poloha maxima.

**Rozhodovací pravidlo (předem):** interval je „informativní“, pokud splní všechny čtyři podmínky:
1. dolní mez simultánního pásu > 0 aspoň v jedné škále (M-A nebo M-B) a směr se shoduje i ve druhé;
2. coverage ≥ 80 % při limitu 0,25 a nezávisí kvalitativně na limitu;
3. stejný směr v M-C a M-E;
4. absolutní Δ ≥ 1 min/den na bin. Hranice 1 min/den je pragmatická, jde o citlivost.

Pokud žádný souvislý rozsah podmínky nesplní, **závěr zůstává deskriptivní**. Přiměřenost jednotlivých hranic posoudí oponent.

### B4. Validační plán (fáze 5)

- **Datové testy:** D1–D12 jako `stopifnot` v kódu, okna jen z celých binů.
- **Modelové:** konvergence / singularita, rezidua vs fitted per bin, vlivná pozorování (Cookova vzdálenost na úrovni pacienta), citlivost na c v log škále.
- **Reprodukční:** stejné výsledky sekvenčně a paralelně (1 vs N workerů) v toleranci 1e-10, dvojí běh = shoda.
- **Vědecký sanity check:** vysvětlitelný rozdíl kumulativních a binových profilů (T_x = Σ binů ≥ x, to lze ověřit matematicky). Nezaměňovat skupinovou precision s individuální reliabilitou.

### B5. Literatura

Následující tvrzení **nejsou zatím ověřena** (doplní se v rámci fáze 2/4 přes PubMed/DOI):
- Rowlands 2018 (IG, MSSE): DOI uveden v předchozích podkladech.
- Rowlands 2019 (MX, Sports Med Open).
- PMID 32920327 (IIM, prospektivní kohorta ENMO).
- Literatura k zero-handling / kompozičním datům (PMID 33023619).
- Literatura ke Gamma/log modelům PA dat (PMID 26881287).

---

## C. Další vědecké otázky (návrh)

| # | Otázka | Hypotéza | Data | Model | Omezení | Priorita |
|---|---|---|---|---|---|---|
| C1 | Liší se **tvar** profilu odpovědi mezi M3, M6 a M18? | Časná odpověď (M3) je posunutá k vyšším intenzitám než pozdní, na relativní škále | bin_long M0–M18 | Společný LMM / Gamma GLMM `minutes ~ visit × bin + (1 \| patient) + (1 \| patient:bin)`, test interakce visit × bin (LRT / bootstrap) | Korelace binů v rámci záznamu; kompoziční závislost | **1** |
| C2 | Zkresluje chybějící sledování (drop-out) odhad změny? | Pacienti bez M18 se baseline profilem neliší (hrubě D20) | record_level, bin_long | Srovnání baseline profilu completers vs non-completers; LMM (MAR) vs complete-pair; pattern-mixture citlivost | Malé n (24 non-completers M18) | **1** |
| C3 | Je změna stabilní posun nebo vzácná vysokointenzitní aktivita u menšiny? | Ve středním pásmu roste většina pacientů, ve vysokém malá podskupina | Párová Δ per pacient | Distribuce individuálních Δ (podíl s Δ > 0, kvantily), koncentrace (podíl celkové Δ připadající na top-10 % pacientů) | Bez test-retest nelze oddělit šum od skutečné změny | 2 |
| C4 | Citlivost na počet validních dnů a sezónu | Závěry se při ≥ 5 dnech a s adjustací na sezónu nemění | record_level (`Nvaliddays`, `season`) | Refit primární analýzy na podmnožině, s opraveným anchor testem | Mírné confounding sezóna × visit (D19) | 2 |
| C5 | Konzistence s IG a MX metrikami | Změna IG/M-metrik odpovídá lokalizaci změny ve spektru | record_level (IG, M5, p-percentily) | Deskriptivní srovnání Δ IG a Δ percentilů s binovým profilem; žádný kauzální model | IG je funkcí stejných dat (není nezávislá validace) | 3 |

Budoucí otázky, mimo současný scope a bez souhlasu neimplementovat:
- Vztah změn ve spektru k MMT8 / FI2 / HAQ.
- Rozdíly mezi subtypy DM / IMNM / ASyS (vyžaduje legendu kódů, D16).

---

## D. Otevřená rozhodnutí pro uživatele (před fází 3–4)

1. **Běhové prostředí:**
   - (a) obnovit Linux renv knihovnu ve WSL (`renv::restore()`, desítky minut, možné systémové závislosti);
   - (b) analytické běhy spouštět ve Windows R 4.5.2, kde vznikly dosavadní výstupy.
2. **Fáze 3 — nezávislý oponent:**
   - Jiný model je k dispozici přes subagenta (např. Fable 5.1 nebo Sonnet 5.5, tedy jiná modelová identita než Opus 5.5).
   - Alternativně externí oponent na straně uživatele.
3. **D16:** legenda kódů subtypu 0–5 a zařazení pacienta „UCTD/SSc?“.
4. **Umístění výstupů:**
   - skripty `scripts/WRK_07_…` a dál;
   - objekty `output/RData/`;
   - figury a tabulky `output/thresholds/`;
   - report `reports/YYMMDD_thresholds_report.Rmd` (dle RULE §1 a §5).
