# MYO-II / Paper A — Kritické posouzení návrhů A07 (revize reportu)

- **Datum:** 2026-10-16
- **Podklad:** `data/threshold/MYO_A07_narhy_na_zmeny.txt`
- **Autor revize:** Claude Opus 5.5 (pro Veroniku Balounovou)
- **Rozhodnutí uživatelky (2026-10-16):**
  - Klinická data pro blok E: `data/threshold/klinicka_data_12102026.xlsx`.
  - Hlavním ukazatelem závažnosti je **ph_VAS**, citlivostní analýza **MMT8**.
  - Nezávislá oponentura (Claude Fable 5.1) proběhne v jednom kole nad výsledky.
- **Princip:** původní skripty 01–06, jejich výstupy a report `261013_thresholds_report.Rmd` zůstávají beze změny. Revize je ve skriptech 07–10, v helperu `R/06_helpers_revision.R` a v reportu `reports/261016_thresholds_report_v2.Rmd`.

Výsledky (čísla) jsou v reportu v2, v `docs/261016_A07_validation_change_log_summary.md` a v tabulkách `output/tables/261016_07…10_*.xlsx`. Tento dokument zachycuje **rozhodnutí a jejich zdůvodnění**.

---

## 1. Souhrnná tabulka

| # | Návrh | Rozhodnutí | Statistické zdůvodnění | Změna estimandu? |
|---|---|---|---|---|
| A | Normalizovat popisný profil binů na 10 mg | **Přijato** | `01_QC` kreslí `value` (surové min/den v binu). Biny mají šířku 10 / 20 / 25 mg, takže při zdvojnásobení šířky v 120 mg vzniká skok nahoru. U E1 je normalizace správně už od plánu v2 (MO2). Normalizuje se po záznamech (`minutes_per_10mg`, kontrola D11). Kvantily jsou vůči kladnému škálování ekvivariantní, proto je výsledek totožný s normalizací souhrnů. Nic se neinterpoluje. | Ne (jen popis) |
| B1 | Zachovat původní podíl top 10 % na čistém zisku | **Přijato, s označením limitu** | Vzorec: `sum(k největších Δ) / sum(Δ)`, kde k = ceiling(0,1·n). Když se zisky a ztráty téměř vyruší, jmenovatel jde k 0 a hodnoty přesahují 100 %. Při `sum(Δ) ≤ 0` je výsledek NA, a proto křivky začínají na různých intenzitách. Nově se počítá index `net/gross = ΣΔ / Σ|Δ|`. Buňky s indexem < 0,20 jsou označené jako nestabilní. | Ne |
| B2 | Nový estimand: podíl top 10 % na kladných ziscích | **Přijato, upraveno** | Omezení na 0–100 % a přímá interpretace jako koncentrace zisků. Doplněno o **referenční hodnotu** při rovnoměrném rozdělení, `min(k, n_pos) / n_pos`, protože bez ní se koncentrace posoudit nedá. k se počítá ze **všech** párů buňky, takže top 10 % pacientů je stejná skupina jako v původním estimandu. Přidány 95% intervaly z pacientského cluster bootstrapu (B = 2000, seed 20261013). | **Ano, nový estimand** (doplněk) |
| B3 | Výběr top 10 %, shody, neceločíselné k | **Upřesněno** | Výběr probíhá zvlášť v každém binu a návštěvě. k = ceiling(0,1·n), skutečný podíl je k/n (např. 5/47 = 10,6 %). Shody na k-té hodnotě součet nemění. Pokud je kladných změn méně než k, sčítají se všechny kladné a podíl je 100 %. | — |
| B4 | Nízké a vysoké intenzity zvlášť | **Přijato** | Panely 20–120 mg (biny 10 mg) a 120–400 mg (biny 20/25 mg) mají vlastní osy y. Původní estimand je ve zvláštní figuře s log osou a s vyznačenými nestabilními buňkami. | — |
| B5 | Rozlišit % se ziskem, koncentraci a podíl na čisté změně | **Přijato** | Tři metriky jsou ve třech oddělených panelech a sloupcích tabulky. Text reportu (H3 „menšina silných responderů“) se přepisuje podle nových čísel. Navíc je uvedena koncentrace pro souhrnné T100 a T200 (čas nad prahem na pacienta). Ta je robustnější než jednotlivé vysoké biny s malým objemem. | — |
| C1 | Ověřit CLR, náhradu nul, 21dílnou kompozici | **Provedeno** | Uzavření (Σclr = 0, součet 1440 min) a počty nahrazených nul podle části a návštěvy. | — |
| C2 | Vliv nestejné šířky binů na trend a lokální odchylky | **Posouzeno analyticky i numericky** | Párová změna log(minuty·10/šířka) se rovná párové změně log(minuty), protože konstanta šířky se odečte. Centrum CLR je v rámci návštěvy konstantní a multiplikativní přeškálování při náhradě nul je konstanta na záznam. Obojí pohltí intercept splinu, takže **rezidua E3 nezávisí na normalizaci ani na centrování** (ověřeno na 1e-16). Šířka se projeví jen dvěma cestami: (i) **polohou a hustotou bodů** na ose log-intenzity, kde je 10 bodů v 20–120 mg a 9 v 120–300 mg, a (ii) **rozlišením**, protože širší bin průměruje přes větší rozsah intenzit. | — |
| C3 | Alternativní seskupení / šířkově vážená reprezentace | **Přijato, se třídou varianty** | *Citlivost stejného estimandu:* spline vážený log-šířkou binu (každý bin váží úměrně rozsahu log-intenzity, který pokrývá). *Změna rozlišení:* jednotná mřížka 20 mg (sloučení dvojic nativních 10mg binů, 14 binů). *Změna rozsahu kompozice:* ohraničené biny až do 400 mg (23 binů + ≥ 400). Poslední dvě mění estimand, protože pravidlo „≥ 2 sousední biny“ znamená jinou minimální šířku hrbu. Reportují se proto jako robustnost závěru, ne jako sensitivita stejného čísla. | Částečně (označeno) |
| C4 | Odlišit sensitivitu od arbitrárního přerozdělení | **Přijato** | Každá varianta má sloupec `class`. Arbitrární přerozdělení (např. kvantilové biny) se **nedělá**, protože vyžaduje interpolaci a mění populaci binů podle dat. | — |
| C5 | Spliny odděleně pro M3/M6/M18 se stejnou specifikací | **Ověřeno** | `est_E3()` fituje `lm(y ~ ns(x, df = 3))` pro každou návštěvu zvlášť. Báze je u všech návštěv identická, protože x se nemění (kontrola `same_basis`). | — |
| C6 | Volba df = 3 (simulace / citlivost) | **Přijato** (df 3 ponecháno jako předem zvolené, ne jako prokazatelně optimální; po oponentuře F3) | (i) Pozorovaný verdikt pro df 1–5. (ii) Plasmode „true null“ pro každé df testu 2–5 se správnou nulou (df 3) i s **misspecifikovanou nulou** (skutečný hladký průběh má df 5). Měří se míra falešných poplachů a síla pro hrby +0,1 / +0,2 log na 2 a 6 binech. | — |
| C7 | Jak se počítají pásy a odchylky | **Popsáno + ověřeno** | Max-t pás nad 19 rezidui návštěvy z pacientského cluster bootstrapu. Pásy jsou nezávisle přepočítány ze 2000 uložených replikací (base R). | — |
| C8 | Panely A/B v reportu | **Přepsáno** | Panel A: body, CI a viditelná hladká křivka návštěvy. Šedé pole vyznačuje „bin zvětšuje podíl dne vůči geometrickému průměru všech 21 částí“. Panel B: rezidua, simultánní pás a explicitní kritéria lokální odchylky. V textu je výslovně uvedeno, že nepřítomnost lokální odchylky **neznamená**, že jsou všechny oblasti klinicky stejně informativní. | — |
| D1 | Jsou rozdíly hladkých trendů M3/M6/M18 testované? | **Ne, jen C1 (M6 − M3 po binech)** | Doplněno. | — |
| D2 | Párové srovnání M3, M6, M18 | **Přijato** | U pacientů s M0 a oběma návštěvami se počítá Δ = clr změna (pozdější) − clr změna (dřívější). Hodnotí se: (i) biny s max-t pásem, (ii) rozdíl lineárního sklonu, (iii) globální test, že rozdíl hladkých křivek je plochý (bootstrap Wald na 3 nekonstantních koeficientech ns(df 3); konstantní posun je identifikovaný, odráží výměnu s částmi <20 a ≥300 mg, a záměrně se netestuje, protože otázkou je tvar). P-hodnoty jsou Holm-korigované přes 3 páry. | **Ano, nové estimandy** |
| D3 | Různé vzorky u návštěv | **Přijato** | Hlavní analýza: všichni pacienti s M0 a oběma návštěvami. Citlivost: společný vzorek pacientů se všemi 4 návštěvami. | — |
| D4 | Sezóna, překryv oken | **Popsáno, neadjustováno** | Tabulka měsíců od M0 a podílu stejné sezóny v páru. Sezóna je s návštěvou strukturálně svázaná (žádný pár M0→M6 není ve stejné sezóně), proto ji nelze od návštěvy oddělit modelem. Adjustace by byla jen extrapolací. | — |
| D5 | Vizuální rozdíl ≠ rozdílná trajektorie | **Přijato** | Panel A figury D1 je označený jako deskriptivní. Závěr se opírá jen o testy z D2. | — |
| E1 | Dostupnost ph_VAS, MMT8, FI2 | **Ověřeno** | Klinický soubor je mimo master. Propojeno 71/71 pacientů, kódy podtypů souhlasí. Dostupnost v M0 viz tabulka `availability`. FI2 má nejvíc chybějících hodnot. | — |
| E2 | Reprezentuje ph_VAS zamýšlený konstrukt? | **Částečně** | Physician global activity (MDAAT) zahrnuje i extramuskulární aktivitu (kůže, plíce, klouby). Ta se mezi podtypy systematicky liší: IMNM má téměř nulovou extramuskulární aktivitu, DM a ASyS vysokou. Stejná hodnota ph_VAS proto u DM a ASyS znamená menší svalové postižení než u IMNM. To snižuje schopnost ph_VAS adjustovat „počáteční svalové postižení“ napříč podtypy. Hlavní je ph_VAS podle rozhodnutí uživatelky, MMT8 je citlivost. MDAAT muscle VAS a extramuskulární VAS se reportují deskriptivně. | — |
| E3 | Kolinearita | **Posouzeno** | Spearman matice, VIF závažnosti vůči podtypu (R² severity ~ subtype). FI2 se do modelu s MMT8 nedává (ρ ≈ 0,85). | — |
| E4 | Úsporné modely | **Přijato** | Jen tři předem dané souhrny změny (ΔT100, sklon clr, ΔIG), žádný výběr binu. Random-intercept LMM změny: `visit * subtype` (+ `sev`) (+ `sev:visit`) (+ výchozí hodnota výsledku). Nejvíc 12 fixních parametrů na 47–48 pacientů a 126–127 pozorování. Kontrasty se počítají při průměrné závažnosti a mají 95% pacientský cluster bootstrap. LRT (ML) pro subtype × visit, globální efekt podtypu a severity × visit. | **Ano, nové (exploratorní)** |
| E5 | subtype × visit, severity × visit | **Přijato** | Jako LRT a varianta E_c. | — |
| E6 | Zachovat F7 | **Ano** | F7 je beze změny. Nová F7b ukazuje kontrasty. Řádky 1 a 2 oddělují vliv **omezení vzorku** (chybějící ph_VAS, hlavně u ASyS) od vlivu **adjustace**. | — |
| E7 | Exploratorní povaha | **Ano** | Žádná korekce na multiplicitu, výsledky jako hypotézy. Bez dalších kovariát (trvání nemoci, GK) kvůli overfittingu. Ty jsou jen v tabulce rozložení. | — |

---

## 2. Co jsem nepřijala a proč

1. **Interpolace syntetických 10mg binů** (výslovně zakázáno v A07). Nepoužito ani pro E3. Slučují se jen nativní hrany GGIR.
2. **Kvantilové / datově řízené přerozdělení binů pro E3.** Mění estimand podle dat a vyžaduje interpolaci, proto nebylo provedeno.
3. **Adjustace na sezónu v D.** Při strukturální vazbě sezóny na návštěvu není identifikovatelná. Je jen popsána.
4. **Rozšíření modelů podtypů na celé spektrum (23 binů × 3 návštěvy × 3 podtypy).** Při zhruba 15 pacientech na podtyp by šlo o mnohonásobné testování s nulovou silou. Zůstávají předem dané tři souhrny.
5. **Přepočet primárního bootstrapu (02).** Primární estimandy E1 až E4 se nemění. Skript 10 ověřuje, že pozorované hodnoty se reprodukují.

## 3. Rizika a otevřené body

- **Klinický soubor** (`klinicka_data_12102026.xlsx`) není součástí uzamčeného masteru. Jeho md5 je v logu skriptu 09. Při změně souboru se blok E musí přepočítat.
- **ph_VAS chybí** u části pacientů s párem, nejvíc u ASyS. Adjustované modely proto stojí na menším a vybraném vzorku (řádek 2 F7b to ukazuje).
- **Plasmode** testuje rozlišovací schopnost pravidla jen pro zvolené tvary hrbů (2 a 6 binů, +0,1 / +0,2 log). Jiné tvary nejsou pokryté.
- **Denní data GGIR** (`d11_part2daysummary` v `data/260123_backup.RData`) by mohla otevřít modul reliability (ICC/SEM/SDC), který byl v plánu v2 odložen. Mimo rozsah A07, doporučeno jako další krok.
