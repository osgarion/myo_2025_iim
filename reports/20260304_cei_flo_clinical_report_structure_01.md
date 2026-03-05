# Návrh Struktury Klinického Reportu

## 1. Shrnutí pro klinika (1 strana)
- `Cíl:` zda lze `ig_gradient` použít u pacientů, kteří dosahují stropu/dna klinických škál.
- `Hlavní závěr:` krátká odpověď pro `mmt8`, `mmt10`, `fi2`, `fi3`, `haq` (ano / částečně / ne).
- `Praktické doporučení:` kdy použít klasickou škálu, kdy doplnit `ig_gradient`, kdy je vhodný jako náhradní marker.

## 2. Klinický problém
- Proč jsou strop/dno efekt problém (ztráta citlivosti, nemožnost zachytit zlepšení/zhoršení).
- Co řeší `ig_gradient` a proč je kandidát na doplnění/náhradu.

## 3. Data a kohorta
- Popis dat: pacient, návštěva (`M0, M3, M6, M18, M30, M48`), klinické škály, `ig_gradient`.
- Definice stropu/dna pro každou škálu.
- Počty pacientů a měření:
  - celá kohorta,
  - kohorta bez stropu/dna,
  - kohorty se stropem.

## 4. Analytický plán (stručně, klinicky srozumitelně)
- EDA:
  - rozdělení hodnot v čase,
  - Kendall korelace (`ig_gradient` vs klinické škály),
  - speciální grafy pro pacienty se stropem.
- Převodní modely (`ig_gradient -> klinická škála`):
  - hurdle přístup (část strop/dno + část mimo strop/dno),
  - out-of-fold vyhodnocení.
- Disease activity modely (PhVAS):
  - modely A/B/C,
  - porovnání výkonu,
  - post-hoc rozdíly mezi návštěvami pro model 3.

## 5. Výsledky podle klinické otázky

### 5.1 Jak velký je problém stropu/dna?
- Podíl strop/dno podle škály a návštěvy.
- Kde škály přestávají být citlivé.

### 5.2 Má `ig_gradient` informaci mimo strop/dno?
- Korelace a směrovost vztahu.
- Stabilita napříč návštěvami.

### 5.3 Dokáže `ig_gradient` predikovat klinickou škálu?
- MAE / RMSE / kalibrace / BA-like chyba.
- Srovnání s baseline (průměr).

### 5.4 Co se děje u pacientů se stropem?
- Trendy `ig_gradient` v čase u ceiling/floor pacientů.
- Vztah k PhVAS (model 3 + post-hoc mezi návštěvami).

## 6. Rozhodovací tabulka pro praxi (nejdůležitější část)
- 1 řádek = 1 škála (`mmt8`, `mmt10`, `fi2`, `fi3`, `haq`).
- Sloupce:
  - `Míra stropu/dna`,
  - `Síla vazby na ig_gradient`,
  - `Výkon převodu`,
  - `Vazba na PhVAS`,
  - `Doporučení`.

Příklad doporučení:
- `Nahradit nelze, ale vhodný doplněk`.
- `Použitelné jen u podskupiny`.
- `Potenciální náhradní marker ve stropové části`.

## 7. Co by měl klinik dělat v ordinaci (implementace)
- Jak kombinovat škálu + `ig_gradient` při stropu.
- Jak interpretovat změnu `ig_gradient` mezi návštěvami.
- Jednoduchý algoritmus rozhodnutí (3 kroky).

## 8. Omezení a jistota závěrů
- Velikost vzorku a počet měření na pozdních návštěvách.
- Citlivost na složení kohorty.
- Co je robustní závěr a co je jen exploratorní.

## 9. Přílohy
- Přehled exportovaných tabulek (`output/ceil_flo/tables`).
- Přehled klíčových obrázků (`output/ceil_flo/figures`).
- Definice proměnných a zkratek.

## 10. Doporučená struktura finálního sdělení (na 5 minut)
- `1 minuta:` problém stropu/dna.
- `2 minuty:` co ukázal `ig_gradient` mimo strop a ve stropu.
- `1 minuta:` lze / nelze nahradit konkrétní škály.
- `1 minuta:` praktický postup pro klinickou práci.
