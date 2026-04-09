dir.create(file.path("output"), showWarnings = FALSE, recursive = TRUE)

normalize_rel <- function(path) {
  gsub("\\\\", "/", path)
}

dir_paths <- list.dirs("output", recursive = TRUE, full.names = TRUE)
dir_paths <- c("output", setdiff(dir_paths, "output"))
dir_paths <- unique(normalize_rel(dir_paths))

direct_files <- function(dir_rel) {
  full <- file.path(getwd(), dir_rel)
  files <- list.files(full, full.names = FALSE, recursive = FALSE)
  files <- files[tolower(files) != "readme.md"]
  info <- file.info(file.path(full, files))
  files[!is.na(info$isdir) & !info$isdir]
}

direct_subdirs <- function(dir_rel) {
  full <- file.path(getwd(), dir_rel)
  dirs <- list.dirs(full, recursive = FALSE, full.names = FALSE)
  setdiff(dirs, basename(full))
}

format_paths <- function(x) {
  x <- unique(x)
  if (length(x) == 0) {
    return("- žádné specifické cesty nejsou pro tuto složku evidovány")
  }
  paste0("- `", x, "`", collapse = "\n")
}

format_files <- function(x, limit = 12) {
  if (length(x) == 0) {
    return("- v této složce nejsou přímo uložené žádné soubory")
  }
  shown <- head(sort(x), limit)
  lines <- paste0("- `", shown, "`")
  if (length(x) > limit) {
    lines <- c(lines, paste0("- ... a dalších ", length(x) - limit, " souborů stejného typu"))
  }
  paste(lines, collapse = "\n")
}

format_subdirs <- function(x) {
  if (length(x) == 0) {
    return("- tato složka nemá další podsložky")
  }
  paste0("- `", sort(x), "`", collapse = "\n")
}

count_extensions <- function(files) {
  if (length(files) == 0) {
    return("žádné přímé soubory")
  }
  ext <- tools::file_ext(files)
  ext[ext == ""] <- "[bez_přípony]"
  tab <- sort(table(ext), decreasing = TRUE)
  paste(paste0(names(tab), ": ", as.integer(tab)), collapse = ", ")
}

common_core_inputs <- c(
  "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Data/ (soubory vyhledávané podle názvů obsahujících `Myokiny`, `IMMET_terapie`, `Treatment Response`)",
  "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Data/MSTN normalizácia na kreat, BCM.xlsx",
  "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/GENEActiv/GGIR/merged_table_final_pro_zaverecnou_zpravu_AZV_100226.xlsx",
  "data/processed/Myokiny sérum_pro stat_processed.xlsx",
  "data/processed/selected_colnames_01.xlsx"
)

common_core_objects <- c(
  "Objekty z `scripts/functions/OBJ_01.R`: `d01_myokiny`, `d01_myokiny_2`, `d02_clinics`, `d05_norm`, `d10_ostatni`, `d10_cas`, navazující spojované objekty jako `d11_clinics_imp`, `d12_join`, `d10_cas_b`, `d10_aktivita`",
  "Pracovní workspace a tabulkové objekty ukládané do `reports/markD_01.RData`",
  "Záloha obrazu prostředí v `data/<datum>_backup.RData` nebo `data/260123_backup.RData`"
)

ggir_report_inputs <- c(
  "Objekt `d18_ggir_data` dostupný v prostředí po načtení `OBJ_01.R` nebo z reportového workspace",
  "Reportový balík `reports/260124_ggir_clinics_01.RData`",
  "R Markdown report `reports/260124_ggir_clinics_01.Rmd` a jeho výstup `reports/260124_ggir_clinics_01.html`"
)

tis_inputs <- c(
  "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/GENEActiv/GGIR/merged_table_final_klinická.xlsx",
  "Objekt `d18_tis_comp_01` ze skriptu `scripts/Script_myo_2025_iim_tis_01.R`",
  "Exportovaný TIS scale dataset `data/processed/260310_d18_tis_scale_01.xlsx`",
  "Soubor `data/processed/260310_d18_tis_imp_01.xlsx` v repozitáři existuje, ale v aktuální verzi skriptu jsem nenašel jeho explicitní export; je pravděpodobné, že jde o starší nebo ručně uložený mezikrok"
)

ceil_inputs <- function() {
  c(
    "R:/MYOZITIDY_VÝZKUM/_IMMET-GRANT/IMMET štatistika/Data/Raw/Geneactiv/20260302_GGIR_klinika/260302_data_geneactive_clinics.xlsx",
    "Načítání zajišťuje `scripts/functions/OBJ_ceil_flo_01.R`, které vytváří objekt `d18_ggir_data`",
    "Konfigurační override je možný přes `options(myo_ceil_flo_data_path_01 = \"...\")`"
  )
}

cluster_inputs <- c(
  "Objekty z `OBJ_01.R` a navazujících částí `Script_myo_2025_iim_01.R`, zejména `d13_join`, `d13_mofa` a odvozené matice pro FMDA/UMAP/PLS",
  "Skripty `scripts/cluster_plsda_pipeline.R` a `scripts/cluster_plsda_pipeline_2.R` zapisují další výstupy do stejného prostoru `output/figures/cluster analyses`"
)

describe_dir <- function(dir_rel) {
  dir_rel <- normalize_rel(dir_rel)

  meta <- list(
    title = paste("Dokumentace složky", dir_rel),
    analysis = "Technická dokumentace výstupní složky.",
    audience = "Složka slouží jako součást analytického archivu projektu a obsahuje buď finální výstupy, nebo tematicky sdružené mezivýstupy.",
    methods = c("Metoda zde není vykonávána přímo; tato složka je kontejner pro výstupy vytvořené jinými skripty."),
    scripts = character(),
    inputs = character(),
    objects = character(),
    naming = "Názvy souborů v projektu typicky začínají datem ve formátu `YYMMDD`, po kterém následuje stručný popis analýzy, hlavní outcome/predictor a pořadové číslo verze.",
    ai_steps = c(
      "Nejprve zkontroluj obsah této složky a přečti README v nadřazené i podřízených složkách.",
      "Teprve potom otevírej příslušné zdrojové skripty a ověřuj, které objekty se do výstupů propisují.",
      "Pokud chceš výstupy obnovit, preferuj znovuspuštění zdrojového skriptu místo ruční editace tabulek nebo obrázků."
    ),
    caveats = character()
  )

  if (dir_rel == "output") {
    meta$title <- "Dokumentace kořenové složky output"
    meta$analysis <- "Kořenová složka `output` sdružuje veškeré analytické výstupy projektu: starší průřezové a longitudinální tabulky, korelační a asociační obrázky, cluster analýzy, radarové porovnání modelů, TIS workflow i novější ceiling/floor analýzy."
    meta$audience <- "Pro nového čtenáře je to hlavní rozcestník k finálním i pracovním výstupům. Pro AI agenta je to místo, odkud se vyplatí začít mapování toho, které skripty generují které složky."
    meta$methods <- c(
      "Samotná složka nic nepočítá; agreguje výsledky z několika nezávislých analytických větví.",
      "Historicky starší výstupy jsou soustředěny hlavně v `output/figures` a `output/tables`.",
      "Tematicky čistší novější workflow používají vlastní podsložky jako `ggir`, `sf36`, `tis`, `*_ggir_eff`, `ceil_flo*` a `exam_incl`."
    )
    meta$scripts <- c(
      "scripts/Script_myo_2025_iim_01.R",
      "scripts/Script_myo_2025_iim_new_cas_01.R",
      "scripts/Script_myo_2025_iim_tis_01.R",
      "scripts/Script_myo_2025_iim_ceilingEffect_01.R a wrappery pro jednotlivé GGIR prediktory",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_01.R až scripts/Script_myo_2025_iim_clinics_ggir_vis_08.R",
      "scripts/cluster_plsda_pipeline.R a scripts/cluster_plsda_pipeline_2.R"
    )
    meta$inputs <- unique(c(common_core_inputs, ggir_report_inputs, tis_inputs, ceil_inputs()))
    meta$objects <- unique(c(common_core_objects, "Objekty `d18_tis_comp_01` a `d18_ggir_data` používané v TIS a GGIR radarových workflow"))
    meta$caveats <- c(
      "Kořen `output` obsahuje směs historických a novějších výstupů; jednotný naming existuje, ale analytická logika se mezi větvemi liší.",
      "Pokud chceš něco reprodukovat přesně, nepostačí sledovat jen datum v názvu souboru; musíš dohledat odpovídající skript a jeho vstupní verze dat."
    )
    return(meta)
  }

  if (dir_rel == "output/figures") {
    meta$analysis <- "Smíšený archiv obrázků z dřívějších analytických větví. Najdeš zde starší PDF grafy z `Script_myo_2025_iim_01.R`, obrázky generované z `Script_myo_2025_iim_new_cas_01.R` a tematické podsložky `cluster analyses`, `final` a `tis`."
    meta$audience <- "Pro člověka jde hlavně o pracovní grafický archiv. Pro AI agenta je důležité vědět, že tahle složka není jedna metodická větev, ale několik rozdílných generátorů výstupů."
    meta$methods <- c(
      "PDF soubory pocházejí hlavně ze starších explorativních a modelových výstupů v `Script_myo_2025_iim_01.R` a `Script_myo_2025_iim_new_cas_01.R`.",
      "PNG/TIFF soubory reprezentují korelační heatmapy, asociační grafy, kanonické mapy, myokine axes a další finální figury."
    )
    meta$scripts <- c(
      "scripts/Script_myo_2025_iim_01.R",
      "scripts/Script_myo_2025_iim_new_cas_01.R",
      "scripts/Script_myo_2025_iim_tis_01.R"
    )
    meta$inputs <- common_core_inputs
    meta$objects <- common_core_objects
    meta$caveats <- c("Část souborů zde existuje paralelně i v podadresářích `final` a `tis`, které představují kurátorovanější podmnožiny.")
    return(meta)
  }

  if (dir_rel == "output/tables") {
    meta$analysis <- "Smíšený archiv tabulek z několika hlavních analytických skriptů. Tato složka obsahuje jak rané přehledové exporty z lineárních smíšených modelů, tak novější korelační, asociační, rCCA a radarové tabulky."
    meta$audience <- "Pro člověka je to hlavní tabulkový archiv projektu. Pro AI agenta je to užitečný index toho, které výsledky byly exportovány mimo reporty."
    meta$methods <- c(
      "Starší tabulky vytváří hlavně `Script_myo_2025_iim_01.R` pomocí exportu modelových koeficientů a kontrastů.",
      "Novější tabulky vytváří `Script_myo_2025_iim_new_cas_01.R` a `Script_myo_2025_iim_clinics_ggir_vis_01.R`, zejména korelace Kendall tau, asociační modely, rCCA a shrnutí radarového porovnání modelů."
    )
    meta$scripts <- c(
      "scripts/Script_myo_2025_iim_01.R",
      "scripts/Script_myo_2025_iim_new_cas_01.R",
      "scripts/Script_myo_2025_iim_tis_01.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_01.R"
    )
    meta$inputs <- common_core_inputs
    meta$objects <- common_core_objects
    meta$caveats <- c(
      "Podsložka `final` obsahuje ručně vybrané finální tabulky; ne všechny tabulky z kořene `output/tables` jsou určeny k publikaci.",
      "Podsložka `tis` je specifická pro TIS workflow a má vlastní zdrojový skript."
    )
    return(meta)
  }

  if (dir_rel == "output/figures/cluster analyses") {
    meta$analysis <- "Výstupy klastrovacích a latentně-prostorových analýz: FMDA, UMAP a PLS(-DA) vizualizace pro baseline i změnové profily."
    meta$audience <- "Novému čtenáři tato složka ukazuje, jak byly hledány vzory v mnohorozměrných datech. AI agent zde najde prostor pro navazující latentní nebo klastrovací reprodukce."
    meta$methods <- c(
      "FMDA grafy vznikají z vícerozměrného mapování klinických a biomarkerových profilů v `Script_myo_2025_iim_01.R`.",
      "UMAP obrázky reprezentují nelineární redukci rozměru pro průzkum podobnosti pacientů.",
      "PLS/PLS-DA pipeline navazuje ve skriptech `cluster_plsda_pipeline.R` a `cluster_plsda_pipeline_2.R`."
    )
    meta$scripts <- c(
      "scripts/Script_myo_2025_iim_01.R",
      "scripts/cluster_plsda_pipeline.R",
      "scripts/cluster_plsda_pipeline_2.R",
      "scripts/functions/FUN_01.R"
    )
    meta$inputs <- common_core_inputs
    meta$objects <- cluster_inputs
    meta$caveats <- c("Výstupy jsou výrazně závislé na tom, které proměnné byly v konkrétní verzi skriptu zařazeny a jak byly filtrovány případy s chybějícími hodnotami.")
    return(meta)
  }

  if (dir_rel == "output/figures/final") {
    meta$analysis <- "Kurátorovaná sada finálních TIFF obrázků z `Script_myo_2025_iim_new_cas_01.R`. Jde zejména o korelační heatmapy, kanonické vztahy, myokine axes, subtype/response ploty a další publikačně či reportově významné figury."
    meta$audience <- "Tato složka je nejblíže finální prezentaci výsledků. Člověk zde typicky hledá grafy použitelné do zprávy nebo rukopisu."
    meta$methods <- c(
      "Korelace: Kendall tau heatmapy pro aktivitu a čas.",
      "Multiblokové / kanonické vztahy: rCCA a příbuzné přístupy nad myokiny, aktivitou a časovými parametry.",
      "Výběr finálních obrázků je tvrdě zakódovaný přímo ve `Script_myo_2025_iim_new_cas_01.R`."
    )
    meta$scripts <- c("scripts/Script_myo_2025_iim_new_cas_01.R")
    meta$inputs <- common_core_inputs
    meta$objects <- c(common_core_objects, "Objekty typu `corr_res_*`, `rcc_fit_*`, `cor_df_*`, `iim_myo_rcca_*` vytvářené přímo ve `Script_myo_2025_iim_new_cas_01.R`")
    return(meta)
  }

  if (dir_rel == "output/figures/tis" || dir_rel == "output/tables/tis") {
    meta$analysis <- "Starší, širší TIS workflow z `Script_myo_2025_iim_tis_01.R`. Nejde o tři radarové exporty z `output/tis`, ale o rozsáhlejší sadu explorativních a validačních výstupů kolem výpočtu TIS."
    meta$audience <- "Pro člověka je to technická dokumentace výpočtu TIS a jeho vztahu k aktivitě, klinice a biomarkerům. Pro AI agenta jde o hlavní zdroj pro pochopení, jak je TIS v projektu definován."
    meta$methods <- c(
      "KNN imputace chybějících klinických hodnot pomocí `recipes::step_impute_knn`.",
      "Převod klinických změn vůči M0 na TIS škálu podle pravidel ve skriptu.",
      "Korelační a regresní analýzy TIS vs GGIR / další klinické parametry."
    )
    meta$scripts <- c("scripts/Script_myo_2025_iim_tis_01.R")
    meta$inputs <- tis_inputs
    meta$objects <- c(
      "Objekty `d18_ggir_tis_data`, `d18_tis_knn_imp`, `d18_tis_comp_01`, `d18_ggir_tis_data_02`",
      "TIS scale export `data/processed/260310_d18_tis_scale_01.xlsx`"
    )
    meta$caveats <- c("Protože soubor `260310_d18_tis_imp_01.xlsx` není v aktuálním skriptu explicitně exportován, ber jeho původ jako informovanou inferenci.")
    return(meta)
  }

  if (startsWith(dir_rel, "output/ceil_flo")) {
    path_parts <- strsplit(dir_rel, "/", fixed = TRUE)[[1]]
    branch_name <- if (length(path_parts) >= 2) path_parts[2] else "ceil_flo"
    predictor_label <- if (branch_name == "ceil_flo") {
      "výchozí variantu s prediktorem `ig_gradient_enmo`"
    } else {
      paste0("alternativní variantu s hlavním GGIR prediktorem `", sub("^ceil_flo_", "", branch_name), "`")
    }

    meta$analysis <- paste0(
      "Ceiling/floor workflow pro ", predictor_label, 
      ". Tato větev testuje, jak se mění interpretace klinických outcome proměnných, pokud část pacientů leží na stropu/podlaze škály a pokud je do modelů zařazen konkrétní GGIR prediktor."
    )
    meta$audience <- "Pro klinika jde o interpretaci toho, zda hraniční skóry deformují vztah mezi aktivitou z akcelerometru a klinickými škálami. Pro AI agenta je to izolované workflow s jasně definovaným vstupem a exporty."
    meta$methods <- c(
      "Identifikace ceiling/floor případů pro `mmt8`, `mmt10`, `fi2`, `fi3` a `haq`.",
      "EDA: hustoty, boundary/non-boundary grafy, Kendallovy asociace, grafy podle exam.",
      "Modelování disease activity a prediction performance na celé kohortě, na boundary subsets a na datech bez ceiling/floor případů.",
      "Wrapper skripty pouze nastavují hlavní prediktor a výstupní složku; samotnou logiku drží `Script_myo_2025_iim_ceilingEffect_01.R`."
    )
    wrapper <- switch(
      branch_name,
      ceil_flo = "scripts/Script_myo_2025_iim_ceilingEffect_01.R",
      ceil_flo_acc_day_spt_wei = "scripts/Script_myo_2025_iim_ceilingEffect_acc_day_spt_wei_01.R",
      ceil_flo_ig_intercept_enmo = "scripts/Script_myo_2025_iim_ceilingEffect_ig_intercept_enmo_01.R",
      ceil_flo_m5_enmo = "scripts/Script_myo_2025_iim_ceilingEffect_m5_enmo_01.R",
      ceil_flo_mvpa_t100_time_min = "scripts/Script_myo_2025_iim_ceilingEffect_mvpa_t100_time_min_01.R",
      ceil_flo_p9931_enmo = "scripts/Script_myo_2025_iim_ceilingEffect_p9931_enmo_01.R",
      "scripts/Script_myo_2025_iim_ceilingEffect_01.R"
    )
    meta$scripts <- c(
      wrapper,
      "scripts/Script_myo_2025_iim_ceilingEffect_01.R",
      "scripts/functions/FUN_ceil_flo_01.R",
      "scripts/functions/OBJ_ceil_flo_01.R"
    )
    meta$inputs <- ceil_inputs()
    meta$objects <- c(
      "Objekt `d18_ggir_data` vytvořený v `OBJ_ceil_flo_01.R` načtením legendy a sheetu `sheet1` z externího XLSX",
      "Odvozený long objekt `d18_ggir_data_long` s indikátorem `cei_flo`",
      "Modelové výsledky a EDA tabulky vytvářené přímo ve `Script_myo_2025_iim_ceilingEffect_01.R`"
    )
    meta$caveats <- c("Výstupní složky `ceil_flo_*` jsou paralelní alternativy nad stejnou datovou základnou; liší se hlavním GGIR prediktorem, nikoli základní klinickou logikou.")
    return(meta)
  }

  if (startsWith(dir_rel, "output/exam_incl")) {
    meta$analysis <- "Varianta radarových a modelových workflow, kde je návštěva (`exam`) explicitně ponechána v modelování. Větve bez `fi3` používají jako hlavní klinický adjustor FI-2, větve pod `exam_incl/fi3` používají FI-3."
    meta$audience <- "Tato složka je důležitá pro porovnání, co se stane po explicitním zohlednění longitudinální návštěvy. Pro klinika jde o citlivostní analýzu. Pro AI agenta o alternativní model specification."
    meta$methods <- c(
      "Smíšené lineární modely s covariate `exam` nebo analogickou faktorovou reprezentací návštěvy.",
      "Stejná radarová logika jako v hlavních složkách `ggir`, `ph_vas`, `sf36`, `tis`, ale s jinou specifikací modelu.",
      "Soubory končící `_exam_` patří větvi s FI-2; soubory `_exam_b_` patří variantě s FI-3."
    )
    meta$scripts <- c(
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_01_exam.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_03_exam.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_04_exam.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_05_exam.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_01_exam_b.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_03_exam_b.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_04_exam_b.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_05_exam_b.R"
    )
    meta$inputs <- c(ggir_report_inputs, tis_inputs)
    meta$objects <- c("Primárně `d18_ggir_data`; pro TIS větve navíc `d18_tis_comp_01`")
    return(meta)
  }

  if (dir_rel == "output/ggir") {
    meta$analysis <- "Modely, kde GGIR parametry vystupují jako outcome a klinické proměnné jako testované prediktory."
    meta$methods <- c(
      "Smíšené lineární modely ve tvaru `ggir_parameter ~ klinický_prediktor + age_m0 + sex + (1 | immet_id)`.",
      "Výsledky jsou agregovány do radarového skórování kombinující velikost efektu, CI, R2 a AIC."
    )
    meta$scripts <- c("scripts/Script_myo_2025_iim_clinics_ggir_vis_04.R")
    meta$inputs <- ggir_report_inputs
    meta$objects <- c("Objekt `d18_ggir_data` s klinickými a GGIR proměnnými")
    return(meta)
  }

  if (dir_rel == "output/sf36") {
    meta$analysis <- "Radarové porovnání modelů pro SF-36 outcome proměnné. Pro každou SF-36 doménu se porovnává model s FI-2 adjustací, model s MMT8 adjustací a single-predictor modely."
    meta$methods <- c(
      "Smíšené lineární modely nad jednotlivými doménami `sf36_*`.",
      "Skórování modelů využívá standardizovaný beta, přesnost intervalu spolehlivosti, semi-partial R2, marginal/conditional R2 a obrácený AIC."
    )
    meta$scripts <- c("scripts/Script_myo_2025_iim_clinics_ggir_vis_03.R")
    meta$inputs <- ggir_report_inputs
    meta$objects <- c("Objekt `d18_ggir_data` připravený funkcí `prepare_ggir_sf36_data_vis03()`")
    return(meta)
  }

  if (dir_rel == "output/tis") {
    meta$analysis <- "Radarové porovnání modelů pro outcome `tis_total` po sloučení TIS výsledků s GGIR daty."
    meta$methods <- c(
      "Vstupní GGIR data se spojí s TIS daty přes `immet_id` a `exam`.",
      "Porovnávají se FI-2 adjusted, MMT8 adjusted a single-predictor modely."
    )
    meta$scripts <- c("scripts/Script_myo_2025_iim_clinics_ggir_vis_05.R")
    meta$inputs <- c(ggir_report_inputs, tis_inputs)
    meta$objects <- c("Objekty `d18_ggir_data` a `d18_tis_comp_01`")
    return(meta)
  }

  if (dir_rel == "output/tis_ggir_eff") {
    meta$analysis <- "Srovnání přínosu TIS a GGIR pro vysvětlení SF-36 PCS. Pro každý GGIR parametr se staví trojice modelů: `tis_only`, `ggir_only`, `tis_ggir`."
    meta$methods <- c(
      "Outcome je `sf36_pcs`, v log variantě `log(sf36_pcs)`.",
      "Covariates: `exam + age_m0 + sex + sf36_mcs + (1 | immet_id)`.",
      "Ke každému GGIR parametru se exportuje tabulka, hlavní graf a diagnostický graf."
    )
    meta$scripts <- c("scripts/Script_myo_2025_iim_clinics_ggir_vis_06.R")
    meta$inputs <- c(ggir_report_inputs, tis_inputs)
    meta$objects <- c("Spojený objekt z `merge_tis_into_ggir_vis06()` a model bundles pro raw/log outcome")
    return(meta)
  }

  if (dir_rel == "output/fi2_ggir_eff") {
    meta$analysis <- "Srovnání přínosu FI-2 a GGIR pro vysvětlení SF-36 PCS. Pro každý GGIR parametr se staví modely `fi2_only`, `ggir_only` a `fi2_ggir`."
    meta$methods <- c(
      "Outcome je `sf36_pcs`, plus paralelní log verze `log(sf36_pcs)`.",
      "Covariates: `exam + age_m0 + sex + sf36_mcs + (1 | immet_id)`.",
      "Export obsahuje tabulku, hlavní porovnávací graf a diagnostiku."
    )
    meta$scripts <- c("scripts/Script_myo_2025_iim_clinics_ggir_vis_07.R")
    meta$inputs <- ggir_report_inputs
    meta$objects <- c("Objekt připravený funkcí `prepare_fi2_ggir_eff_data_vis07()`")
    return(meta)
  }

  if (dir_rel == "output/fi3_ggir_eff") {
    meta$analysis <- "Srovnání přínosu FI-3 a GGIR pro vysvětlení SF-36 PCS. Struktura je analogická k `fi2_ggir_eff`, ale základní klinický prediktor je FI-3."
    meta$methods <- c(
      "Outcome je `sf36_pcs`, plus log verze.",
      "Covariates: `exam + age_m0 + sex + sf36_mcs + (1 | immet_id)`.",
      "Výstupní grafy a tabulky jsou generovány po jednom pro každý GGIR parametr."
    )
    meta$scripts <- c("scripts/Script_myo_2025_iim_clinics_ggir_vis_08.R")
    meta$inputs <- ggir_report_inputs
    meta$objects <- c("Objekt připravený funkcí `prepare_fi3_ggir_eff_data_vis08()`")
    return(meta)
  }

  if (grepl("/ph_vas$", dir_rel)) {
    meta$analysis <- "Radarové porovnání modelů pro outcome physician VAS / PhVAS s GGIR prediktory a klinickými adjustory."
    meta$methods <- c(
      "Srovnávají se modely adjustované o FI-2 nebo FI-3/MMT8 podle konkrétní větve.",
      "Výstupem je radarový graf a XLSX s dílčími metrikami modelů."
    )
    meta$scripts <- c(
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_01.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_01_exam.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_01_exam_b.R"
    )
    meta$inputs <- ggir_report_inputs
    meta$objects <- c("Objekt `d18_ggir_data`; příprava probíhá funkcemi `prepare_ggir_model_data()` nebo exam variantami")
    return(meta)
  }

  if (grepl("/ggir$", dir_rel)) {
    meta$analysis <- "Vizualizační větev, kde GGIR parametry jsou outcome a klinické proměnné jsou kandidátní prediktory."
    meta$methods <- c(
      "Smíšené lineární modely přes návštěvy pacienta.",
      "Radarové skórování umožňuje porovnat, který klinický prediktor nejlépe vysvětluje daný GGIR parametr."
    )
    meta$scripts <- c(
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_04.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_04_exam.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_04_exam_b.R"
    )
    meta$inputs <- ggir_report_inputs
    meta$objects <- c("Objekt `d18_ggir_data`")
    return(meta)
  }

  if (grepl("/sf36$", dir_rel)) {
    meta$analysis <- "SF-36 workflow pro radarové porovnání různých modelových specifikací nad doménami kvality života."
    meta$methods <- c(
      "Pro každou SF-36 doménu vzniká trojice exportů: `fi2adj`, `mmt8adj` a `single`, případně exam/fi3 varianta.",
      "Každý export kombinuje tabulku s metrikami modelů a graf s normalizovaným radarovým skórováním."
    )
    meta$scripts <- c(
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_03.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_03_exam.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_03_exam_b.R"
    )
    meta$inputs <- ggir_report_inputs
    meta$objects <- c("Objekt `d18_ggir_data` připravený SF-36 helper funkcemi")
    return(meta)
  }

  if (grepl("/tis$", dir_rel)) {
    meta$analysis <- "TIS radar workflow navázané na spojená GGIR a TIS data."
    meta$methods <- c(
      "Spojení `d18_ggir_data` a `d18_tis_comp_01` přes `immet_id` a `exam`.",
      "Porovnání FI-2 adjusted, MMT8 adjusted a single modelů; v exam větvích navíc explicitní role návštěvy."
    )
    meta$scripts <- c(
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_05.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_05_exam.R",
      "scripts/Script_myo_2025_iim_clinics_ggir_vis_05_exam_b.R"
    )
    meta$inputs <- c(ggir_report_inputs, tis_inputs)
    meta$objects <- c("Objekty `d18_ggir_data` a `d18_tis_comp_01`")
    return(meta)
  }

  if (grepl("/figures$", dir_rel)) {
    meta$analysis <- "Tato podsložka obsahuje pouze obrazové exporty dané analytické větve."
    meta$audience <- "Nový čtenář zde najde vizuální interpretaci výsledků. AI agent by měl tuto složku číst společně se souběžnou `tables`, protože obrázky bývají shrnutím tabulkových metrik."
    meta$methods <- c(
      "Soubory jsou typicky generovány funkcemi `ggsave()`, `png()` nebo `tiff()`.",
      "Obrázky obvykle nevznikají samostatně; navazují na stejný objekt/model jako souběžná XLSX tabulka."
    )
    return(meta)
  }

  if (grepl("/tables$", dir_rel)) {
    meta$analysis <- "Tato podsložka obsahuje tabulkové exporty dané analytické větve."
    meta$audience <- "Pro klinika je to hlavní místo, kde lze dohledat přesné číselné výsledky za obrázky. Pro AI agenta jsou XLSX soubory nejpřímější stopou k reprodukci jednotlivých grafů."
    meta$methods <- c(
      "Tabulky jsou typicky exportované pomocí `rio::export()` nebo `export()` z balíku `rio`.",
      "U radarových workflow bývá jeden XLSX soubor svázán s jedním obrázkem stejného prefixu."
    )
    return(meta)
  }

  if (grepl("/ppt$", dir_rel)) {
    meta$analysis <- "Prezentační exporty odvozené z ceiling/floor workflow."
    meta$methods <- c(
      "Složka momentálně obsahuje PowerPoint určený pro komunikaci výsledků, nikoli primární analytický mezikrok.",
      "Obsah prezentace vychází z tabulek a obrázků ve stejném rodičovském adresáři."
    )
    meta$scripts <- c("scripts/Script_myo_2025_iim_ceilingEffect_01.R", "reports/cei_flo/*.Rmd")
    meta$inputs <- ceil_inputs()
    meta$objects <- c("Vizuální a tabulkové výstupy z rodičovské složky ceiling/floor")
    return(meta)
  }

  meta
}

refine_leaf_meta <- function(meta, dir_rel) {
  dir_rel <- normalize_rel(dir_rel)
  leaf <- basename(dir_rel)

  if (leaf == "figures") {
    meta$analysis <- paste(
      meta$analysis,
      "Tato konkrétní podsložka obsahuje pouze obrazové exporty dané větve, tedy grafy, heatmapy, radarové plochy, diagnostické obrázky nebo prezentační vizualizace."
    )
    meta$audience <- paste(
      meta$audience,
      "Nový čtenář by měl tuto složku číst společně se souběžnou `tables`, protože obrázek bývá interpretací stejného modelového výstupu."
    )
    meta$methods <- unique(c(
      meta$methods,
      "Obrázky vznikají typicky přes `ggsave()`, `png()` nebo `tiff()` a bývají navázané na stejné objektové výsledky jako XLSX soubory ve vedlejší složce `tables`."
    ))
  }

  if (leaf == "tables") {
    meta$analysis <- paste(
      meta$analysis,
      "Tato konkrétní podsložka obsahuje pouze tabulkové exporty, tedy číselné výsledky, souhrny metrik modelů, pomocné listy pro grafy a případně diagnostické přehledy."
    )
    meta$audience <- paste(
      meta$audience,
      "Pro klinika je tato podsložka důležitá proto, že zde leží přesná čísla stojící za obrázky."
    )
    meta$methods <- unique(c(
      meta$methods,
      "Tabulky jsou obvykle exportované pomocí `rio::export()` nebo `export()`; v radarových workflow bývá jeden XLSX soubor svázán s jedním obrázkem stejného prefixu."
    ))
  }

  if (leaf == "ppt") {
    meta$analysis <- paste(
      meta$analysis,
      "Tato podsložka je vyhrazená pro prezentační soubory a není primárním místem, kde vznikají analytické výsledky."
    )
    meta$methods <- unique(c(
      meta$methods,
      "Prezentace obvykle skládá již hotové grafy a tabulky z ostatních podsložek stejné analytické větve."
    ))
  }

  meta
}

build_readme <- function(dir_rel) {
  files <- direct_files(dir_rel)
  subdirs <- direct_subdirs(dir_rel)
  meta <- refine_leaf_meta(describe_dir(dir_rel), dir_rel)

  lines <- c(
    paste0("# ", meta$title),
    "",
    paste0("**Složka:** `", dir_rel, "`"),
    "",
    "## 1. Pro člověka",
    "",
    "### Co tato složka obsahuje",
    meta$analysis,
    "",
    "### Jak tuto složku číst",
    meta$audience,
    "",
    "### Přímý obsah složky",
    paste0("- počet přímých souborů: ", length(files)),
    paste0("- typy přímých souborů: ", count_extensions(files)),
    paste0("- počet přímých podsložek: ", length(subdirs)),
    "",
    "### Přímé podsložky",
    format_subdirs(subdirs),
    "",
    "### Příklady přímých souborů",
    format_files(files),
    "",
    "### Hlavní použité metody",
    paste0("- ", meta$methods, collapse = "\n"),
    "",
    "### Zdrojové skripty",
    format_paths(meta$scripts),
    "",
    "### Vstupní soubory a jejich umístění",
    format_paths(meta$inputs),
    "",
    "### Vstupní objekty a jak vznikají",
    format_paths(meta$objects),
    "",
    "### Jak číst názvy souborů",
    meta$naming,
    "",
    "## 2. Pro AI agenta",
    "",
    "### Doporučený reprodukční postup",
    paste0("1. ", meta$ai_steps[1]),
    if (length(meta$ai_steps) >= 2) paste0("2. ", meta$ai_steps[2]) else NULL,
    if (length(meta$ai_steps) >= 3) paste0("3. ", meta$ai_steps[3]) else NULL,
    "",
    "### Co zkontrolovat před opakováním analýzy",
    paste0("- ověř, že existují všechny vstupní cesty uvedené výše"),
    paste0("- ověř, že objekty v prostředí odpovídají očekávaným názvům ve zdrojových skriptech"),
    paste0("- porovnej naming nově vytvořených souborů s existujícími exporty v této složce"),
    "",
    "### Rizika a nejasnosti",
    if (length(meta$caveats) == 0) "- bez zvláštních poznámek nad rámec výše uvedeného" else paste0("- ", meta$caveats, collapse = "\n"),
    ""
  )

  enc2utf8(lines)
}

for (dir_rel in dir_paths) {
  readme_path <- file.path(dir_rel, "README.md")
  writeLines(build_readme(dir_rel), readme_path, useBytes = TRUE)
}
