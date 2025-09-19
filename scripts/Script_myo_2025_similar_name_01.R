library(tidyverse)
library(janitor)
library(stringi)
library(fuzzyjoin)

# Example objects
# df  ... your data frame with already-"cleaned" column names (e.g., snake_case)
# full_names ... a character vector of the verbose/full variable names

# 1) One consistent normalization for BOTH sides
norm <- function(x) {
  x |> 
    stri_trans_general("Latin-ASCII") |>      # drop accents (č -> c)
    str_to_lower() |> 
    str_replace_all("[^a-z0-9]+", "_") |>     # collapse to underscores
    str_replace_all("^_|_$", "")               # trim edge underscores
}

full_names <- c("odpověď na terapii M0 vs M6", "Mi-2", "TIF-1γ", "MDA-5", "SAE", "NXP-2", "SRP", "Jo-1", "PM-Scl", "RNP", "Ku", "Ro52", "Ro60", "La", "PL-12", "Pl-7", "Ej", "cN-1A", "anti HMGCR", "Muscle weakness", "Oesophageus motiliy disorder", "Skin rash", "Mechanic´s hands", "Raynaud phenomenon", "Arthritis", "ILD", "Cardiac involvement", "denná dávka GK (mg/kg BW)", "GK", "MTX", "AZA", "CSA", "CPA", "LEF", "MFF", "Tacrolimus", "SAS", "PQ", "RTX", "ABA", "IVIg", " Patient activity VAS", "Physician activity VAS", " Constitutional Disease Activity", "Cutaneous Disease Activity", "Skeletal Disease Activity", "Gastrointestinal Disease Activity", "Pulmonary Disease Activity", "Cardiovascular Disease Activity", "Muscle Disease Activity", "Extramuscular Global Assessment", "Podtyp nemoci_zjednodušený", "Statiny", "DM", "PAD", "Inzulin (terapie)", "Věk", "Pohlaví")


cols <- tibble(col = names(d03), key = norm(names(d03)))
full <- tibble(full_name = full_names,
               key = norm(full_names))

# 2) First try exact matches on the normalized key
m1 <- cols |>
  left_join(full, by = "key")                  # full_name will be NA if no direct match

# 3) Fuzzy match anything still unmatched (Jaro-Winker distance)
todo <- m1 |> filter(is.na(full_name)) |> select(col, key)
cand <- full |> select(full_name, key)

m2 <- todo |>
  stringdist_left_join(cand,
                       by = c("key" = "key"),
                       method = "jw",          # Jaro–Winkler works well for names
                       max_dist = 0.15,        # tighten/loosen as needed (0.10–0.20)
                       distance_col = "dist") |>
  group_by(col) |>
  slice_min(dist, n = 1, with_ties = FALSE) |>
  ungroup()

# 4) Final mapping table (one full_name per column)
map <- m1 |>
  filter(!is.na(full_name)) |>
  bind_rows(m2 |> select(col, full_name)) |>
  distinct(col, .keep_all = TRUE)

# 5) Inspect anything still unmatched
unmatched <- setdiff(names(df), map$col)
unmatched
# -> review these manually or add specific overrides

# 6) (Optional) Rename df's columns to the full names
#    dplyr::rename expects new = old, so build that mapping:
rename_map <- setNames(map$col, make.unique(map$full_name))  # make.unique guards against duplicates
df_renamed <- dplyr::rename(df, !!!rename_map)
