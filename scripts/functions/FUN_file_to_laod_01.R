file_to_load <- function(folder_path, phrase) {
  # Získání všech souborů (včetně cesty)
  files <- list.files(folder_path, recursive = TRUE, full.names = TRUE)
  
  # Filtrování souborů podle fráze v názvu
  matching_files <- files[grepl(phrase, basename(files), ignore.case = TRUE)]
  
  return(matching_files)
}
