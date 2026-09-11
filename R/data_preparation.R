# Data preparation and provenance helpers -------------------------------------

required_packages(c("readxl", "dplyr", "tidyr", "purrr", "stringr", "tibble"))

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
})

CORE_AREA_M2 <- 0.002827

resolve_excel_file <- function(directory, stem) {
  candidates <- file.path(directory, paste0(stem, c(".xlsx", ".xls")))
  existing <- candidates[file.exists(candidates)]

  if (!length(existing)) {
    stop(
      "No Excel file found for '", stem, "'. Expected one of:\n",
      paste(" -", candidates, collapse = "\n")
    )
  }
  if (length(existing) > 1L) {
    warning(
      "Both .xlsx and .xls files found for '", stem,
      "'. Using: ", basename(existing[[1]])
    )
  }
  existing[[1]]
}

REFERENCE_FILES <- c(
  RMQS_2024 = resolve_excel_file(RAW_DIR, "RMQS_2024_COLLEMBOLA"),
  RMQS_2021 = resolve_excel_file(RAW_DIR, "RMQS_2021_COLLEMBOLA"),
  RMQS_Bretagne = resolve_excel_file(RAW_DIR, "RMQS_Bretagne"),
  ANDRA = resolve_excel_file(RAW_DIR, "ANDRA_ARTHROPODA"),
  TIGA_rural = resolve_excel_file(RAW_DIR, "TIGA_rural_MESOFAUNA")
)

BIOINDICATEUR_FILE <- resolve_excel_file(RAW_DIR, "Bioindicateur2_ARTHROPODA")
TAXREF_FILE <- file.path(RAW_DIR, "TaxRef18_Collembola.csv")
HISTORICAL_MERGED_FILE <- file.path(AUDIT_DATA_DIR, "Stage_M2_Helio_Suarez.xlsx")

FIXED_REFERENCE_COLUMNS <- c(
  "Site", "Projet", "Strategie", "Echelle", "Annee",
  "X_WGS84", "Y_WGS84", "X_L93", "Y_L93",
  "CLC_niveau1", "CLC_niveau2", "CLC_niveau3"
)

EXPECTED_ALL_COUNTS <- c(
  RMQS_Biodiversite = 96L,
  RMQS_BioDiv_Bretagne = 98L,
  ANDRA = 132L,
  TIGA_rural = 430L
)

EXPECTED_AGRICULTURAL_COUNTS <- c(
  RMQS_Biodiversite = 58L,
  RMQS_BioDiv_Bretagne = 89L,
  ANDRA = 90L,
  TIGA_rural = 345L
)

read_excel_checked <- function(path, sheet, ...) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  sheets <- readxl::excel_sheets(path)
  if (!sheet %in% sheets) {
    stop(
      "Sheet '", sheet, "' not found in ", basename(path),
      ". Available sheets: ", paste(sheets, collapse = ", ")
    )
  }
  readxl::read_excel(path, sheet = sheet, ...)
}

assert_fixed_reference_columns <- function(df, label = "dataset") {
  missing <- setdiff(FIXED_REFERENCE_COLUMNS, names(df))
  if (length(missing)) {
    stop(
      label, " is missing fixed reference columns: ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}

read_curated_reference_sheets <- function(files = REFERENCE_FILES) {
  assert_files_exist(files, context = "Reference workbook")
  out <- lapply(names(files), function(nm) {
    dat <- read_excel_checked(files[[nm]], "Feuille stats")
    assert_fixed_reference_columns(dat, label = nm)
    dat$.source_workbook <- basename(files[[nm]])
    dat$.source_dataset <- nm
    dat
  })
  names(out) <- names(files)
  out
}

harmonise_reference_sheets <- function(datasets) {
  # This reproduces the final compilation block in Helio's historical
  # `Mise en forme.R`: fixed metadata are cast to character, each dataset is
  # pivoted long, datasets are row-bound, then taxa are pivoted back to columns.
  stripped <- lapply(datasets, function(df) {
    df <- as.data.frame(df)
    df$.source_workbook <- NULL
    df$.source_dataset <- NULL
    assert_fixed_reference_columns(df)
    df %>% mutate(across(all_of(FIXED_REFERENCE_COLUMNS), as.character))
  })

  long <- lapply(stripped, function(df) {
    df %>%
      pivot_longer(
        cols = -all_of(FIXED_REFERENCE_COLUMNS),
        names_to = "Espece",
        values_to = "Abondance"
      )
  })

  bind_rows(long) %>%
    pivot_wider(
      names_from = "Espece",
      values_from = "Abondance",
      values_fill = NA
    )
}

read_bioindicateur_curated <- function(path = BIOINDICATEUR_FILE) {
  assert_files_exist(path, context = "Bioindicateur2 workbook")
  read_excel_checked(path, "Feuille stats")
}

reference_count_table <- function(df) {
  if (!all(c("Projet", "CLC_niveau1") %in% names(df))) {
    stop("Projet and CLC_niveau1 are required for count audit.")
  }
  all_counts <- df %>% count(.data$Projet, name = "n_all")
  ag_counts <- df %>%
    filter(as.character(.data$CLC_niveau1) == "2") %>%
    count(.data$Projet, name = "n_agricultural")

  full_join(all_counts, ag_counts, by = "Projet") %>%
    arrange(.data$Projet)
}

assert_expected_reference_counts <- function(df, strict = TRUE) {
  counts <- reference_count_table(df)
  obs_all <- setNames(counts$n_all, counts$Projet)
  obs_ag <- setNames(counts$n_agricultural, counts$Projet)

  problems <- character()
  for (nm in names(EXPECTED_ALL_COUNTS)) {
    if (is.na(obs_all[[nm]]) || obs_all[[nm]] != EXPECTED_ALL_COUNTS[[nm]]) {
      problems <- c(problems, sprintf("%s all sites: observed %s, expected %s", nm, obs_all[[nm]], EXPECTED_ALL_COUNTS[[nm]]))
    }
    if (is.na(obs_ag[[nm]]) || obs_ag[[nm]] != EXPECTED_AGRICULTURAL_COUNTS[[nm]]) {
      problems <- c(problems, sprintf("%s agricultural sites: observed %s, expected %s", nm, obs_ag[[nm]], EXPECTED_AGRICULTURAL_COUNTS[[nm]]))
    }
  }

  if (length(problems)) {
    msg <- paste("Reference count audit failed:", paste0(" - ", problems, collapse = "\n"), sep = "\n")
    if (strict) stop(msg) else warning(msg)
  }
  invisible(counts)
}

row_key <- function(df, preferred = c("Projet", "Site", "Annee")) {
  cols <- preferred[preferred %in% names(df)]
  if (!length(cols)) return(as.character(seq_len(nrow(df))))
  do.call(paste, c(lapply(df[cols], function(x) ifelse(is.na(x), "<NA>", as.character(x))), sep = "||"))
}

compare_merged_to_historical <- function(rebuilt, historical) {
  # Compare by Projet/Site/Annee rather than row order, because the long/wide
  # reconstruction may legitimately reorder taxon columns.
  key_a <- row_key(rebuilt)
  key_b <- row_key(historical)

  duplicate_a <- sum(duplicated(key_a))
  duplicate_b <- sum(duplicated(key_b))
  common_cols <- intersect(names(rebuilt), names(historical))
  missing_in_rebuilt <- setdiff(names(historical), names(rebuilt))
  extra_in_rebuilt <- setdiff(names(rebuilt), names(historical))

  common_keys <- intersect(key_a, key_b)
  idx_a <- match(common_keys, key_a)
  idx_b <- match(common_keys, key_b)

  metadata_cols <- intersect(FIXED_REFERENCE_COLUMNS, common_cols)
  numeric_cols <- setdiff(common_cols, metadata_cols)

  metadata_mismatches <- 0L
  if (length(common_keys) && length(metadata_cols)) {
    for (nm in metadata_cols) {
      aa <- ifelse(is.na(rebuilt[[nm]][idx_a]), "<NA>", as.character(rebuilt[[nm]][idx_a]))
      bb <- ifelse(is.na(historical[[nm]][idx_b]), "<NA>", as.character(historical[[nm]][idx_b]))
      metadata_mismatches <- metadata_mismatches + sum(aa != bb)
    }
  }

  numeric_mismatches <- 0L
  max_abs_diff <- 0
  compared_numeric_cells <- 0L
  if (length(common_keys) && length(numeric_cols)) {
    for (nm in numeric_cols) {
      aa <- suppressWarnings(as.numeric(rebuilt[[nm]][idx_a]))
      bb <- suppressWarnings(as.numeric(historical[[nm]][idx_b]))
      both_na <- is.na(aa) & is.na(bb)
      one_na <- xor(is.na(aa), is.na(bb))
      finite <- is.finite(aa) & is.finite(bb)
      diff <- rep(0, length(aa))
      diff[finite] <- abs(aa[finite] - bb[finite])
      mismatch <- one_na | (finite & diff > 1e-8)
      numeric_mismatches <- numeric_mismatches + sum(mismatch)
      if (any(finite)) max_abs_diff <- max(max_abs_diff, diff[finite], na.rm = TRUE)
      compared_numeric_cells <- compared_numeric_cells + sum(!both_na)
    }
  }

  tibble(
    metric = c(
      "rebuilt_rows", "historical_rows", "rebuilt_columns", "historical_columns",
      "duplicate_rebuilt_keys", "duplicate_historical_keys",
      "common_keys", "missing_keys_in_rebuilt", "extra_keys_in_rebuilt",
      "missing_columns_in_rebuilt", "extra_columns_in_rebuilt",
      "metadata_cell_mismatches", "numeric_cell_mismatches",
      "compared_numeric_cells", "max_abs_numeric_difference"
    ),
    value = c(
      nrow(rebuilt), nrow(historical), ncol(rebuilt), ncol(historical),
      duplicate_a, duplicate_b,
      length(common_keys), length(setdiff(key_b, key_a)), length(setdiff(key_a, key_b)),
      length(missing_in_rebuilt), length(extra_in_rebuilt),
      metadata_mismatches, numeric_mismatches,
      compared_numeric_cells, max_abs_diff
    ),
    detail = c(
      "", "", "", "", "", "", "", "", "",
      paste(missing_in_rebuilt, collapse = " | "),
      paste(extra_in_rebuilt, collapse = " | "),
      "", "", "", ""
    )
  )
}

# ---- Upstream reconstruction audit -------------------------------------------
# These functions reproduce only transformations explicitly present in the
# historical Mise en forme.R. They deliberately do NOT guess undocumented Excel
# edits or taxonomic corrections.

rename_if_present <- function(df, mapping) {
  for (old in names(mapping)) {
    new <- mapping[[old]]
    if (old %in% names(df) && !new %in% names(df)) names(df)[names(df) == old] <- new
  }
  df
}

prepare_andra_from_temporary <- function(path = REFERENCE_FILES[["ANDRA"]]) {
  d <- read_excel_checked(path, "temporary")
  d <- rename_if_present(d, c("Reference Site" = "Site", "Reference_Site" = "Site", "Année" = "Annee"))
  if (!all(c("Site", "Annee") %in% names(d))) stop("ANDRA temporary sheet lacks Site/Annee after conservative renaming.")
  d <- d %>%
    group_by(.data$Site, .data$Annee) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  scale_cols <- setdiff(names(d)[vapply(d, is.numeric, logical(1))], "Annee")
  d %>% mutate(across(all_of(scale_cols), ~ .x / CORE_AREA_M2))
}

prepare_tiga_from_temporary <- function(path = REFERENCE_FILES[["TIGA_rural"]]) {
  # The preserved workbook has formatting extending to Excel row 1,048,576;
  # n_max prevents readxl from materialising an enormous blank tail.
  d <- read_excel_checked(path, "temporary", n_max = 10000)
  if (!"Site" %in% names(d)) stop("TIGA temporary sheet lacks Site.")
  d <- d %>%
    mutate(Site = stringr::str_remove(as.character(.data$Site), "_[1-3]$")) %>%
    group_by(.data$Site) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(across(where(is.numeric), ~ .x / CORE_AREA_M2))
  d
}

prepare_bretagne_from_temporary <- function(path = REFERENCE_FILES[["RMQS_Bretagne"]]) {
  d <- read_excel_checked(path, "temporary")
  if (!"Site" %in% names(d)) stop("RMQS Bretagne temporary sheet lacks Site.")
  d %>%
    mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0))) %>%
    group_by(.data$Site) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(across(where(is.numeric), ~ .x / CORE_AREA_M2))
}

prepare_bioindicateur_from_temporary <- function(path = BIOINDICATEUR_FILE) {
  d <- read_excel_checked(path, "temporary")
  needed <- c("Site", "Annee", "Nom_site", "Usage", "Traitement", "Contamination")
  d <- rename_if_present(d, c("Année" = "Annee"))
  if (!all(needed %in% names(d))) {
    stop("Bioindicateur2 temporary sheet lacks: ", paste(setdiff(needed, names(d)), collapse = ", "))
  }
  d$Site <- substr(as.character(d$Site), 1, pmax(nchar(as.character(d$Site)) - 1, 0))
  d <- d %>%
    group_by(across(all_of(needed))) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  scale_cols <- setdiff(names(d)[vapply(d, is.numeric, logical(1))], "Annee")
  d %>% mutate(across(all_of(scale_cols), ~ .x / CORE_AREA_M2))
}

get_taxon_columns_curated <- function(df) {
  anchor <- match("CLC_niveau3", names(df))
  if (is.na(anchor) || anchor >= ncol(df)) return(character())
  names(df)[seq.int(anchor + 1L, ncol(df))]
}

compare_preparation_to_curated <- function(rebuilt, curated, key_cols, label) {
  key_cols <- key_cols[key_cols %in% names(rebuilt) & key_cols %in% names(curated)]
  if (!length(key_cols)) {
    return(tibble(dataset = label, status = "no_common_key", detail = "No requested key columns available"))
  }

  key_a <- do.call(paste, c(lapply(rebuilt[key_cols], as.character), sep = "||"))
  key_b <- do.call(paste, c(lapply(curated[key_cols], as.character), sep = "||"))
  common <- intersect(key_a, key_b)
  idx_a <- match(common, key_a)
  idx_b <- match(common, key_b)

  curated_taxa <- get_taxon_columns_curated(curated)
  candidate_taxa <- intersect(curated_taxa, names(rebuilt))

  mismatch <- 0L
  compared <- 0L
  max_diff <- 0
  for (nm in candidate_taxa) {
    aa <- suppressWarnings(as.numeric(rebuilt[[nm]][idx_a]))
    bb <- suppressWarnings(as.numeric(curated[[nm]][idx_b]))
    one_na <- xor(is.na(aa), is.na(bb))
    finite <- is.finite(aa) & is.finite(bb)
    dd <- rep(0, length(aa))
    dd[finite] <- abs(aa[finite] - bb[finite])
    mismatch <- mismatch + sum(one_na | (finite & dd > 1e-8))
    compared <- compared + sum(one_na | finite)
    if (any(finite)) max_diff <- max(max_diff, dd[finite], na.rm = TRUE)
  }

  status <- if (length(common) == nrow(curated) && mismatch == 0L) "matches_for_shared_taxa" else "differences_detected"

  tibble(
    dataset = label,
    status = status,
    rebuilt_rows = nrow(rebuilt),
    curated_rows = nrow(curated),
    common_keys = length(common),
    curated_keys_missing_from_rebuilt = length(setdiff(key_b, key_a)),
    rebuilt_keys_missing_from_curated = length(setdiff(key_a, key_b)),
    shared_taxa_compared = length(candidate_taxa),
    compared_cells = compared,
    mismatched_cells = mismatch,
    max_abs_difference = max_diff,
    detail = if (length(setdiff(curated_taxa, names(rebuilt)))) {
      paste0("Curated taxa absent from reconstructed temporary sheet: ", length(setdiff(curated_taxa, names(rebuilt))))
    } else ""
  )
}

run_upstream_preparation_audit <- function(curated_refs, curated_bio) {
  rows <- list()

  rows[["RMQS_2024"]] <- tibble(
    dataset = "RMQS_2024",
    status = "not_reconstructable_from_preserved_script",
    rebuilt_rows = NA_integer_, curated_rows = nrow(curated_refs[["RMQS_2024"]]),
    common_keys = NA_integer_, curated_keys_missing_from_rebuilt = NA_integer_, rebuilt_keys_missing_from_curated = NA_integer_,
    shared_taxa_compared = NA_integer_, compared_cells = NA_integer_, mismatched_cells = NA_integer_, max_abs_difference = NA_real_,
    detail = "Historical script begins from a Feuille stats sheet that has since been overwritten by its own derived output; upstream manual Excel edits are not recoverable from the script alone."
  )

  rows[["RMQS_2021"]] <- tibble(
    dataset = "RMQS_2021",
    status = "curated_input_only",
    rebuilt_rows = NA_integer_, curated_rows = nrow(curated_refs[["RMQS_2021"]]),
    common_keys = NA_integer_, curated_keys_missing_from_rebuilt = NA_integer_, rebuilt_keys_missing_from_curated = NA_integer_,
    shared_taxa_compared = NA_integer_, compared_cells = NA_integer_, mismatched_cells = NA_integer_, max_abs_difference = NA_real_,
    detail = "Mise en forme.R contains coordinate conversion for RMQS 2021 but no preserved block reconstructing Feuille stats biological densities from raw counts."
  )

  attempts <- list(
    ANDRA = list(fun = prepare_andra_from_temporary, curated = curated_refs[["ANDRA"]], key = c("Site", "Annee")),
    TIGA_rural = list(fun = prepare_tiga_from_temporary, curated = curated_refs[["TIGA_rural"]], key = "Site"),
    RMQS_Bretagne = list(fun = prepare_bretagne_from_temporary, curated = curated_refs[["RMQS_Bretagne"]], key = "Site"),
    Bioindicateur2 = list(fun = prepare_bioindicateur_from_temporary, curated = curated_bio, key = "Site")
  )

  for (nm in names(attempts)) {
    a <- attempts[[nm]]
    rows[[nm]] <- tryCatch(
      compare_preparation_to_curated(a$fun(), a$curated, a$key, nm),
      error = function(e) tibble(
        dataset = nm, status = "audit_error",
        rebuilt_rows = NA_integer_, curated_rows = nrow(a$curated),
        common_keys = NA_integer_, curated_keys_missing_from_rebuilt = NA_integer_, rebuilt_keys_missing_from_curated = NA_integer_,
        shared_taxa_compared = NA_integer_, compared_cells = NA_integer_, mismatched_cells = NA_integer_, max_abs_difference = NA_real_,
        detail = conditionMessage(e)
      )
    )
  }

  bind_rows(rows)
}

# ---- TAXREF audit -------------------------------------------------------------

historical_taxref_check <- function(df, taxref, dataset_label) {
  taxa <- get_taxon_columns_curated(df)
  if (!length(taxa)) {
    return(tibble(dataset = dataset_label, taxon = character(), status = character()))
  }
  lb <- iconv(as.character(taxref$LB_NOM), from = "", to = "UTF-8", sub = "byte")
  genus <- sub(" .*", "", lb)
  direct <- taxa %in% lb
  genus_sp <- grepl("^[A-Z][a-z]+ sp\\.$", taxa)
  genus_name <- sub(" sp\\.$", "", taxa)
  valid_genus_sp <- genus_sp & genus_name %in% genus
  status <- ifelse(direct, "direct_taxref_match", ifelse(valid_genus_sp, "valid_genus_sp", "not_matched_by_historical_rule"))
  tibble(dataset = dataset_label, taxon = taxa, status = status)
}

run_taxref_audit <- function(curated_refs, curated_bio, taxref_path = TAXREF_FILE) {
  if (!file.exists(taxref_path)) {
    return(tibble(dataset = NA_character_, taxon = NA_character_, status = "TaxRef file absent"))
  }
  taxref <- utils::read.csv(taxref_path, sep = ";", stringsAsFactors = FALSE, check.names = FALSE)
  if (!"LB_NOM" %in% names(taxref)) stop("TaxRef file does not contain LB_NOM.")
  bind_rows(
    lapply(names(curated_refs), function(nm) historical_taxref_check(curated_refs[[nm]], taxref, nm)),
    list(historical_taxref_check(curated_bio, taxref, "Bioindicateur2"))
  )
}
