#!/usr/bin/env Rscript

# ==============================================================================
# 00 - Audit provenance and reconstructability of H. Suarez's input preparation
# ==============================================================================
# This script answers two separate questions:
#   A. Can the historical merged file Stage_M2_Helio_Suarez.xlsx be recreated
#      from the current curated "Feuille stats" sheets without using it as input?
#   B. How much of each curated "Feuille stats" can be recreated from the
#      upstream sheets using only transformations explicitly recorded in
#      legacy/Mise_en_forme_Helio_2025.R?
#
# A is required for analysis-level reproducibility. B is a provenance audit and
# is expected to expose undocumented/manual Excel curation where it exists.
# ============================================================================== 

# Locate the project even when Rscript is launched outside the repository root.
.project_paths_candidates <- c(file.path(getwd(), "R", "project_paths.R"))
.args_all <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("^--file=", .args_all, value = TRUE)
if (length(.file_arg)) {
  .script_file <- normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)
  .project_paths_candidates <- c(.project_paths_candidates, file.path(dirname(.script_file), "..", "R", "project_paths.R"))
}
.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (!is.null(.ofile)) {
  .project_paths_candidates <- c(.project_paths_candidates, file.path(dirname(normalizePath(.ofile, mustWork = FALSE)), "..", "R", "project_paths.R"))
}
.project_paths_candidates <- unique(normalizePath(.project_paths_candidates, mustWork = FALSE))
.project_paths_file <- .project_paths_candidates[file.exists(.project_paths_candidates)][1]
if (is.na(.project_paths_file)) stop("Could not locate R/project_paths.R. Open the .Rproj or run inside the repository.")
source(.project_paths_file)
rm(.project_paths_candidates, .args_all, .file_arg, .ofile, .project_paths_file)
if (exists(".script_file")) rm(.script_file)
source_repo("R/io_helpers.R")
source_repo("R/data_preparation.R")

message("[00] Auditing supplied workbooks and historical preprocessing...")

assert_files_exist(c(REFERENCE_FILES, BIOINDICATEUR_FILE), context = "Core curated input")

# ---- 1. Workbook inventory ----------------------------------------------------
workbook_paths <- c(REFERENCE_FILES, Bioindicateur2 = BIOINDICATEUR_FILE)
workbook_inventory <- bind_rows(lapply(names(workbook_paths), function(nm) {
  path <- workbook_paths[[nm]]
  sheets <- readxl::excel_sheets(path)
  tibble(
    dataset = nm,
    file = basename(path),
    sheet = sheets,
    has_feuille_stats = "Feuille stats" %in% sheets,
    has_temporary = "temporary" %in% sheets
  )
}))
write_csv_safe(workbook_inventory, file.path(AUDIT_OUTPUT_DIR, "audit_workbook_inventory.csv"))

input_paths <- c(workbook_paths, TaxRef = TAXREF_FILE,
                 historical_merged = if (file.exists(HISTORICAL_MERGED_FILE)) HISTORICAL_MERGED_FILE else character())
input_checksums <- tibble(
  label = names(input_paths),
  file = basename(input_paths),
  bytes = as.numeric(file.info(input_paths)$size),
  md5 = unname(tools::md5sum(input_paths))
)
write_csv_safe(input_checksums, file.path(AUDIT_OUTPUT_DIR, "audit_input_checksums.csv"))

# ---- 2. Read current curated inputs ------------------------------------------
curated_refs <- read_curated_reference_sheets()
curated_bio <- read_bioindicateur_curated()

curated_summary <- bind_rows(
  lapply(names(curated_refs), function(nm) {
    d <- curated_refs[[nm]]
    tibble(
      dataset = nm,
      rows = nrow(d),
      columns = ncol(d),
      taxon_columns_after_CLC3 = length(get_taxon_columns_curated(d)),
      project_labels = paste(unique(as.character(d$Projet)), collapse = " | ")
    )
  }),
  list(tibble(
    dataset = "Bioindicateur2",
    rows = nrow(curated_bio),
    columns = ncol(curated_bio),
    taxon_columns_after_CLC3 = length(get_taxon_columns_curated(curated_bio)),
    project_labels = if ("Projet" %in% names(curated_bio)) paste(unique(as.character(curated_bio$Projet)), collapse = " | ") else NA_character_
  ))
)
write_csv_safe(curated_summary, file.path(AUDIT_OUTPUT_DIR, "audit_curated_input_summary.csv"))

# ---- 3. Rebuild merged reference dataset from the five curated sheets --------
rebuilt <- harmonise_reference_sheets(curated_refs)
count_table <- reference_count_table(rebuilt)
write_csv_safe(count_table, file.path(AUDIT_OUTPUT_DIR, "audit_rebuilt_reference_counts.csv"))
assert_expected_reference_counts(rebuilt, strict = FALSE)

saveRDS(rebuilt, file.path(AUDIT_OUTPUT_DIR, "rebuilt_reference_from_curated_sheets.rds"))

# Compare to historical merged workbook when available. Crucially, the
# historical workbook is an audit target only, never an input to reconstruction.
if (file.exists(HISTORICAL_MERGED_FILE)) {
  historical <- read_excel_checked(HISTORICAL_MERGED_FILE, "Data")
  merged_compare <- compare_merged_to_historical(rebuilt, historical)
  write_csv_safe(merged_compare, file.path(AUDIT_OUTPUT_DIR, "audit_merged_vs_historical.csv"))

  # Compact machine-readable verdict.
  vals <- setNames(merged_compare$value, merged_compare$metric)
  exact_content_match <-
    vals[["missing_keys_in_rebuilt"]] == 0 &&
    vals[["extra_keys_in_rebuilt"]] == 0 &&
    vals[["missing_columns_in_rebuilt"]] == 0 &&
    vals[["extra_columns_in_rebuilt"]] == 0 &&
    vals[["metadata_cell_mismatches"]] == 0 &&
    vals[["numeric_cell_mismatches"]] == 0

  write_csv_safe(
    tibble(
      audit = "curated_sheets_to_historical_merged",
      status = if (exact_content_match) "PASS" else "CHECK_DIFFERENCES",
      historical_file_used_as_analysis_input = FALSE,
      note = "Column order is ignored; values are compared after matching Projet/Site/Annee keys."
    ),
    file.path(AUDIT_OUTPUT_DIR, "audit_merged_verdict.csv")
  )
} else {
  write_csv_safe(
    tibble(
      audit = "curated_sheets_to_historical_merged",
      status = "NOT_RUN_HISTORICAL_FILE_ABSENT",
      historical_file_used_as_analysis_input = FALSE,
      note = "Stage_M2_Helio_Suarez.xlsx is optional and used only as an audit target."
    ),
    file.path(AUDIT_OUTPUT_DIR, "audit_merged_verdict.csv")
  )
}

# ---- 4. Attempt upstream reconstruction from historical 'temporary' sheets ---
upstream <- run_upstream_preparation_audit(curated_refs, curated_bio)
write_csv_safe(upstream, file.path(AUDIT_OUTPUT_DIR, "audit_upstream_preparation.csv"))

# ---- 5. Reproduce historical TAXREF QC rule ----------------------------------
taxref_audit <- run_taxref_audit(curated_refs, curated_bio)
write_csv_safe(taxref_audit, file.path(AUDIT_OUTPUT_DIR, "audit_taxref_all_taxa.csv"))
write_csv_safe(
  taxref_audit %>% filter(.data$status == "not_matched_by_historical_rule"),
  file.path(AUDIT_OUTPUT_DIR, "audit_taxref_unmatched.csv")
)

# ---- 6. Provenance summary ----------------------------------------------------
provenance <- tibble(
  layer = c(
    "Raw/original field-lab exports",
    "Curated workbook sheets (Feuille stats)",
    "Merged historical Stage_M2_Helio_Suarez.xlsx",
    "Revised analysis dataset"
  ),
  role = c(
    "Upstream source where preserved; historical script documents only part of the transformation and mentions manual Excel edits.",
    "Primary reproducibility input for the revision workflow.",
    "Historical audit target only; no longer required as an input.",
    "Rebuilt automatically by analysis/01_build_analysis_dataset.R from curated Feuille stats sheets."
  ),
  reproducibility_status = c(
    "PARTIAL / dataset-specific",
    "SUPPORTED",
    "DERIVED / AUDIT ONLY",
    "FULLY SCRIPTED FROM CURATED INPUTS"
  )
)
write_csv_safe(provenance, file.path(AUDIT_OUTPUT_DIR, "audit_provenance_layers.csv"))

capture.output(sessionInfo(), file = file.path(AUDIT_OUTPUT_DIR, "sessionInfo_audit.txt"))
message("[00] Audit complete: ", AUDIT_OUTPUT_DIR)
