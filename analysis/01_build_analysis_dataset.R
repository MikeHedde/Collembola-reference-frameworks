#!/usr/bin/env Rscript

# ==============================================================================
# 01 - Build analysis-ready datasets from the separate curated workbooks
# ==============================================================================
# Stage_M2_Helio_Suarez.xlsx is NOT used as input. The five current curated
# Feuille stats sheets are harmonised and merged programmatically.
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

message("[01] Building analysis dataset from separate curated workbooks...")

assert_files_exist(c(REFERENCE_FILES, BIOINDICATEUR_FILE), context = "Core curated input")

curated_refs <- read_curated_reference_sheets()
reference_merged <- harmonise_reference_sheets(curated_refs)
validation_curated <- read_bioindicateur_curated()

# Guardrails based on the manuscript dataset and the supplied historical merge.
counts <- assert_expected_reference_counts(reference_merged, strict = TRUE)
write_csv_safe(counts, file.path(DERIVED_DIR, "reference_counts.csv"))

saveRDS(reference_merged, file.path(DERIVED_DIR, "reference_frameworks_merged.rds"))
saveRDS(validation_curated, file.path(DERIVED_DIR, "bioindicateur2_curated.rds"))

# Small text manifest: avoids writing another opaque Excel file as a required
# intermediate while making provenance explicit.
manifest <- c(
  paste0("Built: ", format(Sys.time(), tz = "UTC"), " UTC"),
  "Reference source sheets:",
  paste0(" - ", names(REFERENCE_FILES), ": ", basename(REFERENCE_FILES), " [Feuille stats]"),
  paste0("Validation source: ", basename(BIOINDICATEUR_FILE), " [Feuille stats]"),
  "Historical Stage_M2_Helio_Suarez.xlsx used as input: NO",
  paste0("Merged reference rows: ", nrow(reference_merged)),
  paste0("Merged reference columns: ", ncol(reference_merged))
)
write_text_lines(manifest, file.path(DERIVED_DIR, "build_manifest.txt"))

# Optional audit against the historical merged workbook, never a dependency.
if (file.exists(HISTORICAL_MERGED_FILE)) {
  historical <- read_excel_checked(HISTORICAL_MERGED_FILE, "Data")
  comparison <- compare_merged_to_historical(reference_merged, historical)
  write_csv_safe(comparison, file.path(DERIVED_DIR, "historical_merge_comparison.csv"))
}

message("[01] Built:")
message("     ", file.path(DERIVED_DIR, "reference_frameworks_merged.rds"))
message("     ", file.path(DERIVED_DIR, "bioindicateur2_curated.rds"))
