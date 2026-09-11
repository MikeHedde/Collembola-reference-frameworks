#!/usr/bin/env Rscript

# Run complete audit -> build -> baseline -> reviewer -> manuscript -> supplementary workflow.
# Recommended use: open the .Rproj, then source("analysis/run_all.R").

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

steps <- c(
  "analysis/00_audit_input_reconstruction.R",
  "analysis/01_build_analysis_dataset.R",
  "analysis/02_reproduce_original_results.R",
  "analysis/reviewer_robustness_analysis.R",
  "analysis/03_regenerate_manuscript_figures.R",
  "analysis/04_regenerate_supplementary_material.R"
)

for (step in steps) {
  message("\n============================================================")
  message("Running ", step)
  message("============================================================")
  source(file.path(ROOT_DIR, step), chdir = FALSE)
}

message("\nComplete workflow finished.")
