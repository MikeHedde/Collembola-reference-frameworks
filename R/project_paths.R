# Project path helpers ---------------------------------------------------------

find_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  starts <- c(
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)) else character(),
    getwd()
  )
  starts <- unique(normalizePath(starts, mustWork = FALSE))

  for (start in starts) {
    current <- start
    repeat {
      if (file.exists(file.path(current, "collembola-reference-frameworks.Rproj"))) {
        return(normalizePath(current, mustWork = TRUE))
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }

  stop(
    "Could not identify repository root. Open collembola-reference-frameworks.Rproj ",
    "or run R/Rscript from inside the repository."
  )
}

ROOT_DIR <- find_repo_root()
DATA_DIR <- file.path(ROOT_DIR, "data")
RAW_DIR <- file.path(DATA_DIR, "raw")
AUDIT_DATA_DIR <- file.path(DATA_DIR, "audit")
DERIVED_DIR <- file.path(DATA_DIR, "derived")
OUTPUT_DIR <- file.path(ROOT_DIR, "output")
AUDIT_OUTPUT_DIR <- file.path(OUTPUT_DIR, "audit")
BASELINE_OUTPUT_DIR <- file.path(OUTPUT_DIR, "baseline")
REVIEWER_OUTPUT_DIR <- file.path(OUTPUT_DIR, "reviewer")
MANUSCRIPT_OUTPUT_DIR <- file.path(OUTPUT_DIR, "manuscript")
SUPPLEMENTARY_OUTPUT_DIR <- file.path(OUTPUT_DIR, "supplementary")
CONFIG_DIR <- file.path(ROOT_DIR, "config")

for (d in c(DERIVED_DIR, OUTPUT_DIR, AUDIT_OUTPUT_DIR, BASELINE_OUTPUT_DIR, REVIEWER_OUTPUT_DIR, MANUSCRIPT_OUTPUT_DIR, SUPPLEMENTARY_OUTPUT_DIR, CONFIG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

source_repo <- function(relative_path, local = parent.frame()) {
  source(file.path(ROOT_DIR, relative_path), local = local, chdir = FALSE)
}
