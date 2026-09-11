#!/usr/bin/env Rscript

# ==============================================================================
# Run the complete reproducible workflow
# ==============================================================================
#
# Each analysis step is executed in a fresh R session. This prevents objects
# created by one script from becoming hidden dependencies of later scripts.
#
# ==============================================================================

source("R/project_paths.R")

steps <- c(
  "analysis/01_validate_analysis_data.R",
  "analysis/02_main_analysis.R",
  "analysis/03_sensitivity_analyses.R",
  "analysis/04_manuscript_figures.R",
  "analysis/05_supplementary_material.R"
)

rscript <- file.path(
  R.home("bin"),
  "Rscript"
)

for (step in steps) {

  message(
    "\n============================================================"
  )

  message(
    "Running ",
    step
  )

  message(
    "============================================================"
  )

  status <- system2(
    rscript,
    args = file.path(
      ROOT_DIR,
      step
    )
  )

  if (status != 0L) {
    stop(
      "Workflow stopped because ",
      step,
      " returned exit status ",
      status,
      "."
    )
  }
}

message(
  "\n============================================================"
)

message(
  "Complete workflow finished successfully."
)

message(
  "============================================================"
)
