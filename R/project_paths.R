# Project path helpers ---------------------------------------------------------

find_repo_root <- function() {

  args <- commandArgs(
    trailingOnly = FALSE
  )

  file_arg <- grep(
    "^--file=",
    args,
    value = TRUE
  )

  starts <- c(
    if (length(file_arg)) {
      dirname(
        normalizePath(
          sub(
            "^--file=",
            "",
            file_arg[1]
          ),
          mustWork = FALSE
        )
      )
    } else {
      character()
    },
    getwd()
  )

  starts <- unique(
    normalizePath(
      starts,
      mustWork = FALSE
    )
  )

  for (start in starts) {

    current <- start

    repeat {

      if (
        file.exists(
          file.path(
            current,
            "collembola-reference-frameworks.Rproj"
          )
        )
      ) {
        return(
          normalizePath(
            current,
            mustWork = TRUE
          )
        )
      }

      parent <- dirname(
        current
      )

      if (
        identical(
          parent,
          current
        )
      ) {
        break
      }

      current <- parent
    }
  }

  stop(
    "Could not identify repository root. ",
    "Open collembola-reference-frameworks.Rproj ",
    "or run R/Rscript from inside the repository."
  )
}


ROOT_DIR <- find_repo_root()

DATA_DIR <- file.path(
  ROOT_DIR,
  "data"
)

ANALYSIS_DATA_DIR <- file.path(
  DATA_DIR,
  "analysis"
)

OUTPUT_DIR <- file.path(
  ROOT_DIR,
  "output"
)

ANALYSIS_OUTPUT_DIR <- file.path(
  OUTPUT_DIR,
  "analysis"
)

SENSITIVITY_OUTPUT_DIR <- file.path(
  OUTPUT_DIR,
  "sensitivity"
)

MANUSCRIPT_OUTPUT_DIR <- file.path(
  OUTPUT_DIR,
  "manuscript"
)

SUPPLEMENTARY_OUTPUT_DIR <- file.path(
  OUTPUT_DIR,
  "supplementary"
)


for (d in c(
  ANALYSIS_DATA_DIR,
  OUTPUT_DIR,
  ANALYSIS_OUTPUT_DIR,
  SENSITIVITY_OUTPUT_DIR,
  MANUSCRIPT_OUTPUT_DIR,
  SUPPLEMENTARY_OUTPUT_DIR
)) {

  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


source_repo <- function(
  relative_path,
  local = parent.frame()
) {

  source(
    file.path(
      ROOT_DIR,
      relative_path
    ),
    local = local,
    chdir = FALSE
  )
}
