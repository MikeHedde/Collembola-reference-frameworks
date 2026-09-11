# Input/output helpers ---------------------------------------------------------

required_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required R package(s): ", paste(missing, collapse = ", "),
      "\nInstall them before running the workflow."
    )
  }
  invisible(TRUE)
}

assert_files_exist <- function(paths, context = "Required input") {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      context, " file(s) not found:\n",
      paste0(" - ", missing, collapse = "\n"),
      "\nSee data/README.md."
    )
  }
  invisible(paths)
}

write_csv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file = path, row.names = FALSE, na = "")
  invisible(path)
}

write_text_lines <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(x, con = path, useBytes = TRUE)
  invisible(path)
}

as_clean_character <- function(x) {
  out <- as.character(x)
  out[is.na(x)] <- NA_character_
  trimws(out)
}

normalise_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))
