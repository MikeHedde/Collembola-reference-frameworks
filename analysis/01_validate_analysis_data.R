#!/usr/bin/env Rscript

# ==============================================================================
# 01 - Validate analysis-ready data
# ==============================================================================
# The public workflow starts from site-level indicator values.
# Taxon-by-site abundance matrices and original source workbooks are not
# redistributed with this repository.
# ==============================================================================

source("R/project_paths.R")

message("[01] Validating analysis-ready datasets...")

# ------------------------------------------------------------------------------
# Expected input files
# ------------------------------------------------------------------------------

input_files <- c(
  reference =
    file.path(ANALYSIS_DATA_DIR, "reference_indicator_values.csv"),
  bioindicateur2 =
    file.path(ANALYSIS_DATA_DIR, "bioindicateur2_indicator_values.csv"),
  contrasts =
    file.path(ANALYSIS_DATA_DIR, "bioindicateur2_contrasts.csv"),
  controls =
    file.path(ANALYSIS_DATA_DIR, "bioindicateur2_controls.csv"),
  map =
    file.path(ANALYSIS_DATA_DIR, "map_points.csv")
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing analysis-ready file(s):\n",
    paste(" -", missing_files, collapse = "\n")
  )
}

# ------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------

reference <- read.csv(
  input_files["reference"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

bioindicateur2 <- read.csv(
  input_files["bioindicateur2"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

contrasts <- read.csv(
  input_files["contrasts"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

controls <- read.csv(
  input_files["controls"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

map_points <- read.csv(
  input_files["map"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------------------------
# Column checks
# ------------------------------------------------------------------------------

assert_columns <- function(data, required, data_name) {

  missing <- setdiff(required, names(data))

  if (length(missing) > 0) {
    stop(
      data_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
}

assert_columns(
  reference,
  c(
    "site_id",
    "framework",
    "sampling_strategy",
    "spatial_extent",
    "year",
    "longitude",
    "latitude",
    "CLC_niveau1",
    "RSr",
    "Shannon",
    "Pielou",
    "ab",
    "strict_total",
    "below_rarefaction_n"
  ),
  "reference_indicator_values.csv"
)

assert_columns(
  bioindicateur2,
  c(
    "site_id",
    "longitude",
    "latitude",
    "RSr",
    "Shannon",
    "Pielou",
    "ab"
  ),
  "bioindicateur2_indicator_values.csv"
)

assert_columns(
  contrasts,
  c(
    "Nom_site",
    "treatment_site",
    "treatment",
    "control_site",
    "control_treatment",
    "indicator",
    "raw_treatment",
    "raw_control",
    "LRR"
  ),
  "bioindicateur2_contrasts.csv"
)

assert_columns(
  controls,
  c(
    "Nom_site",
    "control_Site",
    "control_Traitement",
    "control_Contamination",
    "include"
  ),
  "bioindicateur2_controls.csv"
)

assert_columns(
  map_points,
  c(
    "framework",
    "Site",
    "longitude",
    "latitude"
  ),
  "map_points.csv"
)

# ------------------------------------------------------------------------------
# Structural checks
# ------------------------------------------------------------------------------

# site_id is not globally unique:
# some RMQS sites occur in both national and Brittany frameworks,
# and some ANDRA sites were sampled in more than one year.
reference_key <- paste(
  reference$framework,
  reference$site_id,
  ifelse(
    is.na(reference$year),
    "NA",
    reference$year
  ),
  sep = "::"
)

if (anyDuplicated(reference_key)) {
  stop(
    "Duplicate framework + site_id + year combination detected ",
    "in reference_indicator_values.csv"
  )
}

if (anyDuplicated(bioindicateur2$site_id)) {
  stop("Duplicate site_id detected in bioindicateur2_indicator_values.csv")
}

if (nrow(reference) != 756) {
  stop(
    "Expected 756 reference observations; found ",
    nrow(reference)
  )
}

if (nrow(bioindicateur2) != 26) {
  stop(
    "Expected 26 Bioindicateur2 observations; found ",
    nrow(bioindicateur2)
  )
}

if (nrow(contrasts) != 28) {
  stop(
    "Expected 28 treatment-control contrasts; found ",
    nrow(contrasts)
  )
}

if (nrow(controls) != 4) {
  stop(
    "Expected 4 explicit control definitions; found ",
    nrow(controls)
  )
}

# ------------------------------------------------------------------------------
# Framework counts
# ------------------------------------------------------------------------------

expected_total <- c(
  RMQS_Biodiversite = 96,
  RMQS_BioDiv_Bretagne = 98,
  ANDRA = 132,
  TIGA_rural = 430
)

observed_total <- table(reference$framework)

for (fw in names(expected_total)) {

  observed <- if (fw %in% names(observed_total)) {
    as.integer(observed_total[fw])
  } else {
    0L
  }

  if (observed != expected_total[fw]) {
    stop(
      "Unexpected total n for ",
      fw,
      ": expected ",
      expected_total[fw],
      ", found ",
      observed
    )
  }
}

# Agricultural land cover = CLC level 1 class 2
reference_agri <- reference[
  reference$CLC_niveau1 == 2 &
    !is.na(reference$CLC_niveau1),
]

expected_agricultural <- c(
  RMQS_Biodiversite = 58,
  RMQS_BioDiv_Bretagne = 89,
  ANDRA = 90,
  TIGA_rural = 345
)

observed_agricultural <- table(reference_agri$framework)

for (fw in names(expected_agricultural)) {

  observed <- if (fw %in% names(observed_agricultural)) {
    as.integer(observed_agricultural[fw])
  } else {
    0L
  }

  if (observed != expected_agricultural[fw]) {
    stop(
      "Unexpected agricultural n for ",
      fw,
      ": expected ",
      expected_agricultural[fw],
      ", found ",
      observed
    )
  }
}

# ------------------------------------------------------------------------------
# Indicator checks
# ------------------------------------------------------------------------------

for (indicator in c("RSr", "ab")) {

  if (any(!is.finite(reference[[indicator]]))) {
    stop(
      "Non-finite ",
      indicator,
      " values detected in reference data."
    )
  }

  if (any(!is.finite(bioindicateur2[[indicator]]))) {
    stop(
      "Non-finite ",
      indicator,
      " values detected in Bioindicateur2 data."
    )
  }
}

# Contrast indicators used in the manuscript
if (!all(contrasts$indicator %in% c("RSr", "ab"))) {
  stop(
    "Unexpected indicator in bioindicateur2_contrasts.csv"
  )
}

# ------------------------------------------------------------------------------
# Report
# ------------------------------------------------------------------------------

message("")
message("Reference observations: ", nrow(reference))
message("Bioindicateur2 observations: ", nrow(bioindicateur2))
message("Treatment-control contrasts: ", nrow(contrasts))
message("")

message("Reference framework totals:")
print(observed_total)

message("")
message("Agricultural reference populations:")
print(observed_agricultural)

message("")
message("[01] Analysis-ready data validation PASSED.")
