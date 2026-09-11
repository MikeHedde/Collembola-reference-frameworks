# Indicator scoring functions --------------------------------------------------
#
# The public workflow starts from site-level indicator values already computed
# from the original taxonomic data. This file contains only the functions needed
# to derive empirical quantile boundaries and ordinal scores.

required_packages(c("dplyr", "tibble"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

FRAMEWORKS <- c(
  "RMQS_Biodiversite",
  "RMQS_BioDiv_Bretagne",
  "ANDRA",
  "TIGA_rural"
)

FRAMEWORK_LABELS <- c(
  RMQS_Biodiversite = "RMQS-Biodiversité",
  RMQS_BioDiv_Bretagne = "RMQS BioDiv Bretagne",
  ANDRA = "ANDRA",
  TIGA_rural = "TIGA rural"
)

# RSr and density are the two indicators analysed in the main manuscript.
# Shannon and Pielou are retained in the analysis-ready dataset for
# supplementary descriptive analyses.
SUPPLEMENTARY_INDICATORS <- c(
  "RSr",
  "Shannon",
  "Pielou",
  "ab"
)

MANUSCRIPT_INDICATORS <- c(
  "RSr",
  "ab"
)

N_SCORE_CLASSES <- 7L
RAREFACTION_N <- 500L
AGRICULTURAL_CLC1 <- 2L


# Empirical score boundaries ---------------------------------------------------

score_breaks <- function(
  x,
  nb_class = N_SCORE_CLASSES
) {

  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    stop(
      "At least two finite reference values are required ",
      "to build scores."
    )
  }

  if (nb_class < 3L) {
    stop("nb_class must be >= 3.")
  }

  n_internal_intervals <- nb_class - 2L

  empirical_boundaries <- as.numeric(
    stats::quantile(
      x,
      probs = seq(
        0,
        1,
        length.out = n_internal_intervals + 1L
      ),
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )

  c(
    -Inf,
    empirical_boundaries,
    Inf
  )
}


score_values <- function(
  x,
  breaks
) {

  out <- rep(
    NA_integer_,
    length(x)
  )

  ok <- is.finite(x)

  if (!any(ok)) {
    return(out)
  }

  out[ok] <- findInterval(
    x[ok],
    breaks,
    rightmost.closed = FALSE,
    all.inside = TRUE
  ) - 1L

  out
}


breaks_to_table <- function(
  breaks,
  framework,
  indicator,
  scenario
) {

  internal <- breaks[
    is.finite(breaks)
  ]

  probabilities <- seq(
    0,
    1,
    length.out = length(internal)
  )

  tibble(
    scenario = scenario,
    framework = framework,
    indicator = indicator,
    boundary = paste0(
      "q",
      sprintf(
        "%02d",
        round(100 * probabilities)
      )
    ),
    probability = probabilities,
    value = internal
  )
}


build_break_table <- function(
  data,
  frameworks = FRAMEWORKS,
  indicators = MANUSCRIPT_INDICATORS,
  scenario = "reference_framework"
) {

  out <- list()
  k <- 1L

  for (fw in frameworks) {

    ref_fw <- data[
      as.character(data$framework) == fw,
      ,
      drop = FALSE
    ]

    if (!nrow(ref_fw)) {
      stop(
        "No observations found for framework: ",
        fw
      )
    }

    for (ind in indicators) {

      out[[k]] <- breaks_to_table(
        score_breaks(
          ref_fw[[ind]]
        ),
        framework = fw,
        indicator = ind,
        scenario = scenario
      )

      k <- k + 1L
    }
  }

  bind_rows(out)
}


get_breaks_from_table <- function(
  table,
  framework,
  indicator
) {

  values <- table %>%
    filter(
      .data$framework == .env$framework,
      .data$indicator == .env$indicator
    ) %>%
    arrange(
      .data$probability
    ) %>%
    pull(
      "value"
    )

  if (!length(values)) {
    stop(
      "No boundaries found for ",
      framework,
      " / ",
      indicator
    )
  }

  if (
    anyNA(values) ||
    any(!is.finite(values))
  ) {
    stop(
      "Non-finite empirical boundaries found for ",
      framework,
      " / ",
      indicator
    )
  }

  if (
    is.unsorted(
      values,
      strictly = FALSE
    )
  ) {
    stop(
      "Empirical boundaries are not ordered for ",
      framework,
      " / ",
      indicator
    )
  }

  c(
    -Inf,
    values,
    Inf
  )
}


score_validation_data <- function(
  validation_data,
  break_table,
  frameworks = unique(
    break_table$framework
  ),
  indicators = unique(
    break_table$indicator
  ),
  scenario_name = "reference_framework"
) {

  if (
    !"site_id" %in%
      names(validation_data)
  ) {
    stop(
      "validation_data must contain a site_id column."
    )
  }

  out <- list()
  k <- 1L

  for (fw in frameworks) {

    for (ind in indicators) {

      boundaries <-
        get_breaks_from_table(
          break_table,
          fw,
          ind
        )

      out[[k]] <- tibble(
        site_id =
          as.character(
            validation_data$site_id
          ),
        indicator = ind,
        framework = fw,
        scenario = scenario_name,
        raw_value =
          validation_data[[ind]],
        score =
          score_values(
            validation_data[[ind]],
            boundaries
          )
      )

      k <- k + 1L
    }
  }

  bind_rows(out)
}
