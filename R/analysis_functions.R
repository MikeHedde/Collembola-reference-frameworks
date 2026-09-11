# Sensitivity-analysis helper functions ----------------------------------------

required_packages(c(
  "dplyr",
  "tidyr",
  "tibble"
))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})


# Extreme-value sensitivity ----------------------------------------------------

winsorise <- function(
  x,
  probs = c(0.01, 0.99)
) {

  q <- stats::quantile(
    x,
    probs = probs,
    na.rm = TRUE,
    type = 7
  )

  pmax(
    pmin(x, q[[2]]),
    q[[1]]
  )
}


extreme_iqr_flag <- function(
  x,
  k = 3
) {

  q <- stats::quantile(
    x,
    probs = c(0.25, 0.75),
    na.rm = TRUE,
    type = 7
  )

  iqr <- q[[2]] - q[[1]]

  x < (q[[1]] - k * iqr) |
    x > (q[[2]] + k * iqr)
}


# Equal-size resampling --------------------------------------------------------

balanced_subsampling_breaks <- function(
  ref_data,
  n,
  B,
  seed,
  frameworks = FRAMEWORKS,
  indicators = MANUSCRIPT_INDICATORS
) {

  set.seed(seed)

  n_by_framework <- table(
    as.character(
      ref_data$framework
    )
  )

  if (
    any(
      n_by_framework[frameworks] < n
    )
  ) {
    stop(
      "At least one framework has fewer than n = ",
      n,
      " eligible observations."
    )
  }

  result <- vector(
    "list",
    B *
      length(frameworks) *
      length(indicators)
  )

  k <- 1L

  for (b in seq_len(B)) {

    sampled_rows <- lapply(
      frameworks,
      function(fw) {

        idx <- which(
          as.character(
            ref_data$framework
          ) == fw
        )

        sample(
          idx,
          size = n,
          replace = FALSE
        )
      }
    )

    names(sampled_rows) <-
      frameworks

    for (fw in frameworks) {

      d <- ref_data[
        sampled_rows[[fw]],
        ,
        drop = FALSE
      ]

      for (ind in indicators) {

        tmp <- breaks_to_table(
          score_breaks(
            d[[ind]]
          ),
          fw,
          ind,
          "balanced_n58"
        )

        tmp$replicate <- b

        result[[k]] <- tmp
        k <- k + 1L
      }
    }
  }

  bind_rows(result)
}


# Cross-framework concordance --------------------------------------------------

pairwise_score_agreement <- function(
  score_long
) {

  frameworks <- unique(
    score_long$framework
  )

  pairs <- combn(
    frameworks,
    2,
    simplify = FALSE
  )

  out <- list()
  k <- 1L

  for (ind in unique(
    score_long$indicator
  )) {

    dat_ind <- score_long %>%
      filter(
        .data$indicator == ind
      )

    for (pair in pairs) {

      wide <- dat_ind %>%
        filter(
          .data$framework %in% pair
        ) %>%
        select(
          site_id,
          framework,
          score
        ) %>%
        distinct() %>%
        pivot_wider(
          names_from = "framework",
          values_from = "score"
        )

      a <- wide[[pair[1]]]
      b <- wide[[pair[2]]]

      ok <- complete.cases(a, b)

      a <- a[ok]
      b <- b[ok]

      broad_a <- cut(
        a,
        breaks = c(
          -Inf,
          2,
          4,
          Inf
        ),
        labels = c(
          "0-2",
          "3-4",
          "5-6"
        )
      )

      broad_b <- cut(
        b,
        breaks = c(
          -Inf,
          2,
          4,
          Inf
        ),
        labels = c(
          "0-2",
          "3-4",
          "5-6"
        )
      )

      out[[k]] <- tibble(
        indicator = ind,
        framework_1 = pair[1],
        framework_2 = pair[2],
        n = length(a),
        spearman_rho =
          if (length(a) >= 3L) {
            suppressWarnings(
              stats::cor(
                a,
                b,
                method = "spearman"
              )
            )
          } else {
            NA_real_
          },
        exact_agreement =
          if (length(a)) {
            mean(a == b)
          } else {
            NA_real_
          },
        within_one_class =
          if (length(a)) {
            mean(
              abs(a - b) <= 1
            )
          } else {
            NA_real_
          },
        broad_class_agreement =
          if (length(a)) {
            mean(
              broad_a == broad_b
            )
          } else {
            NA_real_
          }
      )

      k <- k + 1L
    }
  }

  bind_rows(out)
}


# Equal-weight pooled reference ------------------------------------------------

build_balanced_pooled_breaks <- function(
  ref_data,
  n,
  B,
  seed,
  indicators = MANUSCRIPT_INDICATORS
) {

  set.seed(seed)

  boot <- vector(
    "list",
    B * length(indicators)
  )

  k <- 1L

  for (b in seq_len(B)) {

    idx <- unlist(
      lapply(
        FRAMEWORKS,
        function(fw) {

          candidate <- which(
            as.character(
              ref_data$framework
            ) == fw
          )

          sample(
            candidate,
            n,
            replace = FALSE
          )
        }
      )
    )

    d <- ref_data[
      idx,
      ,
      drop = FALSE
    ]

    for (ind in indicators) {

      tmp <- breaks_to_table(
        score_breaks(
          d[[ind]]
        ),
        "POOLED_BALANCED",
        ind,
        "pooled_balanced"
      )

      tmp$replicate <- b

      boot[[k]] <- tmp
      k <- k + 1L
    }
  }

  boot <- bind_rows(boot)

  summary <- boot %>%
    group_by(
      .data$framework,
      .data$indicator,
      .data$boundary,
      .data$probability,
      .data$scenario
    ) %>%
    summarise(
      value =
        median(
          .data$value,
          na.rm = TRUE
        ),
      q025 =
        quantile(
          .data$value,
          0.025,
          na.rm = TRUE
        ),
      q975 =
        quantile(
          .data$value,
          0.975,
          na.rm = TRUE
        ),
      .groups = "drop"
    )

  list(
    draws = boot,
    summary = summary
  )
}
