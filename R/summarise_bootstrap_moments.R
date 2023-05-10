#' Summarise Bootstrap traits
#' @description Summarizes the mean and confidence interval for each
#' trait moment.
#' @param bootstrap_moments trait moments from trait_np_bootstrap or
#' trait_parametric_bootstrap
#' @param parametric logical; default is TRUE.
#' Should Confidence Intervals be calculated parametrically
#' (using the mean and SD) or nonparametrically (using quantiles).
#' @param sd_mult Number of standard deviations around each moment,
#' defaults to one
#' @param ci Desired confidence level for use when parametric is false.
#' Defaults to 0.95.

#'
#' @return tibble with the grouping variables and the mean of each moment
#'  (+/- sd_mult * SD)
#'
#' @importFrom dplyr n group_by any_of summarise across rename
#' @importFrom rlang .data

#' @examples
#' library(dplyr)
#' data(community)
#' data(trait)
#'
#' # Filter community data to make example faster
#' community <- community |>
#'   filter(PlotID %in% c("A", "B"))
#'
#' filled_traits <- trait_fill(
#'   comm = community,
#'   traits = trait,
#'   scale_hierarchy = c("Site", "PlotID"),
#'   taxon_col = "Taxon", value_col = "Value",
#'   trait_col = "Trait", abundance_col = "Cover"
#' )
#'
#' # Note that more replicates and a greater sample size are advisable
#' # Here we set them low to make the example run quickly
#' boot_traits <- trait_np_bootstrap(filled_traits,
#'   nrep = 20,
#'   sample_size = 100
#' )
#'
#' trait_summarise_boot_moments(boot_traits)
#' @export

trait_summarise_boot_moments <- function(bootstrap_moments,
                                         parametric = TRUE,
                                         sd_mult = 1, ci = 0.95) {
  attrib <- attr(bootstrap_moments, "attrib")
  groups <- c(
    as.character(attrib$scale_hierarchy),
    attrib$trait_col,
    attrib$other_col
  )

  # add treatment_col to groups if used
  if (!is.null(attrib$treatment_col)) {
    groups <- c(groups, paste0(attrib$treatment_col, "_comm"))
  }


  # calculate means of moments
  summ_bootstrap_moments <- bootstrap_moments |>
    ungroup() |>
    group_by(across(any_of(groups))) |>
    rename(MEAN = "mean") |> # capitalise to simplify code below
    summarise(
      n = n(),
      mean = mean(.data$MEAN),
      ci_low_mean = get_ci(
        data = .data$MEAN, sd_mult = sd_mult, ci = ci,
        which = "low", parametric = parametric
      ),
      ci_high_mean = get_ci(
        data = .data$MEAN, sd_mult = sd_mult, ci = ci,
        which = "high", parametric = parametric
      ),
      var = mean(.data$variance),
      ci_low_var = get_ci(
        data = .data$variance, sd_mult = sd_mult, ci = ci,
        which = "low", parametric = parametric
      ),
      ci_high_var = get_ci(
        data = .data$variance, sd_mult = sd_mult, ci = ci,
        which = "high", parametric = parametric
      ),
      skew = mean(.data$skewness),
      ci_low_skew = get_ci(
        data = .data$skewness, sd_mult = sd_mult, ci = ci,
        which = "low", parametric = parametric
      ),
      ci_high_skew = get_ci(
        data = .data$skewness, sd_mult = sd_mult, ci = ci,
        which = "high", parametric = parametric
      ),
      kurt = mean(.data$kurtosis),
      ci_low_kurt = get_ci(
        data = .data$kurtosis, sd_mult = sd_mult, ci = ci,
        which = "low", parametric = parametric
      ),
      ci_high_kurt = get_ci(
        data = .data$kurtosis, sd_mult = sd_mult, ci = ci,
        which = "high", parametric = parametric
      )
    )



  return(summ_bootstrap_moments)
}
