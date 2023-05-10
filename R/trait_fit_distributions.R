#' Fit trait distributions
#' @description Function to fit parametric distributions for each
#' species-by-trait combination at the finest scale of the
#' user-supplied hierarchy.
#' This function returns a tibble containing the fitted parameters.
#' @param filled_traits output from the trait_fill function.
#' @param distribution_type the type of statistical distribution to use.
#' Character or named list. Currently accepts "normal","lognormal", and "beta".
#'
#' @details The distributions can either be a single distribution type which is
#' used for all traits, or traits can be assigned specific distributions types
#' by supplying the function with a named vector of traits, e.g.
#' `c(height = "normal", mass = "lognormal"))`.
#'
#' @return a tibble containing fitted distribution parameters for each trait in each species for each plot.
#'
#' @importFrom stats var
#' @importFrom e1071 skewness kurtosis
#' @importFrom dplyr slice_sample group_by summarise
#' summarize select group_by_at
#' @importFrom stats var

#' @examples
#' library(dplyr)
#' data(community)
#' data(trait)
#'
#' filled_traits <- trait_fill(
#'   comm = community |>
#'     filter(PlotID %in% c("A", "B")),
#'   traits = trait,
#'   scale_hierarchy = c("Site", "PlotID"),
#'   taxon_col = "Taxon", value_col = "Value",
#'   trait_col = "Trait", abundance_col = "Cover"
#' )
#'
#' fitted_distributions <- trait_fit_distributions(
#'   filled_traits = filled_traits,
#'   distribution_type = "normal"
#' )
#' @export
trait_fit_distributions <- function(filled_traits,
                                    distribution_type = "normal") {
  # Check filled traits
  if (!inherits(filled_traits, "filled_trait")) {
    stop("filled traits are not appropriately formatted.
    Please use trait_fill() ")
  }


  # Pull useful information from filled traits object
  value_col <- attributes(filled_traits)$attrib$value_col
  trait_col <- attributes(filled_traits)$attrib$trait_col
  taxon_col <- attributes(filled_traits)$attrib$taxon_col
  abundance_col <- attributes(filled_traits)$attrib$abundance_col
  scale_hierarchy <- attributes(filled_traits)$attrib$scale_hierarchy

  # If only a single type of distribution was supplied, apply it to all traits
  if (length(distribution_type) == 1) {
    distributions <- unique(filled_traits[[trait_col]])
    distribution_type <- rep(x = distribution_type, 
                             times = length(distributions))
    names(distribution_type) <- distributions
  } # if statement

  # Data checks

  # check distribution type is supported
  if (!all(unique(distribution_type) %in% c("beta", "normal", "lognormal"))) {
    x <- distribution_type[!distribution_type %in%
      c("beta", "normal", "lognormal")]
    stop(glue("Distribution type '{x}' not supported "))
  }


  # Beta checks
  if ("beta" %in% distribution_type) {
    # make sure every species x hierarchy combination has at least 2 data points
    beta_counts <- filled_traits |>
      filter(.data[[trait_col]] %in%
        names(distribution_type)[distribution_type == "beta"])

    if (any(beta_counts$n_sample < 2)) {
      stop("Fitting a beta distrbution requires 2+ points per distribution.
         We suggest re-imputing traits with
         a minimum sample size of (at least) 2 traits.")
    }

    # check that values are between 0 and 1
    beta_vals <- filled_traits |>
      ungroup() |>
      filter(.data[[trait_col]] %in%
        names(distribution_type)[distribution_type == "beta"]) |>
      select(all_of(value_col))

    if (any(beta_vals > 1 | beta_vals < 0)) {
      stop("For a Beta distribution values must be between 0 and 1.")
    }
  } # end beta checks


  # lognormal checks
  if ("lognormal" %in% distribution_type) {
    ln_vals <- filled_traits |>
      ungroup() |>
      filter(.data[[trait_col]] %in%
        names(distribution_type)[distribution_type == "lognormal"]) |>
      select(all_of(value_col))

    if (any(ln_vals <= 0)) {
      stop("For a lognormal distribution values must be positive.")
    }
  } # end lognormal checks

  # Main body

  distribution_parms <- filled_traits |>
    group_by(.data[[c(taxon_col)]],
      .data[[c(abundance_col)]], .data$n_sample,
      .add = TRUE
    ) |>
    mutate(distribution_type =
      distribution_type[.data[[trait_col]]]
    ) |>
    summarize(
      quiet(suppressWarnings(get_dist_parms(
        data = .data[[value_col]],
        distribution_type = .data$distribution_type[1]
      ))),
      .groups = "keep"
    )

  # set arguments as attributes so next functions have access to them
  attr(distribution_parms, "attrib") <- attributes(filled_traits)$attrib

  class(distribution_parms) <- c(
    "parametric_distributions",
    class(distribution_parms)
  )

  return(distribution_parms)
} # fx
