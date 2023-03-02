#' Bootstrap traits parametrically
#' @description Function for parametric bootstrap resampling to calculate
#' community weighted trait mean and higher moments.
#' @param fitted_distributions
#' Fitted distribution object returned by trait_fit_distributions
#' @param nrep number of bootstrap replicates
#' @param sample_size bootstrap size
#' @param raw logical; argument to extract the raw data of the
#' trait distributions.
#' The default is `raw = FALSE`.
#' If `raw = TRUE`, `nrep` is restricted to 1 to avoid
#' memory issues.
#'
#' @details `trait_parametric_bootstrap()` is a parametric analogue of the
#' `trait_np_bootstrap()`.
#' It randomly samples from among the fitted distributions
#' proportionally to species abundance.
#' The number of samples per replicated are drawn
#' specified with the parameter sample_size,
#' and the number of replicates is specified
#' by the parameter` nrep`.
#' From these distributions the function estimates the mean
#' and the higher moments including variance, skewness and kurtosis.
#'
#' The output of `trait_parametric_bootstrap()` can be summarized using
#' `trait_summarize_boot_moments()`.

#' @return a tibble
#'
#' @importFrom stats var
#' @importFrom e1071 skewness kurtosis
#' @importFrom dplyr slice_sample group_by summarise
#' @importFrom tidyr unnest
#' @importFrom purrr map list_rbind
#' @examples
#' library(dplyr)
#' data(community)
#' data(trait)
#'
#' # Filter trait and community data to make example faster
#'
#' community <- community |>
#'   filter(
#'     PlotID %in% c("A", "B"),
#'     Site == 1
#'   )
#'
#' trait <- trait |>
#'   filter(Trait %in% c("Plant_Height_cm"))
#'
#' filled_traits <- trait_fill(
#'   comm = community,
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
#'
#' # Note that more replicates and a greater sample size are advisable
#' # Here we set them low to make the example run quickly
#' parametric_distributions <- trait_parametric_bootstrap(
#'   fitted_distributions = fitted_distributions,
#'   nrep = 5,
#'   sample_size = 100
#' )
#'
#' moment_summary <- trait_summarise_boot_moments(
#'   bootstrap_moments = parametric_distributions,
#'   parametric = FALSE
#' )
#' @export
trait_parametric_bootstrap <- function(fitted_distributions,
                                       nrep = 100,
                                       sample_size = 200,
                                       raw = FALSE) {
  if (raw) {
    nrep <- 1
  }
  # Check that inputs makes sense

  if (!inherits(fitted_distributions, "parametric_distributions")) {
    stop("Fitted distributions not properly formatted.
         Please use trait_fit_distributions()")
  }

  if (!all(nrep %% 1 == 0 & sample_size %% 1 == 0)) {
    stop("nrep and sample_size should be integers")
  }


  # Pull useful information from filled traits object
  trait_col <- attributes(fitted_distributions)$attrib$trait_col
  abundance_col <- attributes(fitted_distributions)$attrib$abundance_col
  taxon_col <- attributes(fitted_distributions)$attrib$taxon_col
  scale_hierarchy <- attributes(fitted_distributions)$attrib$scale_hierarchy
  attrib <- attr(fitted_distributions, "attrib")

  bootstrap_moments <- map(
    seq_len(nrep),
    ~ {
      raw_dist <- fitted_distributions |>
        group_by_at(c(as.character(scale_hierarchy), trait_col)) |>
        slice_sample(
          n = sample_size,
          replace = TRUE,
          weight_by = .data[[abundance_col]]
        ) |>
        group_by_at(c(
          as.character(scale_hierarchy),
          trait_col, taxon_col,
          "parm1", "parm2", "distribution_type"
        )) |>
        summarise(n_drawn = n(), .groups = "keep") |>
        mutate(
          draw_value =
            list(distribution_handler(
              parm1 = parm1,
              parm2 = parm2,
              n = n_drawn,
              type = distribution_type
            ))
        ) |>
        group_by_at(c(as.character(scale_hierarchy), trait_col)) |>
        unnest(draw_value)

      if (raw) {
        return(raw_dist)
      } else {
        raw_dist |>
          summarise(
            mean = mean(draw_value),
            variance = var(draw_value),
            skewness = skewness(draw_value),
            kurtosis = kurtosis(draw_value),
            .groups = "keep"
          )
      }
    },
    .id = "n"
  ) |>
    list_rbind()

  attr(bootstrap_moments, "attrib") <- attrib

  # make bootstrap_moments an ordinary tibble
  class(bootstrap_moments) <-
    class(bootstrap_moments)[!class(bootstrap_moments) ==
      "parametric_distributions"]
  return(bootstrap_moments)
} # end fx
