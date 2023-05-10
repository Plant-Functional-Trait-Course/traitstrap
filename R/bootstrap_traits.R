#' Bootstrap traits
#' @description Function for nonparametric bootstrap resampling to
#' calculate community weighted trait mean and higher moments.
#' @param filled_traits output from the trait_fill function.
#' @param nrep number of bootstrap replicates
#' @param sample_size bootstrap size
#' @param raw logical; argument to extract the raw data of the trait
#' distributions.
#' The default is `raw = FALSE`. If `raw = TRUE`, `nrep` is restricted to 1 to
#' avoid memory issues.
#'
#' @details The observed traits are re-sampled in proportion to their weights,
#' e.g. the abundance of a species or the biomass.
#' Values across all individuals
#' in a community are resampled `sample_size` times to incorporate the full
#' spectrum of trait variation, generating `nrep` trait distributions.
#' From these distributions the function estimates the mean and the higher
#' moments including variance, skewness and kurtosis.
#'
#' #' The output of `trait_np_bootstrap()` can be summarized using
#' `trait_summarize_boot_moments()`.
#'
#' @return a tibble with columns for each grouping variable of `filled_traits`
#'  (usually the elements of scale_hierarchy and the traits column), 
#'  and the moments mean, variance, skewness, and kurtosis.
#'
#' @importFrom stats var
#' @importFrom e1071 skewness kurtosis
#' @importFrom dplyr slice_sample group_by summarise
#' @importFrom purrr map list_rbind
#' @examples
#' library(dplyr)
#' data(community)
#' data(trait)
#'
#' # Filter community data to make example faster
#' community <- community |>
#'   filter(
#'     PlotID %in% c("A", "B"),
#'     Site == 1
#'   )
#' filled_traits <- trait_fill(
#'   comm = community,
#'   traits = trait,
#'   scale_hierarchy = c("Site", "PlotID"),
#'   taxon_col = "Taxon", value_col = "Value",
#'   trait_col = "Trait", abundance_col = "Cover"
#' )
#'
#' boot_traits <- trait_np_bootstrap(filled_traits,
#'   nrep = 20,
#'   sample_size = 200
#' )
#' @export

trait_np_bootstrap <- function(filled_traits,
                               nrep = 100,
                               sample_size = 200,
                               raw = FALSE) {
  if (isTRUE(raw)) {
    nrep <- 1
  }

  attrib <- attr(filled_traits, "attrib")
  value_col <- attrib$value_col
  bootstrap_moments <- map(
    seq_len(nrep),
    ~ {
      raw_dist <- slice_sample(filled_traits,
        n = sample_size,
        replace = TRUE, weight_by = weight
      )
      if (raw) {
        return(raw_dist)
      } else {
        # get all the happy moments
        raw_dist |>
          summarise(
            mean = mean(.data[[value_col]]),
            variance = var(.data[[value_col]]),
            skewness = skewness(.data[[value_col]]),
            kurtosis = kurtosis(.data[[value_col]])
          )
      }
    },
    .id = "n"
  ) |>
    list_rbind()

  attr(bootstrap_moments, "attrib") <- attrib

  # make bootstrap_moments an ordinary tibble
  class(bootstrap_moments) <-
    class(bootstrap_moments)[!class(bootstrap_moments) == "filled_trait"]
  return(bootstrap_moments)
}
