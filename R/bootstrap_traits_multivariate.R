#' Bootstrap traits
#' @description Function for nonparametric bootstrap resampling to
#' calculate community
#' weighted trait correlations, other bivariate or multivariate statistics
#' @param selected_traits output from the trait_fill function.
#' @param nrep number of bootstrap replicates
#' @param sample_size bootstrap size
#' @param raw logical; argument to extract the raw data of the trait
#' distributions.
#' The default is `raw = FALSE`.
#' If `raw = TRUE`, `nrep` is restricted to 1 to avoid memory issues.
#' @param id column name of unique identifiers of each leaf
#' @param fun bivariate or multivariate function to apply
#'
#' @details The observed and selected leaves are re-sampled in proportion to
#' their weights, e.g. the abundance of a species or the biomass.
#' Values across all individuals
#' in a community are resampled `sample_size` times to incorporate the full
#' spectrum of trait variation, generating `nrep` trait distributions.
#' From these distributions the function estimates the mean and the higher
#' moments including variance, skewness and kurtosis.
#'
#' #' The output of trait_np_bootstrap can be summarized using
#' trait_summarize_boot_moments.
#'
#' @return a tibble
#'
#' @importFrom stats var
#' @importFrom e1071 skewness kurtosis
#' @importFrom magrittr %>%
#' @importFrom dplyr slice_sample group_by summarise n_distinct
#' @importFrom tidyr pivot_wider nest unnest
#' @importFrom purrr map_dfr
#' @examples
#' require(dplyr)
#' require(tidyr)
#' require(ggplot2)
#' require(purrr)
#' data(community)
#' data(trait)
#' selected_traits <- trait_fill(
#'   comm = community, traits = trait,
#'   scale_hierarchy = c("Site", "PlotID"),
#'   taxon_col = "Taxon", value_col = "Value",
#'   trait_col = "Trait", abundance_col = "Cover",
#'   complete_only = TRUE, leaf_id = "ID"
#' )
#' boot_traits <- trait_multivariate_bootstrap(selected_traits, fun = cor)
#' boot_traits_long <- boot_traits |>
#'   mutate(correlations = map(result, ~ cor_to_df(.x))) |>
#'   select(-result) |>
#'   unnest(correlations)
#'
#' boot_traits_long |>
#'   ggplot(aes(x = paste(row, "v", col), y = value)) +
#'   geom_violin() +
#'   facet_grid(Site ~ PlotID) +
#'   coord_flip() +
#'   labs(y = "Correlation", x = "")

#' @export

trait_multivariate_bootstrap <- function(selected_traits,
                                         nrep = 100,
                                         sample_size = 200,
                                         raw = FALSE,
                                         id = "ID",
                                         fun) {
  if (raw) {
    nrep <- 1
  }

  attrib <- attr(selected_traits, "attrib")
  value_col <- attrib$value_col

  n_traits <- n_distinct(selected_traits[[attrib$trait_col]])
  trait_names <- unique(selected_traits[[attrib$trait_col]])

  # check complete traits
  check_n_traits <- selected_traits %>%
    # remove leaves with incomplete data
    ungroup(.data[[attrib$trait_col]]) %>%
    group_by(.data[[id]], .add = TRUE) |>
    mutate(.n = n())

  if (any(check_n_traits$.n != n_traits)) {
    stop("Some leaves with incomplete set of traits.
         Please run trait_fill() with complete_only set to TRUE.")
  }

  # pivot_wider
  selected_traits_wide <- selected_traits %>%
    # remove unneeded columns
    select(
      -.data[[attrib$taxon_col]], -.data[[attrib$abundance_col]],
      -.data$n_sample, -.data$max_n_in_sample,
      -.data$level, -.data$sum_abun
    ) %>%
    # pivot
    pivot_wider(
      names_from = .data[[attrib$trait_col]],
      values_from = .data[[value_col]]
    ) %>%
    ungroup(.data[[id]]) %>%
    select(-.data[[id]])

  bootstrap_moments <- map_dfr(
    1:nrep,
    ~ {
      raw_dist <- slice_sample(selected_traits_wide,
        n = sample_size,
        replace = TRUE, weight_by = weight
      ) %>%
        select(-.data$weight) %>%
        nest(data = all_of(trait_names))
      if (isTRUE(raw)) {
        return(raw_dist)
      } else {
        # get all the happy moments
        raw_dist %>%
          summarise(
            result = map(.data$data, fun)
          )
      }
    },
    .id = "n"
  )

  attr(bootstrap_moments, "attrib") <- attrib

  # make bootstrap_moments an ordinary tibble
  class(bootstrap_moments) <-
    class(bootstrap_moments)[!class(bootstrap_moments) == "selected_trait"]
  return(bootstrap_moments)
}
