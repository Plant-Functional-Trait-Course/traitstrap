#' Bootstrap traits
#' @description Function for nonparametric bootstrap resampling to
#' calculate community
#' weighted trait correlations, other bivariate or multivariate statistics
#' @param filled_traits output from the trait_fill function.
#' @param nrep number of bootstrap replicates
#' @param sample_size bootstrap size
#' @param raw logical; argument to extract the raw data of the trait
#' distributions.
#' The default is `raw = FALSE`.
#' If `raw = TRUE`, `nrep` is restricted to 1 to avoid memory issues.
#' @param id column name of unique identifiers of each leaf
#' @param fun bivariate or multivariate function to apply
#' 
#' @details The observed and filled leaves are re-sampled in proportion to
#' their weights, e.g. the abundance of a species or the biomass.
#' Values across all individuals in a community are
#' resampled `sample_size` times to incorporate the full
#' spectrum of trait variation, generating `nrep` trait distributions.
#' The function `fun` is applied to the trait distribution at the finest level
#'  of the filled trait hierarchy.
#'  
#'  Unexpected columns in `filled_traits` are deleted with a warning.
#'
#' Note that due to the flexibility of this function,
#' the output CAN NOT be summarized using
#' `trait_summarise_boot_moments`.
#'
#' @return a tibble with columns for the elements of the scale_hierarchy, 
#' and a list column result which includes the output of `fun`.
#'
#' @importFrom stats var
#' @importFrom e1071 skewness kurtosis
#' @importFrom dplyr slice_sample group_by summarise n_distinct
#' @importFrom tidyr pivot_wider nest unnest any_of
#' @importFrom purrr map list_rbind
#' @importFrom glue glue
#' 
#' @examples
#' require(dplyr)
#' require(tidyr)
#' require(ggplot2)
#' require(purrr)
#'
#' data(community)
#' data(trait)
#'
#' filled_traits <- trait_fill(
#'   comm = community |>
#'     filter(
#'       PlotID %in% c("A", "B"),
#'       Site == 1
#'     ),
#'   traits = trait,
#'   scale_hierarchy = c("Site", "PlotID"),
#'   taxon_col = "Taxon", value_col = "Value",
#'   trait_col = "Trait", abundance_col = "Cover",
#'   complete_only = TRUE, leaf_id = "ID"
#' )
#'
#' # Note that more replicates and a greater sample size are advisable
#' # Here we set them low to make the example run quickly
#' boot_traits <- trait_multivariate_bootstrap(filled_traits,
#'   fun = cor,
#'   nrep = 10,
#'   sample_size = 100
#' )
#'
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

trait_multivariate_bootstrap <- function(filled_traits,
                                         nrep = 100,
                                         sample_size = 200,
                                         raw = FALSE,
                                         id = "ID",
                                         fun) {
  if (raw) {
    nrep <- 1
  }

  attrib <- attr(filled_traits, "attrib")
  value_col <- attrib$value_col

  n_traits <- n_distinct(filled_traits[[attrib$trait_col]])
  trait_names <- unique(filled_traits[[attrib$trait_col]])

  # check complete traits
  check_n_traits <- filled_traits |>
    # remove leaves with incomplete data
    ungroup(.data[[attrib$trait_col]]) |>
    group_by(.data[[id]], .add = TRUE) |>
    mutate(.n = n())

  if (any(check_n_traits$.n != n_traits)) {
    stop("Some leaves with incomplete set of traits.
         Please run trait_fill() with complete_only set to TRUE.")
  }

  # remove unneeded columns
  filled_traits <- filled_traits |>
    select(
      -.data[[attrib$taxon_col]], -.data[[attrib$abundance_col]],
      -.data$n_sample, -.data$max_n_in_sample,
      -.data$level, -.data$sum_abun
    ) 
  
  # check for extra columns and remove with warning if necessary
  expected_columns <- c(as.character(attrib$scale_hierarchy), attrib$taxon_col, 
                        attrib$trait_col, attrib$value_col,
                        attrib$abundance_col, attrib$other_col,
                        attrib$treatment_col, attrib$leaf_id,
                        "weight") 
  if (!all(colnames(filled_traits) %in% expected_columns)) {
    warning(glue("Columns {colnames(filled_traits)[!(colnames(filled_traits) %in% expected_columns)]} were not expected and have been removed."))
    filled_traits <- filled_traits |> 
      select(any_of(expected_columns))
  }
    
    # pivot
    filled_traits_wide <- filled_traits |> 
    pivot_wider(
      names_from = .data[[attrib$trait_col]],
      values_from = .data[[value_col]]
    ) |>
    ungroup(.data[[id]]) |>
    select(-.data[[id]])

  bootstrap_moments <- map(
    seq_len(nrep),
    ~ {
      raw_dist <- slice_sample(filled_traits_wide,
        n = sample_size,
        replace = TRUE, weight_by = weight
      ) |>
        select(-.data$weight) |>
        nest(data = all_of(trait_names))
      if (isTRUE(raw)) {
        return(raw_dist)
      } else {
        # get all the happy moments
        raw_dist |>
          summarise(
            result = map(.data$data, fun)
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
