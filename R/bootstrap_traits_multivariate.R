#' Bootstrap traits
#' @description Function for nonparametric bootstrap resampling to calculate community 
#' weighted trait correlations, other bivariate or multivariate statistics
#' @param imputed_traits output from the trait_impute function.
#' @param nrep number of bootstrap replicates
#' @param sample_size bootstrap size 
#' @param raw logical; argument to extract the raw data of the trait distributions.
#' The default is raw = FALSE. If raw = TRUE, nrep is restricted to 1 to avoid 
#' memory issues.
#' @param id column name of unique identifiers of each leaf
#' @param fun bivariate or multivate function to apply
#' @description
#' 
#' @details The observed traits are re-sampled in proportion to their weights, 
#' e.g. the abundance of a species or the biomass. Values across all individuals 
#' in a community are resampled n times (nrep) to incorporate the full 
#' spectrum of trait variation, generating n number (sample_size) of trait distributions. 
#' From these distributions the function estimates the mean and the higher moments
#' including variance, skewness and kurtosis.
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
#' @importFrom tidyr pivot_wider nest
#' @importFrom purrr map_dfr
#' @examples 
#' data(community)
#' data(trait)
#' imputed_traits <-trait_impute(comm = community, traits = trait,
#'                  scale_hierarchy = c("Site", "PlotID"),
#'                  taxon_col = "Taxon", value_col = "Value",
#'                  trait_col = "Trait", abundance_col = "Cover")
#' boot_traits <- trait_multivariate_bootstrap(imputed_traits, fun = cor)
#' @export

trait_multivariate_bootstrap <- function(imputed_traits, nrep = 100, sample_size = 200, raw = FALSE, id = "ID", fun) {
  if (raw) {nrep <- 1}
  #  stopifnot(class(traits_com) == "imputed_traits")
  attrib <- attr(imputed_traits, "attrib")
  value_col <- attrib$value_col
  
  n_traits <- n_distinct(imputed_traits[[attrib$trait_col]])
  trait_names <- unique(imputed_traits[[attrib$trait_col]])
  
  # pivot wider
  imputed_traits_wide <- imputed_traits %>%
    # remove leaves with incomplete data
    ungroup(.data[[attrib$trait_col]]) %>% 
    group_by(.data[[id]], .add = TRUE) %>% 
    filter(n() == n_traits) %>% 
    #remove unneeded columns
    select(-.data[[attrib$taxon_col]], -.data[[attrib$abundance_col]], 
           -.data$n_sample, -.data$max_n_in_sample, 
           -.data$level, -.data$sum_abun) %>%
    #pivot
    pivot_wider(names_from = .data[[attrib$trait_col]],
                values_from = .data[[value_col]]) %>% 
    ungroup(.data[[id]]) %>% 
    select(-.data[[id]])
  
  bootstrap_moments <- map_dfr(
    1:nrep,
    ~ {
      raw_dist <- slice_sample(imputed_traits_wide, n = sample_size,
                               replace = TRUE, weight_by = weight) %>% 
        select(-.data$weight) %>% 
        nest(data = all_of(trait_names))
      if (isTRUE(raw)){
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
    class(bootstrap_moments)[!class(bootstrap_moments) == "imputed_trait"]
  return(bootstrap_moments)
}
