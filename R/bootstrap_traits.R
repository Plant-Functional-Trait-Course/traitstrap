#' Bootstrap traits
#' @description Bootstrap impute traits 
#' @param imputed_traits imputed trait and community data in long format
#' @param nrep number of bootstrap replicates
#' @param sample_size number of plants per sample
#' @description 
#' 
#' @return a tibble
#' 
#' @importFrom stats var 
#' @importFrom e1071 skewness kurtosis
#' @importFrom magrittr %>%
#' @importFrom dplyr slice_sample group_by summarise across any_of
#' @importFrom purrr map_df
#' @export

trait_np_bootstrap <- function(imputed_traits, nrep = 100, sample_size = 200){  
  #  stopifnot(class(traits_com) == "imputed_traits")
  attrib <- attr(imputed_traits, "attrib")
  value_col <- attrib$value_col
  bootstrapMoments <- map_df(
    1:nrep,
    ~{slice_sample(imputed_traits, n = sample_size,  
                   replace = TRUE, weight_by = weight) %>% 
        # get all the happy moments
        summarise(
          across(any_of(value_col)), 
          list(mean = mean, variance = var, 
               skewness = skewness, kurtosis = kurtosis))},
    .id = "n"
  )
  
  attr(bootstrapMoments, "attrib") <- attrib
  class(bootstrapMoments)<- class(bootstrapMoments)[!class(bootstrapMoments) == "imputed_trait"]#make bootstrapMoments an ordinary tibble
  return(bootstrapMoments)
}

