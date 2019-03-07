#' Impute traits
#' @description Impute traits 
#' @param comm community data in long format
#' @param traits trait data in long format
#' @param scale_hierarchy character vector of site/block/plot hierarchy (large to small)
#' @param taxon_col character; name of taxon column in comm and traits
#' @param trait_col character; name of trait name column in traits
#' @param value_col character; name of trait value column in traits
#' @param abundance_col character; name of species abundance column in comm
#' @param other_col name of other grouping columns in comm
#' 
#' @description 
#' 
#' @return 
#' 
#' @importFrom stats sd var weighted.mean
#' @importFrom magrittr %>%
#' @importFrom dplyr select one_of mutate group_by filter left_join n
#' @importFrom purrr map_df
#' @importFrom rlang !!! maybe_missing .data
 
trait_impute <- function(comm, traits, scale_hierarchy = c("Country", "Site", "BlockID", "PlotID"), taxon_col = "taxon",  trait_col = "trait", value_col = "Value", abundance_col = "Cover", other_col = ""){
  out <- length(scale_hierarchy):1 %>%  
    map_df(~{#browser()
      scale_level <- .
      scale_drop <- character(0L)
      if(scale_level < length(scale_hierarchy)){
        scale_drop <- scale_hierarchy[(scale_level + 1):length(scale_hierarchy)]
      }
      scale_keep <- scale_hierarchy[1:scale_level]

     traits <- traits %>% select(-one_of(scale_drop))
     comm %>% 
        left_join(traits, by = c(scale_keep, taxon_col)) %>% 
        group_by_at(vars(one_of(c(scale_keep, taxon_col, trait_col, maybe_missing(other_col))))) %>%
      rename(weight = !!abundance_col) %>% 
      mutate(weight = weight/n()) %>% 
#      group_by(!!!keep_scale, !!!trait_col, !!!other_col) %>% 
      mutate(level = factor(scale_hierarchy[scale_level], levels = scale_hierarchy, ordered = TRUE))
  }) %>% 
    filter(!is.na(!!!value_col)) %>% 
    filter(level == min(.data$level)) %>% 
    group_by_at(.vars = vars(one_of(c(scale_hierarchy, trait_col, other_col))))

  class(out) <- c(class(out), "imputed_traits")
  out
  }




#' Bootstrap traits
#' @description Bootstrap impute traits 
#' @param traits_comm imputed trait and community data in long format
#' @param nrep number of bootstrap replicates
#' @param sample_size number of plants per sample
#' @description 
#' 
#' @return 
#' 
#' @importFrom stats var 
#' @importFrom e1071 skewness kurtosis
#' @importFrom magrittr %>%
#' @importFrom dplyr sample_n group_by summarise
#' @importFrom purrr map_df
#' @importFrom rlang !!! maybe_missing

trait_np_bootstrap <- function(traits_comm, nrep = 100, sample_size = 200){  
#  stopifnot(class(traits_com) == "imputed_traits")
  
  bootstrapMoments_All <- map_df(1:nrep, ~sample_n(traits_comm, size = sample_size,  replace = TRUE, weight = weight), .id = "n") %>% 
    group_by(n, add = TRUE) %>% 
    # get all the happy moments
    summarise(mean = mean(Value), variance = var(Value), skewness = skewness(Value), kurtosis = kurtosis(Value))
  
  return(bootstrapMoments_All)
}
