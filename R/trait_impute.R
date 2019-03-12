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
#' @importFrom dplyr select one_of mutate group_by filter left_join n group_by_at
#' @importFrom purrr map_df
#' @importFrom rlang !!! !! .data
#' @export
 
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
        group_by_at(vars(one_of(c(scale_keep, taxon_col, trait_col, other_col)))) %>%
      rename(weight = !!abundance_col) %>% 
      mutate(weight = weight/n()) %>% 
      mutate(level = factor(scale_hierarchy[scale_level], levels = scale_hierarchy, ordered = TRUE))
  }) %>% 
    filter(!is.na(!!!value_col)) %>% 
    filter(level == min(.data$level)) %>% 
    group_by_at(.vars = vars(one_of(c(scale_hierarchy, trait_col, other_col))))

  attr(out, "scale_hierarchy") <- scale_hierarchy
  attr(out, "taxon_col") <- taxon_col
  attr(out, "trait_col") <- trait_col
  attr(out, "value_col") <- value_col
  attr(out, "abundance_col") <- abundance_col
  attr(out, "other_col") <- other_col
  
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
#' @importFrom dplyr sample_n group_by summarise_at vars one_of
#' @importFrom purrr map_df
#' @export

trait_np_bootstrap <- function(traits_comm, nrep = 100, sample_size = 200){  
#  stopifnot(class(traits_com) == "imputed_traits")
  value_col <- attr(traits_comm, "value_col")
  bootstrapMoments <- map_df(1:nrep, ~sample_n(traits_comm, size = sample_size,  replace = TRUE, weight = weight), .id = "n") %>% 
    group_by(n, add = TRUE) %>% 
    # get all the happy moments
    summarise_at(.vars = vars(one_of(value_col)), .funs = list(mean = mean, variance = var, skewness = skewness, kurtosis = kurtosis))
  
  attr(bootstrapMoments, "scale_hierarchy") <- attr(traits_comm, "scale_hierarchy")
  attr(bootstrapMoments, "trait_col") <- attr(traits_comm, "trait_col") 
  attr(bootstrapMoments, "other_col") <- attr(traits_comm, "other_col")
  return(bootstrapMoments)
}


#' Summarise Bootstrap traits
#' @description Bootstraptraits 
#' @param BootstrapMoments trait moments from trait_np_bootstrap

#' @description 
#' 
#' @return 
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr n group_by summarise_at vars one_of group_by_at summarise
#' @importFrom rlang .data
#' @export

SummariseBootMoments <- function(BootstrapMoments){
  groups <- c(attr(BootstrapMoments, "scale_hierarchy"), 
              attr(BootstrapMoments, "trait_col"),
              attr(BootstrapMoments, "other_col"))
              
  # calculate means of moments 
  sBootstrapMoments <- BootstrapMoments %>% 
    group_by_at(vars(one_of(groups))) %>% 
    summarise(n = n(),
              Mean = mean(.data$mean), CIlow.mean = .data$Mean - sd(.data$mean), CIhigh.mean = .data$Mean + sd(.data$mean),
              Var = mean(.data$variance), CIlow.var = .data$Var - sd(.data$variance), CIhigh.var = .data$Var + sd(.data$variance),
              Skew = mean(.data$skewness), CIlow.skew = .data$Skew - sd(.data$skewness), CIhigh.skew = .data$Skew + sd(.data$skewness),
              Kurt = mean(.data$kurtosis), CIlow.kurt = .data$Kurt - sd(.data$kurtosis), CIhigh.Kurt = .data$Kurt + sd(.data$kurtosis)) 
  
  return(sBootstrapMoments)
}

