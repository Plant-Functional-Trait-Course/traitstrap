#' Summarise Bootstrap traits
#' @description Bootstraptraits 
#' @param BootstrapMoments trait moments from trait_np_bootstrap

#' @description 
#' 
#' @return tibble with the grouping variables and the mean of each moment (+/- 1SD)
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr n group_by any_of summarise across
#' @importFrom rlang .data
#' @export

trait_summarise_boot_moments <- function(BootstrapMoments){
  attrib <- attr(BootstrapMoments, "attrib")
  groups <- c(attrib$scale_hierarchy, 
              attrib$trait_col,
              attrib$other_col)
  
  # calculate means of moments 
  sBootstrapMoments <- BootstrapMoments %>% 
    group_by(across(any_of(groups))) %>% 
    summarise(
      n = n(),
      Mean = mean(.data$mean),
      CIlow.mean = .data$Mean - sd(.data$mean),
      CIhigh.mean = .data$Mean + sd(.data$mean),
      
      Var = mean(.data$variance),
      CIlow.var = .data$Var - sd(.data$variance),
      CIhigh.var = .data$Var + sd(.data$variance),
      
      Skew = mean(.data$skewness),
      CIlow.skew = .data$Skew - sd(.data$skewness),
      CIhigh.skew = .data$Skew + sd(.data$skewness),
      
      Kurt = mean(.data$kurtosis),
      CIlow.kurt = .data$Kurt - sd(.data$kurtosis),
      CIhigh.Kurt = .data$Kurt + sd(.data$kurtosis)
    )
  
  return(sBootstrapMoments)
}
