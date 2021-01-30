#' Summarise Bootstrap traits
#' @description Bootstraptraits 
#' @param BootstrapMoments trait moments from trait_np_bootstrap

#' @description 
#' 
#' @return tibble with the grouping variables and the mean of each moment (+/- 1SD)
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr n group_by any_of summarise across rename
#' @importFrom rlang .data
#' @export

trait_summarise_boot_moments <- function(BootstrapMoments){
  attrib <- attr(BootstrapMoments, "attrib")
  groups <- c(as.character(attrib$scale_hierarchy), 
              attrib$trait_col,
              attrib$other_col)
  
  # calculate means of moments 
  sBootstrapMoments <- BootstrapMoments %>% 
    ungroup() %>% 
    group_by(across(any_of(groups))) %>% 
    rename(MEAN = .data$mean) %>% #capitalise to simplify code below
    summarise(
      n = n(),
      mean = mean(.data$MEAN),
      ci_low_mean = .data$mean - sd(.data$MEAN),
      ci_high_mean = .data$mean + sd(.data$MEAN),
      
      var = mean(.data$variance),
      ci_low_var = .data$var - sd(.data$variance),
      ci_high_var = .data$var + sd(.data$variance),
      
      skew = mean(.data$skewness),
      ci_low_skew = .data$skew - sd(.data$skewness),
      ci_high_skew = .data$skew + sd(.data$skewness),
      
      kurt = mean(.data$kurtosis),
      ci_low_kurt = .data$kurt - sd(.data$kurtosis),
      ci_high_Kurt = .data$kurt + sd(.data$kurtosis)
    )
  
  return(sBootstrapMoments)
}
