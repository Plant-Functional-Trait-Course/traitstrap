
######################################################

#' Summarise Bootstrap traits nonparametrically
#' @description Find the mean and confidence interval for each moment 
#' @param bootstrap_moments trait moments from trait_np_bootstrap or trait_bootstrap_parametric
#' @param CI Desired non-parametric confidence level.  Default is 0.95
#' @return tibble with the grouping variables and the mean, low, and high CIs
#' @importFrom magrittr %>%
#' @importFrom dplyr n group_by any_of summarise across rename
#' @importFrom rlang .data
#' @examples 
#' data(community)
#' data(trait)
#' imputed_traits <-trait_impute(comm = community, traits = trait,
#'                  scale_hierarchy = c("Site", "PlotID"),
#'                  taxon_col = "Taxon", value_col = "Value",
#'                  trait_col = "Trait", abundance_col = "Cover")
#' boot_traits <- trait_np_bootstrap(imputed_traits)
#' trait_summarise_boot_moments2(boot_traits)
#' @export
trait_summarise_boot_moments2 <- function(bootstrap_moments,CI = 0.95){
  
 #First, check if the number of replicates makes sense with the CI
    nreps <- max(as.integer(unique(bootstrap_moments$n)))
    tail_size <- ((1-CI)/2)*nreps
  
  #Check that the tail size is a whole number.  Have to hack this a bit because floating point issues I imagine
    if(round(tail_size,digits = 10) != round(tail_size,1)){
    stop("CI specified not compatible with the number of replicates.")
      }
  
  #
  
  
    attrib <- attr(bootstrap_moments, "attrib")
    groups <- c(as.character(attrib$scale_hierarchy), attrib$trait_col, 
                attrib$other_col)
    
    if (!is.null(attrib$treatment_col)) {
      groups <- c(groups, paste0(attrib$treatment_col, "_comm"))
    }
    
      bootstrap_moments %>% ungroup() %>% 
      group_by(across(any_of(groups))) %>% rename(MEAN = .data$mean) %>%
      summarise(n = n(),
                mean = mean(.data$MEAN), 
                ci_low_mean = sort(.data$MEAN)[tail_size+1],
                ci_high_mean = sort(.data$MEAN,decreasing = T)[tail_size+1],
                var=mean(.data$variance), 
                ci_low_var = sort(.data$variance)[tail_size+1],
                ci_high_var = sort(.data$variance,decreasing = T)[tail_size+1],
                skew=mean(.data$skewness), 
                ci_low_skew = sort(.data$skewness)[tail_size+1],
                ci_high_skew = sort(.data$skewness,decreasing = T)[tail_size+1],
                kurt=mean(.data$kurtosis), 
                ci_low_kurt = sort(.data$kurtosis)[tail_size+1],
                ci_high_kurt = sort(.data$kurtosis,decreasing = T)[tail_size+1],
                ) ->  summ_bootstrap_moments

      return(summ_bootstrap_moments)

}#fx

