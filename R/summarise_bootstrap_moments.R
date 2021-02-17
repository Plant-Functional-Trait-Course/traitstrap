#' Summarise Bootstrap traits
#' @description Find the mean and confidence interval for each moment 
#' @param bootstrap_moments trait moments from trait_np_bootstrap
#' @param sd_mult Number of standard deviations around each moment,
#' defaults to one

#' @description
#'
#' @return tibble with the grouping variables and the mean of each moment
#'  (+/- sd_mult * SD)
#'
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
#' trait_summarise_boot_moments(boot_traits)
#' @export

trait_summarise_boot_moments <- function(bootstrap_moments, sd_mult = 1) {
  attrib <- attr(bootstrap_moments, "attrib")
  groups <- c(as.character(attrib$scale_hierarchy),
              attrib$trait_col,
              attrib$other_col)
  
  #add treatment_col to groups if used
  if (!is.null(attrib$treatment_col)) {
    groups <- c(groups, paste0(attrib$treatment_col, "_comm"))
  }


  # calculate means of moments
  summ_bootstrap_moments <- bootstrap_moments %>%
    ungroup() %>%
    group_by(across(any_of(groups))) %>%
    rename(MEAN = .data$mean) %>% #capitalise to simplify code below
    summarise(
      n = n(),
      mean = mean(.data$MEAN),
      ci_low_mean = .data$mean - sd(.data$MEAN) * sd_mult,
      ci_high_mean = .data$mean + sd(.data$MEAN) * sd_mult,

      var = mean(.data$variance),
      ci_low_var = .data$var - sd(.data$variance) * sd_mult,
      ci_high_var = .data$var + sd(.data$variance) * sd_mult,

      skew = mean(.data$skewness),
      ci_low_skew = .data$skew - sd(.data$skewness) * sd_mult,
      ci_high_skew = .data$skew + sd(.data$skewness) * sd_mult,

      kurt = mean(.data$kurtosis),
      ci_low_kurt = .data$kurt - sd(.data$kurtosis) * sd_mult,
      ci_high_kurt = .data$kurt + sd(.data$kurtosis) * sd_mult 
    )

  return(summ_bootstrap_moments)
}
