#' Bootstrap traits parametrically
#' @description Bootstrap fitted distributions
#' @param fitted_distributions
#' Fitted distribution object returned by fit_trait_distributions
#' @param nrep number of bootstrap replicates
#' @param sample_size number of plants per sample
#' @description
#'
#' @return a tibble
#'
#' @importFrom stats var
#' @importFrom e1071 skewness kurtosis
#' @importFrom magrittr %>%
#' @importFrom dplyr slice_sample group_by summarise rowwise
#' @importFrom tidyr unnest
#' @importFrom purrr map_dfr
#' @examples
#' data(community)
#' data(trait)
#' imputed_traits <-trait_impute(comm = community, traits = trait,
#'                  scale_hierarchy = c("Site", "PlotID"),
#'                  taxon_col = "Taxon", value_col = "Value",
#'                  trait_col = "Trait", abundance_col = "Cover")
#' fitted_distributions <- fit_trait_distributions(imputed_traits =
#' imputed_traits, distribution_type = "normal")
#' parametric_dists <- trait_parametric_bootstrap(fitted_distributions =
#' fitted_distributions, nrep = 10, sample_size = 100)
#' moment_summary <-
#' trait_summarise_boot_moments(bootstrap_moments = parametric_dists,
#' parametric = FALSE)
#' @export
trait_parametric_bootstrap <- function(fitted_distributions,
                                       nrep = 100,
                                       sample_size = 200) {


  #Pull useful information from imputed traits object
  trait_col <- attributes(fitted_distributions)$attrib$trait_col
  abundance_col <- attributes(fitted_distributions)$attrib$abundance_col
  taxon_col <- attributes(fitted_distributions)$attrib$taxon_col
  scale_hierarchy <- attributes(fitted_distributions)$attrib$scale_hierarchy
  attrib <- attr(fitted_distributions, "attrib")

  bootstrap_moments <- map_dfr(
    1:nrep,
    ~ {

      fitted_distributions %>%
        group_by_at(c(as.character(scale_hierarchy), trait_col)) %>%
        slice_sample(n = sample_size,
                     replace = T,
                     weight_by = .data[[abundance_col]]) %>%
        group_by_at(c(as.character(scale_hierarchy),
                      trait_col, taxon_col,
                      "parm1", "parm2", "distribution_type")) %>%
        summarise(n_drawn = n(), .groups = "keep") %>%
        rowwise() %>%
        mutate(draw_value =
                 list(distribution_handler(parm1 = parm1,
                                      parm2 = parm2,
                                      n = n_drawn,
                                      type = distribution_type))) %>%
        group_by_at(c(as.character(scale_hierarchy), trait_col)) %>%
        unnest(draw_value) %>%
        summarise(mean = mean(draw_value),
                  variance = var(draw_value),
                  skewness = skewness(draw_value),
                  kurtosis = kurtosis(draw_value),
                  .groups = "keep")

    },
    .id = "n"
  )

  attr(bootstrap_moments, "attrib") <- attrib

  # make bootstrap_moments an ordinary tibble
  class(bootstrap_moments) <-
    class(bootstrap_moments)[!class(bootstrap_moments) ==
                               "parametric_distributions"]
  return(bootstrap_moments)


}# end fx
