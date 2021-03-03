#' Fit trait distributions
#' @description Fits distributions to imputed traits and
#' returns a tibble containing the fitted parameters.
#' @param imputed_traits imputed trait and community data
#' in long format
#' @param distribution_type the type of statistical distribution to use.
#' Character or named list. Currently accepts "normal","lognormal", and "beta".
#' @note Distribution type should be a single character (e.g. "normal")
#' or a named list (e.g. list(height="normal",mass="lognormal"))
#' @description
#'
#' @return a tibble
#'
#' @importFrom stats var
#' @importFrom e1071 skewness kurtosis
#' @importFrom magrittr %>%
#' @importFrom dplyr slice_sample group_by summarise
#' summarize select group_by_at
#' @importFrom purrr map_dfr
#' @importFrom stats var
#' @importFrom fitdistrplus fitdist

#' @examples
#' data(community)
#' data(trait)
#' imputed_traits <- trait_impute(comm = community, traits = trait,
#'                   scale_hierarchy = c("Site", "PlotID"),
#'                   taxon_col = "Taxon", value_col = "Value",
#'                   trait_col = "Trait", abundance_col = "Cover")
#' fitted_distributions <- trait_fit_distributions(
#'                         imputed_traits = imputed_traits,
#'                         distribution_type = "normal")
#'
#' @export
trait_fit_distributions <- function(imputed_traits,
                                    distribution_type = "normal") {

  #Check imputed traits
  if (! inherits(imputed_traits, "imputed_trait")) {
    stop("Imputed traits are not appropriately formatted.
    Please use trait_impute() ")
  }


  #Pull useful information from imputed traits object
  value_col <- attributes(imputed_traits)$attrib$value_col
  trait_col <- attributes(imputed_traits)$attrib$trait_col
  taxon_col <- attributes(imputed_traits)$attrib$taxon_col
  abundance_col <- attributes(imputed_traits)$attrib$abundance_col
  scale_hierarchy <- attributes(imputed_traits)$attrib$scale_hierarchy

  #If only a single type of distribution was supplied, apply it to all traits
  if (length(distribution_type) == 1) {
    distributions <- unique(as.data.frame(imputed_traits)[, trait_col])
    dist_list <- replicate(n = length(distributions), expr = distribution_type)
    names(dist_list) <- distributions
    distribution_type <- as.list(dist_list)
    rm(dist_list, distributions)
  } #if statement

  #Data checks

  #check distribution type is supported
  if (any(!unique(distribution_type) %in% c("beta", "normal", "lognormal"))) {
    x <- distribution_type[!distribution_type %in%
                            c("beta", "normal", "lognormal")]
    stop(glue("distribution type '{x}' not supported "))
  }


  #Beta checks
  if ("beta" %in% distribution_type) {

    #make sure every species x hierarchy combination has at least 2 data points
    beta_counts <- imputed_traits %>%
      filter(.data[[trait_col]] %in%
               names(distribution_type)[distribution_type == "beta"])

    if (any(beta_counts$n_leaves < 2)) {
    stop("Fitting a beta distrbution requires 2+ points per distribution.
         We suggest re-imputing traits with
         a minimum sample size of (at least) 2 traits.")
      }

    #check that values are between 0 and 1
    beta_vals <- imputed_traits %>%
      ungroup() %>%
      filter(.data[[trait_col]] %in%
               names(distribution_type)[distribution_type == "beta"]) %>%
      select(all_of(value_col))

    if (any(beta_vals > 1 | beta_vals < 0)) {
      stop("For a Beta distribution values must be between 0 and 1.")
    }


  }# end beta checks


  #lognormal checks
  if ("lognormal" %in% distribution_type) {

    ln_vals <- imputed_traits %>%
      ungroup() %>%
      filter(.data[[trait_col]] ==
               names(distribution_type)[distribution_type == "lognormal"]) %>%
      select(all_of(value_col))

    if (any(ln_vals <= 0)) {
      stop("For a lognormal distribution values must be positive.")
    }

  }#end lognormal checks

  # Main body

  distribution_parms <- imputed_traits %>%
    group_by(.data[[c(taxon_col)]],
             .data[[c(abundance_col)]], .data$n_leaves, .add = TRUE) %>%
    mutate(distribution_type = unlist(
      distribution_type[.data[[trait_col]]])) %>%
    summarize(quiet(get_dist_parms(data = .data[[value_col]],
                             distribution_type = unique(distribution_type))),
              .groups = "keep")


  #set arguments as attributes so next functions have access to them
  attr(distribution_parms, "attrib") <- attributes(imputed_traits)$attrib

  class(distribution_parms) <- c("parametric_distributions",
                                 class(distribution_parms))

  return(distribution_parms)

}# fx
