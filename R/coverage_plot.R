#' Coverage plot of Imputed Traits
#' @description Function calculates the trait coverage of the community 
#' for each level of the sampling hierarchy and makes a barplot.
#' @param imputed_traits output from trait_np_bootstrap or trait_parametric_bootstrap
#' funtion.
#' @param other_col_how what to do with the other columns in other data.
#'  Options are to `filter` by one of the columns, add them to the x-`axis`,
#' `facet` by them, or to `ignore`.
#' @param ... optional filters for use with `other_col_how = "filter"`
#' @description Shows at which level the data are coming from in each plot.
#'
#' @return a ggplot
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise distinct bind_cols ungroup
#' @importFrom rlang .data
#' @importFrom ggplot2 autoplot ggplot geom_col facet_wrap aes
#' @importFrom ggplot2 scale_y_continuous labs fortify facet_grid
#' @importFrom stats reformulate
#' @importFrom glue glue
#' @examples
#' require("ggplot2")
#' data(community)
#' data(trait)
#' imputed_traits <-trait_impute(comm = community, traits = trait,
#'                  scale_hierarchy = c("Site", "PlotID"),
#'                  taxon_col = "Taxon", value_col = "Value",
#'                  trait_col = "Trait", abundance_col = "Cover")
#' autoplot(imputed_traits)
#' @export

autoplot.imputed_trait <- function(imputed_traits, other_col_how, ...) {
  #get scale_hierarchy and concatenate to make an ID
  attrib <- attr(imputed_traits, "attrib")
  
  if (length(attrib$other_col) > 0 & missing(other_col_how)) {
    stop(
      glue("Autoplot needs to know what to do with the other_col. Options are:
            axis: add to x-axis labels
            facet: add facets in the plot
            filter: filter out one value of interest
            ignore: ignore")
    )
  }
  
  if (!missing(other_col_how)) {
    other_col_how <- match.arg(
      arg = other_col_how, 
      choices = c("axis", "facet", "filter", "ignore")
    )
  }
  
  imputed_traits_summary <- fortify.imputed_trait(imputed_traits, other_col_how, ...)
  
  plot <- ggplot(imputed_traits_summary, aes(x = .data$.id,
                                     y = .data$s, fill = .data$level)) +
    geom_col() +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Plot", y = "Proportion of cover", fill = "Data source")
  
  # add facets
  if (!missing(other_col_how) && other_col_how == "facet"){
    plot <- plot + 
      facet_grid(reformulate(attrib$trait_col, attrib$other_col))
  } else {
   plot <- plot + facet_wrap(~ .data[[attrib$trait_col]])
  }
  
  return(plot)
}
