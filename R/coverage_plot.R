#' Coverage plot of Imputed Traits
#' @param imputed_traits trait moments from trait_np_bootstrap
#' @description Shows at which level the data are coming from in each plot.
#'
#' @return a ggplot
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise distinct bind_cols ungroup
#' @importFrom rlang .data
#' @importFrom ggplot2 autoplot ggplot geom_col facet_wrap aes
#' @importFrom ggplot2 scale_y_continuous labs fortify
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
  
  if (length(attrib$other_col) > 0)

  imputed_traits_summary <- fortify.imputed_trait(imputed_traits)
  ggplot(imputed_traits_summary, aes(x = .data$.id,
                                     y = .data$s, fill = .data$level)) +
    geom_col() +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Plot", y = "Proportion of cover", fill = "Data source") +
    facet_wrap(~ .data[[attrib$trait_col]])
}
