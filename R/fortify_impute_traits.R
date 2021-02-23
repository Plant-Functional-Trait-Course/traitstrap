#' Fortify Imputed Traits
#' @param imputed_traits trait moments from trait_np_bootstrap
#' @description Calculates at which level the data are coming from in each plot.
#'
#' @return a tibble
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise distinct bind_cols ungroup any_of
#' @importFrom rlang .data
#' @importFrom ggplot2 fortify
#' @export


fortify.imputed_trait <- function(imputed_traits) {
  #get scale_hierarchy and concatenate to make an ID
  attrib <- attr(imputed_traits, "attrib")
  scale_hierarchy <- attrib$scale_hierarchy
  scale_hierarchy <- scale_hierarchy[scale_hierarchy != "global"]
  id <- imputed_traits %>%
    ungroup() %>%
    select(any_of(scale_hierarchy)) %>%
    apply(1, paste, collapse = "_")

  imputed_traits_summary <- imputed_traits %>%
    bind_cols(.id = id) %>%
    group_by(
      .data$.id, .data$level,
      .data[[attrib$trait_col]], .data[[attrib$taxon_col]],
      .add = TRUE) %>%
    distinct(.data[[attrib$abundance_col]], .keep_all = TRUE) %>%
    summarise(s = sum(.data[[attrib$abundance_col]]) / .data$sum_abun) %>%
    summarise(s = sum(.data$s))

  attr(imputed_traits_summary, "attrib") <- attrib

  return(imputed_traits_summary)
}
