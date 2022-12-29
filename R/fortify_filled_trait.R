#' Fortify filled Traits
#' @param filled_traits filled traits from trait_fill
#' @param other_col_how what to do with the other columns in other data.
#'  Options are to `filter` by one of the columns, add them to the x-`axis`,
#'  `facet` by them, or to `ignore`.
#' @param ... optional filters for use with `other_col_how = "filter"`
#' @description Calculates at which level the data are coming from in each plot.
#'
#' @return a tibble
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise distinct bind_cols ungroup any_of
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom ggplot2 fortify
#' @export


fortify.filled_trait <- function(filled_traits, other_col_how, ...) {
  # get scale_hierarchy and concatenate to make an ID
  attrib <- attr(filled_traits, "attrib")
  scale_hierarchy <- as.character(attrib$scale_hierarchy)
  scale_hierarchy <- scale_hierarchy[scale_hierarchy != "global"]

  if (!missing(other_col_how)) {
    other_col_how <- match.arg(
      arg = other_col_how,
      choices = c("axis", "facet", "filter", "ignore")
    )
    if (other_col_how == "filter") {
      filled_traits <- filled_traits %>% filter(...)
    }
    if (other_col_how == "axis") {
      scale_hierarchy <- c(scale_hierarchy, attrib$other_col)
    }
  }

  id <- filled_traits %>%
    ungroup() %>%
    select(any_of(scale_hierarchy)) %>%
    apply(1, paste, collapse = "_")

  filled_traits_summary <- filled_traits %>%
    bind_cols(.id = id) %>%
    group_by(
      .data$.id,
      .data$level,
      .data[[attrib$trait_col]], .data[[attrib$taxon_col]],
      .add = TRUE
    ) %>%
    distinct(.data[[attrib$abundance_col]], .keep_all = TRUE) %>%
    summarise(s = sum(.data[[attrib$abundance_col]]) / .data$sum_abun) %>%
    summarise(s = sum(.data$s))

  attr(filled_traits_summary, "attrib") <- attrib

  return(filled_traits_summary)
}
