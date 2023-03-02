#' Which taxa lack traits
#' @description Function gives overview of which taxa are missing traits.
#' @param filled_trait output of trait_fill function.
#' @param comm community data
#' @importFrom dplyr left_join ungroup group_by across all_of summarise distinct
#' @importFrom rlang .data
#'
#' @export
trait_missing <- function(filled_trait, comm) {
  attrib <- attr(filled_trait, "attrib")

  comm |>
    left_join(filled_trait,
      by = attrib$taxon_col,
      suffix = c("", "_trait")
    ) |>
    ungroup() |>
    group_by(.data[[attrib$taxon_col]]) |>
    mutate(max_abun = max(.data[[attrib$abundance_col]])) |>
    ungroup() |>
    distinct(
      across(all_of(c(attrib$trait_col, attrib$taxon_col))),
      .data$max_abun
    ) |>
    group_by(.data[[attrib$taxon_col]], .data$max_abun) |>
    summarise(n_traits = n())
}
