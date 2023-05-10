#' Which taxa lack traits
#' @description Function gives overview of which taxa are missing traits.
#' @param filled_trait output of trait_fill function.
#' @param comm community data
#' @returns A `tibble` with columns
#' #' \itemize{
#'  \item{Taxon }{Species names (actual name depends on `taxon_col`
#'               argument to `trait_fill()`)}
#'  \item{max_abun }{Maximum abundance of that taxa.
#'                Be more concerned about taxa missing traits with high abundances.}
#'  \item{n }{Number of occurrences of the taxon. 
#'    Be more concerned about taxa missing traits with many occurrences.}
#'  \item{n_traits }{Number of traits for each species.
#'    Ideally all should equal the number of traits you have measured.}
#'  }
#' @importFrom dplyr left_join ungroup group_by across all_of summarise distinct
#' @importFrom rlang .data
#' @examples
#' data(community)
#' data(trait)
#' filled_traits <- trait_fill(
#'   comm = community, traits = trait,
#'   scale_hierarchy = c("Site", "PlotID"),
#'   taxon_col = "Taxon", value_col = "Value",
#'   trait_col = "Trait", abundance_col = "Cover"
#' )
#' trait_missing(filled_traits, community)
#' @export
trait_missing <- function(filled_trait, comm) {
  attrib <- attr(filled_trait, "attrib")

  filled_trait <- filled_trait |> 
    group_by(.data[[attrib$taxon_col]]) |> 
    distinct(.data[[attrib$trait_col]])
  
  comm |>
    group_by(.data[[attrib$taxon_col]]) |>
    summarise(
      max_abun = max(.data[[attrib$abundance_col]]),
      n = n()
      ) |>
    left_join(
      filled_trait,
      by = attrib$taxon_col,
      suffix = c("", "_trait")
    ) |>
    group_by(.data[[attrib$taxon_col]], .data$max_abun, .data$n) |>
    summarise(n_traits = sum(!is.na(.data[[attrib$trait_col]])))
}
