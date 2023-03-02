#' @importFrom purrr set_names map
#' @importFrom rlang .data
#' @importFrom dplyr filter mutate rename_with select bind_rows matches groups
#' @importFrom stringr str_remove

trait_fill_multi_level <- function(comm, traits,
                                   scale_hierarchy, global,
                                   taxon_col, trait_col,
                                   value_col, abundance_col,
                                   treatment_col, treatment_level,
                                   other_col, keep_all, min_n_in_sample) {
  # iterate over taxon_col
  result <- taxon_col |>
    set_names() |>
    map(~ trait_fill(
      comm = comm,
      traits = traits,
      scale_hierarchy = scale_hierarchy,
      global = FALSE, # global already added if needed
      taxon_col = .x,
      trait_col = trait_col,
      value_col = value_col,
      abundance_col = abundance_col,
      treatment_col = treatment_col,
      treatment_level = treatment_level,
      other_col = other_col,
      keep_all = keep_all,
      min_n_in_sample = min_n_in_sample
    ))

  # grab attributes
  result_attr <- attr(result[[1]], "attrib")

  result <- result |>
    map(select, -matches(paste0(taxon_col, "_trait"))) |>
    map(
      rename_with, str_remove,
      matches(paste0(taxon_col, "_comm")), "_comm"
    ) |>
    bind_rows(.id = "taxon_level")

  # filter lowest value of taxon hierarchy
  # need to add species to group
  current_groups <- groups(result) |> as.character()
  result <- result |>
    ungroup() |>
    mutate(taxon_level = factor(.data$taxon_level,
      levels = taxon_col, ordered = TRUE
    )) |>
    group_by(across(all_of(c(current_groups, taxon_col[1])))) |>
    filter(.data$taxon_level == min(.data$taxon_level)) |>
    ungroup() |>
    group_by(across(all_of(current_groups)))

  # reinstate attributes
  attr(result, "attrib") <- result_attr

  # add taxon_hierarchy to attributes
  attr(result, "taxon_hierarchy") <- taxon_col

  # add class
  class(result) <- c("filled_trait", class(result))
  # return
  return(result)
}
