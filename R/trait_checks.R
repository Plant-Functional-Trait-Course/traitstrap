# ' Checker functions
#' @importFrom glue glue glue_collapse

trait_check_scales <- function(scale_hierarchy, comm_names) {
  if (!all(scale_hierarchy %in% comm_names)) {
    bad_scales <- glue_collapse(
      x = scale_hierarchy[!scale_hierarchy %in% comm_names],
      sep = ", ", last = ", and "
    )
    stop(glue("scale_hierarchy levels {bad_scales} not in names(comm)"))
  }
}

trait_check_taxon_col <- function(taxon_col, comm_names, traits_names) {
  if (!all(taxon_col %in% comm_names)) {
    bad_taxon <- glue_collapse(
      x = taxon_col[!taxon_col %in% comm_names],
      sep = ", ", last = ", and "
    )
    stop(glue("taxon_col {bad_taxon} not in names(comm)"))
  }

  if (!all(taxon_col %in% traits_names)) {
    bad_taxon <- glue_collapse(
      x = taxon_col[!taxon_col %in% traits_names],
      sep = ", ", last = ", and "
    )
    stop(glue("taxon_col {bad_taxon} not in names(traits)"))
  }
}

trait_check_col <- function(col, names) {
  if (!(length(col) == 1 && col %in% names)) {
    col_name <- deparse(substitute(col))
    df_name <- deparse(substitute(names)) |>
      sub(pattern = "_names", replacement = "")
    stop(glue("{col_name} '{trait_col}' not in names({df_name})"))
  }
}

trait_check_other_col <- function(other_col, comm_names, traits_names) {
  if (!all(other_col %in% comm_names)) {
    bad_other <- glue_collapse(
      x = other_col[!other_col %in% comm_names],
      sep = ", ", last = ", and "
    )
    stop(glue("other_col levels {bad_other} not in names(comm)"))
  }

  if (any(other_col %in% traits_names)) {
    bad_other <- glue_collapse(
      x = other_col[!other_col %in% traits_names],
      sep = ", ", last = ", and "
    )
    warning(glue("other_col levels {bad_other} are in  names(traits). \\
                 These columns will be removed from traits"))
  }
}

trait_check_treatment <- function(
    treatment_col, comm_names, traits_names,
    comm, traits, treatment_level, scale_hierarchy) {
  trait_check_col(treatment_col, comm_names)
  trait_check_col(treatment_col, traits_names)

  # treatment col must be factor
  if (!is.factor(comm[[treatment_col]])) {
    stop(glue("treatment_col {treatment_col} is not a factor in comm"))
  }
  if (!is.factor(traits[[treatment_col]])) {
    stop(glue("treatment_col {treatment_col} is not a factor in traits"))
  }

  # treatment col levels must be identical in traits and comm
  if (!identical(
    levels(traits[[treatment_col]]),
    levels(comm[[treatment_col]])
  )
  ) {
    stop("treatment_col has have different levels in comm and traits")
  }

  # check treatment_level is valid
  if (is.null(treatment_level)) {
    stop("treatment_level must be specified when treatment_col is used")
  }
  if (!treatment_level %in% scale_hierarchy) {
    stop("treatment_level must be in scale_hierarchy")
  }
}

traits_check_global <- function(traits_names) {
  if (any(traits_names == "global")) {
    stop(glue("Cannot add global column as column called\\
                global already exists. Maybe set 'global = FALSE'"))
  }
}
