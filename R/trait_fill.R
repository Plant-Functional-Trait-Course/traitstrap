#' fill traits
#' @description A function for trait filling using a hierarchical sampling
#' design, which allows to account for incomplete trait collections,
#' traits from different spatial or temporal levels
#' (i.e. local traits vs. databases) and experimental designs.
#' @param comm a dataframe in long format with community data
#' @param traits a dataframe in long format with trait data
#' @param scale_hierarchy character vector of the sampling hierarchy from large
#' to small (e.g. site/block/plot)
#' @param taxon_col character; name of taxon column in comm and traits.
#' Can be a vector (e.g. "species", "genus"), in which case if traits cannot be
#' filled for the first taxon column, subsequent columns will be used in order.
#' @param trait_col character; name of trait name column in traits
#' @param value_col character; name of trait value column in traits
#' @param abundance_col character; name of species abundance column in comm
#' @param other_col name of other grouping columns in comm
#' @param treatment_col character; optional name of treatment_col in comm and
#' traits. Must refer to a factor where first level is control.
#' @param treatment_level character; optional name of scale_hierarchy at which
#' treatment should be filtered
#' @param global logical; calculate traits at global scale.
#' Must not be a column called global in the traits data.
#' @param keep_all logical; keep trait data at all available levels
#' or just finest scale available
#' @param min_n_in_sample numeric; minimum number in sample with traits to
#' accept before searching for traits higher up the hierarchy. The default is 5.
#' @param complete_only logical; use only leaves with a full set of traits.
#' Set to TRUE when imputing for \code{trait_multivariate_bootstrap}
#' @param leaf_id character; unique leaf identifiers.
#' Only needed when \code{complete_only} is TRUE.
#'
#'
#' @details The function uses a hierarchical sampling design, which allows it to
#' account for incomplete trait collections, traits from different spatial or
#' temporal levels (i.e. local traits vs. databases) and/or
#' experimental designs.
#'
#' With scale_hierarchy you can define the levels at which the traits have been
#' collected and their order starting with the highest level
#' (e.g. global database, region, site, block, plot).
#'
#' `trait_fill()` will choose if available a trait value from the lowest
#' level, i.e. species X from plot A and if no trait is available from that
#' level, it will move up the hierarchy and choose a trait from species X from
#' plot B at the same site.
#' If there is no trait available from species X in the same site, it
#' will choose a trait value from another site.
#'
#' The argument `min_n_in_samples` allows users to define the minimum number in
#' sample at each level for the trait filling.
#' If the minimum number is not reached,
#' trait values from the next level will also be selected, to avoid sampling the
#' same individual several times, which will result in unrealistic variances.
#' The default value is 5.
#'
#' In the `other_col` argument other grouping variables in the community dataset
#' can be defined and will be kept after the trait filling step.
#'
#' Traitstrap also allows to include taxonomy and experimental design in the
#' trait filling step.
#'
#' With taxon_col a hierarchy for the taxonomy can be defined.
#' If traits for a specific
#' species are not available, traits from next level,
#' e.g. the genus will be selected.
#' For this a list of the taxonomic hierarchy has to be defined (e.g. "taxon",
#' "genus", "family").
#'
#' The argument treatment_col allows to incorporate an experimental design where
#' traits are selected from the same experimental treatment or the first factor
#' level, which is assumed to be the control.
#' Therefore, it is important to order
#' the levels of a treatment in the right order, i.e. the first level has to be
#' the control.
#' If you have two or more treatments and you want filling to be done only
#' within a treatment, and not from a treatment and the control, then
#' make the first level of the factor a level that is not in the data.
#' The filling step can be defined at certain level using the
#' treatment_level argument. Depending on the experimental design trait filling
#' should occur a certain level, e.g. block or site.

#'
#' @return a tibble with extra class `filled_trait`. 
#'
#' @importFrom stats sd var weighted.mean
#' @importFrom dplyr select any_of all_of mutate group_by filter left_join n
#' @importFrom dplyr inner_join across ungroup
#' @importFrom purrr map list_rbind
#' @importFrom rlang !!! !! .data
#' @importFrom glue glue glue_collapse
#' @importFrom tibble lst
#' @examples
#' data(community)
#' data(trait)
#' filled_traits <- trait_fill(
#'   comm = community, traits = trait,
#'   scale_hierarchy = c("Site", "PlotID"),
#'   taxon_col = "Taxon", value_col = "Value",
#'   trait_col = "Trait", abundance_col = "Cover"
#' )
#' @export

trait_fill <- function(comm,
                       traits,
                       scale_hierarchy = c(
                         "Country", "Site",
                         "BlockID", "PlotID"
                       ),
                       global = TRUE,
                       taxon_col = "taxon", trait_col = "trait",
                       value_col = "Value", abundance_col = "Cover",
                       treatment_col = NULL, treatment_level = NULL,
                       other_col = character(0),
                       keep_all = FALSE,
                       min_n_in_sample = 5,
                       complete_only = FALSE,
                       leaf_id) {
  #### sanity checks on input (are columns present etc) ####
  comm_names <- names(comm)
  traits_names <- names(traits)

  # check data have all scales in scale_hierarchy
  trait_check_scales(scale_hierarchy, comm_names)

  # check taxon_col is valid
  trait_check_taxon_col(taxon_col, comm_names, traits_names)

  # check trait_col is valid
  trait_check_col(trait_col, traits_names)

  # check value_col is valid
  trait_check_col(value_col, traits_names)

  # check abundance_col is valid
  trait_check_col(abundance_col, comm_names)

  # if used, check treatment_col is valid
  if (is.null(treatment_col)) {
    use_treat <- FALSE
  } else {
    use_treat <- TRUE
    trait_check_treatment(treatment_col, comm_names, traits_names, comm, traits, treatment_level, scale_hierarchy)
  }

  # check other_cols are valid
  trait_check_other_col(other_col, comm_names, traits_names)

  # remove and other col in traits
  traits <- traits |> select(-any_of(other_col))

  #### prep ####
  # add global to scale_hierarchy if necessary
  if (isTRUE(global)) {
    # check not already a "global" column
    traits_check_global(traits_names)

    # add global column to traits
    traits <- traits |> mutate(global = "global")
    comm <- comm |> mutate(global = "global")
    scale_hierarchy <- c("global", scale_hierarchy)
  }

  # remove NA trait values
  traits <- traits |> filter(!is.na(.data[[value_col]]))

  ## remove leaves with incomplete set of traits if complete_only
  if (isTRUE(complete_only)) {
    n_traits <- n_distinct(traits[[trait_col]])
    if (missing(leaf_id) || !leaf_id %in% names(traits)) {
      stop("leaf_id is missing or not in the traits data")
    }
    traits <- traits |>
      group_by(.data[[leaf_id]]) |>
      mutate(.n = n())
    if (any(traits$.n > n_traits)) {
      stop("Leaves can only have one measurement per trait.
       Check leaf ID are unique.")
    }
    traits <- traits |>
      filter(.data$.n == n_traits) |>
      select(-.data$.n)
  }

  ## check for NA in abundance
  if (any(is.na(comm[[abundance_col]]))) {
    stop(glue("cannot have NA in the {abundance_col}\\
              column of the community data"))
  }

  ##### routine if length(taxon_col) > 1
  if (length(taxon_col) > 1) {
    result <- trait_fill_multi_level(
      comm = comm,
      traits = traits,
      scale_hierarchy = scale_hierarchy,
      global = global,
      taxon_col = taxon_col,
      trait_col = trait_col,
      value_col = value_col,
      abundance_col = abundance_col,
      treatment_col = treatment_col,
      treatment_level = treatment_level,
      other_col = other_col,
      keep_all = keep_all,
      min_n_in_sample = min_n_in_sample
    )
    return(result)
  }

  # calculate plot scale sum of abundances
  comm <- comm |>
    group_by(across(all_of(c(as.character(scale_hierarchy), other_col)))) |>
    # calculate sum abundance
    mutate(sum_abun = sum(.data[[abundance_col]]))

  # make ordered factor of scale hierarchy
  scale_hierarchy <- factor(scale_hierarchy,
    levels = rev(scale_hierarchy),
    ordered = TRUE
  )

  #### iterate over grouping hierarchy####
  out <- scale_hierarchy |>
    map(~ {
      scale_level <- .x

      # drop scales from the hierarchy
      scale_drop <- as.character(
        scale_hierarchy[scale_hierarchy < scale_level]
      )

      # scales to keep
      scale_keep <- as.character(
        scale_hierarchy[scale_hierarchy >= scale_level]
      )

      traits <- traits |> select(-any_of(scale_drop))
      result <- comm |>
        ungroup() |> # avoids dplyr problem #5473
        # join to traits
        inner_join(traits,
          by = c(scale_keep, taxon_col),
          suffix = c("_comm", "_trait"),
          relationship = "many-to-many"
        ) |>
        # group by kept scales
        group_by(
          across(all_of(c(
            scale_keep, other_col,
            trait_col, taxon_col, "sum_abun"
          )))
        )
      # filter if using treatment_col
      if (use_treat) {
        if (treatment_level == scale_level) {
          col_comm <- paste0(treatment_col, "_comm")
          col_trait <- paste0(treatment_col, "_trait")
          result <- result |>
            filter( # same treatment OR first level (must be control)
              .data[[col_comm]] == .data[[col_trait]] |
                .data[[col_trait]] == levels(.data[[col_trait]])[1]
            )
        }
      }

      # calculate weights
      result <- result |>
        ungroup() |>
        group_by(across(all_of(c(
          as.character(scale_hierarchy),
          taxon_col, trait_col, other_col
        )))) |>
        mutate(
          n_sample = n(),
          weight = .data[[abundance_col]] / .data$n_sample,
          level = scale_level
        )

      result
    }) |>
    list_rbind() # end of iterate over grouping hierarchy

  # check some output
  if (nrow(out) == 0) {
    stop("No output.
       Check taxon names and plot/site names in the community and
       trait data frames match.")
  }

  # get max number available
  out <- out |>
    ungroup() |>
    group_by(across(all_of(c(
      as.character(scale_hierarchy),
      taxon_col, trait_col, other_col
    )))) |>
    mutate(max_n_in_sample = max(.data$n_sample))

  #### filter out lowest good level of hierarchy for each taxon & trait ####
  if (!keep_all) { # keep only finest scale trait data available
    out <- out |>
      ungroup() |>
      group_by(across(all_of(c(
        as.character(scale_hierarchy),
        taxon_col, trait_col, other_col
      )))) |>
      filter(
        # group has more than minimum number in sample
        .data$n_sample >= min_n_in_sample |
          # OR maximum number when max is lower than required
          (.data$n_sample == .data$max_n_in_sample &
            .data$max_n_in_sample <= min_n_in_sample)
      ) |>
      # filter lowest level available by group
      filter(.data$level == min(.data$level))
  }

  out <- out |>
    ungroup() |> # avoids dplyr problem #5473
    group_by(across(all_of(c(
      as.character(scale_hierarchy),
      trait_col, other_col
    ))))

  # add treatment group if used
  if (!is.null(treatment_col)) {
    out <- out |>
      group_by(.data[[paste0(treatment_col, "_comm")]], .add = TRUE)
  }

  # set arguments as attributes so next functions have access to them
  attrib <- tibble::lst(
    scale_hierarchy, taxon_col, trait_col,
    value_col, abundance_col, other_col,
    treatment_col, treatment_level
  )
  attr(out, "attrib") <- attrib

  class(out) <- c("filled_trait", class(out))

  return(out)
}
