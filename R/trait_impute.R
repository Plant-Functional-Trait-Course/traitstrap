#' Impute traits
#' @description A function for trait imputation using a hierarchical sampling design, 
#' which allows to account for incomplete trait collections, traits from different 
#' spatial or temporal levels (i.e. local traits vs. databases) and experimental designs.
#' @param comm a datafram in long format with community data
#' @param traits a dataframe in long format with trait data
#' @param scale_hierarchy character vector of the sampling hierarchy from large 
#' to small (e.g. site/block/plot)
#' @param taxon_col character; name of taxon column in comm and traits.
#' Can be a vector (e.g. "species", "genus"), in which case if traits cannot be 
#' imputed for the first taxon column, subsequent columns will be used in order.
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
#' @param min_n_in_sample numeric; minimum number in sample with traits to accept before
#' searching for traits higher up the hierarchy. The default is 5.
#'
#' @description
#' 
#' @details The function uses a hierarchical sampling design, which allows it to 
#' account for incomplete trait collections, traits from different spatial or temporal
#' levels (i.e. local traits vs. databases) and/or experimental designs. 
#' 
#' With scale_hierarchy you can define the levels at which the traits have been 
#' collected and their order starting with the highest level (e.g. global database, 
#' region, site, block, plot).
#' 
#' The trait_impute function will choose if available a trait value from the lowest 
#' level, i.e. species X from plot A and if no trait is available from that level, 
#' it will move up the hierarchy and choose a trait from species X from plot B at 
#' the same site. If there is no trait available from species X in the same site, it 
#' will choose a trait value from another site.
#' 
#' The argument min_n_in_samples allows users to define the minimum number in sample
#' at each level for the trait imputation. If the minimum number is not reached, trait
#' values from the next level will also be imputed, to avoid sampling the same 
#' individual several times, which will result in unrealistic variances. The default 
#' value is 5
#' 
#' In the other_col arumgent other grouping variable in the community dataset can be 
#' defined and will be kept after the trait imputation step.
#' 
#' Traitstrap also allows to include taxonomy and experimental design in the 
#' trait imputation step.
#' 
#' With taxon_col a hierarchy for the taxonomy can be defined. If traits for a specific
#' species are not available, traits from next level, e.g. the genus will be imputed.
#' For this a list of the taxonomic hierarchy has to be defined (e.g. "taxon", 
#' "genus", "family").
#' 
#' The argument treatment_col allows to incorporate an experimental design where 
#' traits are imputed from the same experimental treatment or the first factor 
#' level, which is assumed to be the control. Therefore, it is important to order 
#' the levels of a treatment in the right order, i.e. the first level has to be 
#' the control. The imputation step can be defined at certain level using the
#' treatment_level argument. Depending on the experimental design trait imputation
#' should occur a certain level, e.g. block or site.

#'
#' @return a tibble with extra class \code{imputed_trait}
#'
#' @importFrom stats sd var weighted.mean
#' @importFrom magrittr %>%
#' @importFrom dplyr select any_of all_of mutate group_by filter left_join n
#' @importFrom dplyr inner_join across ungroup
#' @importFrom purrr map_dfr
#' @importFrom rlang !!! !! .data
#' @importFrom glue glue glue_collapse
#' @importFrom tibble lst
#' @examples
#' data(community)
#' data(trait)
#' imputed_traits <-trait_impute(comm = community, traits = trait,
#'                  scale_hierarchy = c("Site", "PlotID"),
#'                  taxon_col = "Taxon", value_col = "Value",
#'                  trait_col = "Trait", abundance_col = "Cover")
#' @export

trait_impute <- function(
  comm,
  traits,
  scale_hierarchy = c("Country", "Site", "BlockID", "PlotID"),
  global = TRUE,
  taxon_col = "taxon",  trait_col = "trait",
  value_col = "Value", abundance_col = "Cover",
  treatment_col = NULL, treatment_level = NULL,
  other_col = character(0),
  keep_all = FALSE,
  min_n_in_sample = 5) {

  #### sanity checks on input (are columns present etc) ####
  #check data have all scales in scale_hierarchy
  if (!all(scale_hierarchy %in% names(comm))) {
    bad_scales <- glue_collapse(
      x = scale_hierarchy[!scale_hierarchy %in% names(comm)],
      sep = ", ", last = ", and ")
    stop(glue("scale_hierarchy levels {bad_scales} not in names(comm)"))
  }

  #check taxon_col is valid
  if (!all(taxon_col %in% names(comm))) {
    bad_taxon <- glue_collapse(
      x = taxon_col[!taxon_col %in% names(comm)],
      sep = ", ", last = ", and ")
    stop(glue("taxon_col {bad_taxon} not in names(comm)"))
  }

  if (!all(taxon_col %in% names(traits))) {
    bad_taxon <- glue_collapse(
      x = taxon_col[!taxon_col %in% names(traits)],
      sep = ", ", last = ", and ")
    stop(glue("taxon_col {bad_taxon} not in names(traits)"))
  }

  #check trait_col is valid
  if (!(length(trait_col) == 1 & trait_col %in% names(traits))) {
    stop(glue("trait_col '{trait_col}' not in names(traits)"))
  }

  #check value_col is valid
  if (!(length(value_col) == 1 & value_col %in% names(traits))) {
    stop(glue("value_col '{value_col}' not in names(traits)"))
  }

  #check abundance_col is valid
  if (!(length(abundance_col) == 1 & abundance_col %in% names(comm))) {
    stop(glue("value_col '{abundance_col}' not in names(comm)"))
  }

  #if used, check treatment_col is valid
  if (is.null(treatment_col)) {
    use_treat <- FALSE
  } else {
    use_treat <- TRUE
    if (!treatment_col %in% names(comm)) {
      stop(glue("treatment_col {treatment_col} not in names(comm)"))
    }
    if (!treatment_col %in% names(traits)) {
      stop(glue("treatment_col {treatment_col} not in names(traits)"))
    }
    if (!is.factor(comm[[treatment_col]])) {
      stop(glue("treatment_col {treatment_col} is not a factor in comm"))
    }
    if (!is.factor(traits[[treatment_col]])) {
      stop(glue("treatment_col {treatment_col} is not a factor in traits"))
    }
    if (!identical(
      levels(traits[[treatment_col]]),
      levels(comm[[treatment_col]]))
    ) {
      stop("treatment_col has have different levels in comm and traits")
    }
    #check treatment_level is valid
    if (is.null(treatment_level)) {
      stop("treatment_level must be specified when treatment_col is used")
    }
    if (!treatment_level %in% scale_hierarchy) {
      stop("treatment_level must be in scale_hierarchy")
    }
  }

  #check other_cols are valid 
  if (!all(other_col %in% names(comm))) {
    bad_other <- glue_collapse(
      x = other_col[!other_col %in% names(comm)],
      sep = ", ", last = ", and ")
    stop(glue("other_col levels {bad_other} not in names(comm)"))
  }
  #check other_cols are not in traits
  if (any(other_col %in% names(traits))) {
    bad_other <- glue_collapse(
      x = other_col[!other_col %in% names(traits)],
      sep = ", ", last = ", and ")
    warning(glue("other_col levels {bad_other} are in  names(traits). \\
                 These columns will be removed from traits"))
    traits <- traits %>% select(-any_of(other_col))
  }

  
  #### prep ####
  #add global to scale_hierarchy if necessary
  if (isTRUE(global)) {
    #check not already a "global" column
    if (any(names(traits) == "global")) {
      stop(glue("Cannot add global column as column called\\
                global already exists. Maybe set 'global = FALSE'"))
    }
    #add global column to traits
    traits <- traits %>% mutate(global = "global")
    comm <- comm %>% mutate(global = "global")
    scale_hierarchy <- c("global", scale_hierarchy)
  }

  #remove NA trait values
  traits <- traits %>% filter(!is.na(.data[[value_col]]))

  ## check for NA in abundance
  if (any(is.na(comm[[abundance_col]]))) {
    stop(glue("cannot have NA in the {abundance_col}\\
              column of the community data"))
  }

  ##### routine if length(taxon_col) > 1
  if (length(taxon_col) > 1) {
    result <- trait_impute_multi_level(
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

  #calculate plot scale sum of abundances
  comm <- comm %>%
    group_by(across(all_of(c(as.character(scale_hierarchy), other_col)))) %>%
    #calculate sum abundance
    mutate(sum_abun = sum(.data[[abundance_col]]))

  # make ordered factor of scale hierarchy
  scale_hierarchy <- factor(scale_hierarchy,
                            levels = rev(scale_hierarchy),
                            ordered = TRUE)

  ####iterate over grouping hierarchy####
  out <- scale_hierarchy %>%
    map_dfr(~{
      scale_level <- .x

      #drop scales from the hierarchy
      scale_drop <- as.character(
        scale_hierarchy[scale_hierarchy < scale_level])

      #scales to keep
      scale_keep <- as.character(
        scale_hierarchy[scale_hierarchy >= scale_level])

      traits <- traits %>% select(-any_of(scale_drop))
      result <- comm %>%
       ungroup() %>% #avoids dplyr problem #5473
       #join to traits
       inner_join(traits,
                  by = c(scale_keep, taxon_col),
                  suffix = c("_comm", "_trait")) %>%
       #group by kept scales
       group_by(
         across(all_of(c(scale_keep, other_col,
                         trait_col, taxon_col, "sum_abun")))
       )
      #filter if using treatment_col
      if (use_treat) {
        if (treatment_level == scale_level) {
          col_comm <- paste0(treatment_col, "_comm")
          col_trait <- paste0(treatment_col, "_trait")
          result <- result %>%
            filter(# same treatment OR first level (must be control)
              .data[[col_comm]] == .data[[col_trait]] |
              .data[[col_trait]] == levels(.data[[col_trait]])[1]
              )
        }
      }

      #calculate weights
      result <- result %>%
        ungroup() %>%
        group_by(across(all_of(c(as.character(scale_hierarchy),
                                 taxon_col, trait_col, other_col)))) %>%
       mutate(
         n_sample = n(),
         weight = .data[[abundance_col]] / .data$n_sample,
         level = scale_level
       )

      result
  }) # end of iterate over grouping hierarchy

  #check some output
  if (nrow(out) == 0) {
    stop("No output. Check taxon names and plot/site names in the community and trait data frames match.")
  }
  
  #get max number available
  out <- out %>%
    ungroup() %>%
    group_by(across(all_of(c(as.character(scale_hierarchy),
                             taxon_col, trait_col, other_col)))) %>%
    mutate(max_n_in_sample = max(.data$n_sample))

  #### filter out lowest good level of hierarchy for each taxon & trait ####
  if (!keep_all) {#keep only finest scale trait data available
    out <- out %>%
      ungroup() %>%
      group_by(across(all_of(c(as.character(scale_hierarchy),
                               taxon_col, trait_col, other_col)))) %>%
      filter(
        #group has more than minimum number in sample
        .data$n_sample >= min_n_in_sample |
          # OR maximum number when max is lower than required
       (.data$n_sample == .data$max_n_in_sample &
          .data$max_n_in_sample <= min_n_in_sample)
      ) %>%
      #filter lowest level available by group
      filter(.data$level == min(.data$level))
  }

  out <- out %>%
    ungroup() %>% #avoids dplyr problem #5473
    group_by(across(all_of(c(as.character(scale_hierarchy),
                             trait_col, other_col))))

  #add treatment group if used
  if (!is.null(treatment_col)) {
    out <- out %>%
      group_by(.data[[paste0(treatment_col, "_comm")]], .add = TRUE)
  }

  #set arguments as attributes so next functions have access to them
  attrib <- tibble::lst(scale_hierarchy, taxon_col, trait_col,
                        value_col, abundance_col, other_col, 
                        treatment_col, treatment_level)
  attr(out, "attrib") <- attrib

  class(out) <- c("imputed_trait", class(out))

  return(out)
}
