#' Impute traits
#' @description Impute traits 
#' @param comm community data in long format
#' @param traits trait data in long format
#' @param scale_hierarchy character vector of site/block/plot hierarchy (large to small)
#' @param taxon_col character; name of taxon column in comm and traits
#' @param trait_col character; name of trait name column in traits
#' @param value_col character; name of trait value column in traits
#' @param abundance_col character; name of species abundance column in comm
#' @param other_col name of other grouping columns in comm
#' @param treatment_col optional name of treatment_col in comm. Must refer to a factor where first level is control.
#' @param global logical; calculate traits at global scale. Must not be a column called global in the traits data.
#' @param keep_all logical; keep trait data at all available levels or just finest scale available
#' 
#' @description 
#' 
#' @return a tibble with extra class \code{imputed_trait}
#' 
#' @importFrom stats sd var weighted.mean
#' @importFrom magrittr %>%
#' @importFrom dplyr select any_of mutate group_by filter left_join n  rename inner_join across
#' @importFrom purrr map_df
#' @importFrom rlang !!! !! .data
#' @importFrom glue glue glue_collapse
#' @importFrom tibble lst
#' @export
 
trait_impute <- function(
  comm, 
  traits, 
  scale_hierarchy = c("Country", "Site", "BlockID", "PlotID"),
  global = TRUE,
  taxon_col = "taxon",  trait_col = "trait", 
  value_col = "Value", abundance_col = "Cover", 
  treatment_col,
  other_col = character(0), 
  keep_all = FALSE){
  
  #check data have all scales in scale_hierarchy
  if(!all(scale_hierarchy %in% names(comm))){
    bad_scales <- glue_collapse(
      x = scale_hierarchy[!scale_hierarchy %in% names(comm)], 
      sep = ", ", last = ", and ")
    stop(glue("scale_hierarchy levels {bad_scales} not in names(comm)"))
  }
  
  #if used, check treatment_col is valid
  if(!missing(treatment_col)){
    if(!treatment_col %in% names(comm)){
      stop(glue("treatment_col {treatment_col} not in names(comm)"))
    }
    if(!is.factor(comm[[treatment_col]])){
      stop(glue("treatment_col {treatment_col} is not a factor in comm"))
    }
  }
  
  
  #add global to scale_hierachy if necessary
  if(isTRUE(global)){
    #check not already a "global" column
    if(any(names(traits) == "global")){
      stop(glue("Cannot add global column as column called\\
                global already exists. Maybe set 'global = FALSE'"))
    }
    #add global column to traits
    traits <- traits %>% mutate(global = "global")
    comm <- comm %>% mutate(global = "global")
    scale_hierarchy <- c("global", scale_hierarchy)
  }
    
  #remove NA trait values
  traits <- traits %>% filter(!is.na(!!value_col))
  
  #calculate plot scale sum of abundances
  comm <- comm %>% 
    group_by(across(any_of(c(scale_hierarchy, other_col)))) %>%
    #calculate sum abundance
    rename(abundance = !!abundance_col) %>% 
    mutate(sum_abun = sum(.data$abundance)) 
  
  # make ordered factor of scale hierarchy
  scale_hierarchy <- factor(scale_hierarchy, 
                            levels = rev(scale_hierarchy), 
                            ordered = TRUE)
  
  #iterate over grouping hierarchy
  out <- scale_hierarchy %>%  
    map_df(~{#browser()
      scale_level <- .x #catch the . -meaning changes within the dplyr chain
      
      #drop scales from the hierarchy
      scale_drop <- character(0L)
      if(scale_level != max(scale_hierarchy)){
        scale_drop <- scale_hierarchy[(scale_level + 1):length(scale_hierarchy)]
      }
      scale_keep <- character(0L)
      if(scale_level > 0){
        scale_keep <- scale_hierarchy[1:scale_level]
      }
     traits <- traits %>% select(-any_of(scale_drop))
     comm %>%
       #group by kept scales
       group_by(across(any_of(c(scale_keep, other_col)))) %>%
       #join to traits 
       inner_join(traits, by = c(scale_keep, taxon_col)) %>% 
       group_by(
         across(any_of(c(trait_col, taxon_col, "sum_abun"))), 
         .add = TRUE
       ) %>%
       #calculate weights
       mutate(
         weight = .data$abundance/n(),
         level = factor(scale_hierarchy[scale_level], 
                        levels = scale_hierarchy, ordered = TRUE)
         )
  }) %>% 
    filter(!is.na(!!!value_col))#remove NA values


  if(!keep_all){#keep only finest scale trait data available
    out <- out %>% 
      filter(.data$level == max(.data$level)) 
  }
  
  out <- out %>%     
    group_by(across(any_of(c(scale_hierarchy, trait_col, other_col))))
  
  
  #set arguments as attributes so next functions have access to them
  attrib <- tibble::lst(scale_hierarchy, taxon_col, trait_col, 
                        value_col, abundance_col, other_col)
  attr(out, "attrib") <- attrib

  class(out) <- c("imputed_trait", class(out))
  
  out
  }











