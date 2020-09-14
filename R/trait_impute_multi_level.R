#' @importFrom purrr set_names map
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr filter mutate rename_with select bind_rows matches
#' @importFrom stringr str_remove
#'  
trait_impute_multi_level <- function(call){
  # Find position of taxon_col
  taxon_position <- match("taxon_col", names(call))
  
  #extract taxon_col value
  taxon_hierarchy <- call[[taxon_position]] 
  taxon_hierarchy <- eval(taxon_hierarchy) #make results of call into vector 
  taxon_hierarchy <- set_names(taxon_hierarchy)

  #iterate over taxon_col
  result <- map(taxon_hierarchy, ~{
    # change taxon_col in call and eval 
    call[[taxon_position]] <- as.character(.x)
    eval(call)
  }) %>% 
    map(select, -matches(paste0(taxon_hierarchy, "_trait"))) %>% 
    map(rename_with, str_remove, matches(paste0(taxon_hierarchy, "_comm")), "_comm") %>% 
    bind_rows(.id = "taxon_level")
   
  
  #filter lowest value of taxon hierarchy 
  # need to add species to group
  result <- result %>% 
  #  group_by(.data[[taxon_hierarchy]], .add = TRUE) %>% 
    mutate(taxon_level = factor(.data$taxon_level, levels = taxon_hierarchy, ordered = TRUE)) %>% 
    filter(.data$taxon_level == min(.data$taxon_level))
  
  # add taxon_hierarchy to attributes
  attr(result, "taxon_hierarchy") <- taxon_hierarchy
  
  #return
  return(result)
}