#' @importFrom purrr set_names map
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr filter mutate rename_with select bind_rows matches groups
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
  }) 
  #grab attributes
  result_attr <- attr(result[[1]], "attrib")

  result <- result  %>% 
    map(select, -matches(paste0(taxon_hierarchy, "_trait"))) %>% 
    map(rename_with, str_remove, matches(paste0(taxon_hierarchy, "_comm")), "_comm") %>% 
    bind_rows(.id = "taxon_level")
   
  
  #filter lowest value of taxon hierarchy 
  # need to add species to group
  current_groups <- groups(result) %>% as.character()
  result <- result %>% 
    ungroup() %>% 
    mutate(taxon_level = factor(.data$taxon_level, levels = taxon_hierarchy, ordered = TRUE)) %>% 
    group_by(across(all_of(c(current_groups, taxon_hierarchy[1])))) %>% 
    filter(.data$taxon_level == min(.data$taxon_level))   %>% 
    ungroup() %>% 
    group_by(across(all_of(current_groups))) 
  
  #reinstate attributes
  attr(result, "attrib") <- result_attr
  
  # add taxon_hierarchy to attributes
  attr(result, "taxon_hierarchy") <- taxon_hierarchy
  
  #add class
  class(result) <- c("imputed_trait", class(result))
  #return
  return(result)
}