#' @importFrom purrr map_df
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr filter mutate
#'  
trait_impute_multi_level <- function(call){
  # Find position of taxon_col
  taxon_position <- match("taxon_col", names(call))
  
  #extract taxon_col value
  taxon_hierarchy <- call[[taxon_position]]

  #iterate over taxon_col
  result <- map_df(taxon_hierarchy, ~{
    # change taxon_col in call and eval 
    call[[taxon_position]] <- as.character(.x)
    eval(call)
  }, .id = "taxon_level")
   
  
  #filter lowest value of taxon hierarchy 
  #warning - depends on grouping structure being correct!
  result <- result %>% 
    mutate(taxon_level = factor(.data$taxon_level, levels = taxon_hierarchy, ordered = TRUE)) %>% 
    filter(.data$taxon_level == min(.data$taxon_level))
  
  # add taxon_hierarchy to attributes
  attr(result, "taxon_hierarchy") <- taxon_hierarchy
  
  #return
  return(result)
}