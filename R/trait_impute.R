#' Impute traits
#' @description Impute traits 
#' 
trait_impute <- function(comm, traits, scale_hierarchy = c("Country", "Site", "BlockID", "PlotID"), taxon_col = "taxon",  trait_col = "trait", value_col = "Value", abundance_col = "Cover", other_col = "Year"){
  
  out <- ((length(scale_hierarchy)):1) %>%  
    map_df({
      scale_drop <- character(0L)
      if(. < length(scale_hierarchy)){
        scale_drop <- scale_hierarchy[.:length(scale_hierarchy)]
      }
      scale_keep <- scale_hierarchy[1:.]

     trait %>% select(-one_of(scale_drop))
     comm %>% 
        left_join(trait, by = c(scale_keep, taxon)) %>% 
        group_by(!!!keep_scale, !!!taxon_col, !!!trait_col, !!!other_col) %>% 
      mutate(weight = !!!cover_col/n()) %>% 
#      group_by(!!!keep_scale, !!!trait_col, !!!other_col) %>% 
      mutate(level = factor(scale_hierarchy[.], levels = scale_hierarchy, ordered = TRUE))
  }) %>% 
    filter(!is.na(!!!value_col)) %>% 
    filter(level == min(level)) %>% 
    group_by(!!!keep_scale, !!!trait_col, !!!other_col)

  class(out) <- "imputed_traits"
  }




  
trait_bootstrap <- function(traits_comm){  
  stopifnot(class(traits_com) == "imputed_traits")
  
  BootstrapMoments_All <- rerun(.n = nrep, sample_n(traits_comm, size = samplesize,  replace = TRUE, weight = traits_comm$weight)) %>%
    bind_rows(.id = "n") %>% 
    group_by(n, add = TRUE) %>% 
    # get all the happy moments
    summarise(Mean = mean(Value), Variance = var(Value), Skewness = skewness(Value), Kurtosis = kurtosis(Value))
  
  return(BootstrapMoments_All)
}
