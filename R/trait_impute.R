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
#' @param global logical; calculate traits at global scale. Must not be a column called global in the traits data.
#' 
#' @description 
#' 
#' @return 
#' 
#' @importFrom stats sd var weighted.mean
#' @importFrom magrittr %>%
#' @importFrom dplyr select one_of mutate group_by filter left_join n group_by_at rename
#' @importFrom purrr map_df
#' @importFrom rlang !!! !! .data
#' @importFrom glue glue glue_collapse
#' @export
 
trait_impute <- function(comm, traits, 
                         scale_hierarchy = c("Country", "Site", "BlockID", "PlotID"),
                         global = TRUE,
                         taxon_col = "taxon",  trait_col = "trait", 
                         value_col = "Value", abundance_col = "Cover", 
                         other_col = character(0)){
  #check data have all scales in scale_hierarchy
  if(!all(scale_hierarchy %in% names(comm))){
    bad_scales <- glue_collapse(
      x = scale_hierarchy[!scale_hierarchy %in% names(comm)], 
      sep = ", ", last = ", and ")
    stop(glue("scale_heirarchy levels {bad_scales} not in names(comm)"))
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
    group_by_at(vars(one_of(c(scale_hierarchy, other_col)))) %>%
    #calculate sum abundance
    rename(abundance = !!abundance_col) %>% 
    mutate(sum_abun = sum(.data$abundance)) 
  
  
  
  #iterate over grouping hierarchy
  out <- length(scale_hierarchy):1 %>%  
    map_df(~{#browser()
      scale_level <- . #catch the . -meaning changes within the dplyr chain
      
      #drop scales from the hierarchy
      scale_drop <- character(0L)
      if(scale_level < length(scale_hierarchy)){
        scale_drop <- scale_hierarchy[(scale_level + 1):length(scale_hierarchy)]
      }
      scale_keep <- character(0L)
      if(scale_level > 0){
        scale_keep <- scale_hierarchy[1:scale_level]
      }
     traits <- traits %>% select(-one_of(scale_drop))
     comm %>%
       #group by kept scales
       group_by_at(vars(one_of(c(scale_keep, other_col)))) %>%
       #join to traits 
       inner_join(traits, by = c(scale_keep, taxon_col)) %>% 
       group_by_at(vars(one_of(c(trait_col, taxon_col, "sum_abun"))), .add = TRUE) %>%
       #calculate weights
       mutate(
         weight = .data$abundance/n(),
         level = factor(scale_hierarchy[scale_level], 
                        levels = scale_hierarchy, ordered = TRUE)
         )
  }) %>% 
    filter(!is.na(!!!value_col)) %>% #remove NA values
    filter(.data$level == max(.data$level)) %>% 
    group_by_at(.vars = vars(one_of(c(scale_hierarchy, trait_col, other_col))))

  #set arguments as attributes so next functions have access to them
  attr(out, "scale_hierarchy") <- scale_hierarchy
  attr(out, "taxon_col") <- taxon_col
  attr(out, "trait_col") <- trait_col
  attr(out, "value_col") <- value_col
  attr(out, "abundance_col") <- abundance_col
  attr(out, "other_col") <- other_col
  
  class(out) <- c(class(out), "imputed_traits")
  
  out
  }




#' Bootstrap traits
#' @description Bootstrap impute traits 
#' @param imputed_traits imputed trait and community data in long format
#' @param nrep number of bootstrap replicates
#' @param sample_size number of plants per sample
#' @description 
#' 
#' @return 
#' 
#' @importFrom stats var 
#' @importFrom e1071 skewness kurtosis
#' @importFrom magrittr %>%
#' @importFrom dplyr sample_n group_by summarise_at vars one_of
#' @importFrom purrr map_df
#' @export

trait_np_bootstrap <- function(imputed_traits, nrep = 100, sample_size = 200){  
#  stopifnot(class(traits_com) == "imputed_traits")
  value_col <- attr(imputed_traits, "value_col")
  bootstrapMoments <- map_df(
    1:nrep,
    ~sample_n(imputed_traits, size = sample_size,  replace = TRUE, weight = weight),
    .id = "n"
    ) %>% 
    group_by(n, add = TRUE) %>% 
    # get all the happy moments
    summarise_at(
      .vars = vars(one_of(value_col)), 
      .funs = list(mean = mean, variance = var, 
                   skewness = skewness, kurtosis = kurtosis))
  
  attr(bootstrapMoments, "scale_hierarchy") <- attr(imputed_traits, "scale_hierarchy")
  attr(bootstrapMoments, "trait_col") <- attr(imputed_traits, "trait_col") 
  attr(bootstrapMoments, "other_col") <- attr(imputed_traits, "other_col")
  return(bootstrapMoments)
}


#' Summarise Bootstrap traits
#' @description Bootstraptraits 
#' @param BootstrapMoments trait moments from trait_np_bootstrap

#' @description 
#' 
#' @return 
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr n group_by summarise_at vars one_of group_by_at summarise
#' @importFrom rlang .data
#' @export

SummariseBootMoments <- function(BootstrapMoments){
  groups <- c(attr(BootstrapMoments, "scale_hierarchy"), 
              attr(BootstrapMoments, "trait_col"),
              attr(BootstrapMoments, "other_col"))
              
  # calculate means of moments 
  sBootstrapMoments <- BootstrapMoments %>% 
    group_by_at(vars(one_of(groups))) %>% 
    summarise(
      n = n(),
      Mean = mean(.data$mean),
      CIlow.mean = .data$Mean - sd(.data$mean),
      CIhigh.mean = .data$Mean + sd(.data$mean),
      
      Var = mean(.data$variance),
      CIlow.var = .data$Var - sd(.data$variance),
      CIhigh.var = .data$Var + sd(.data$variance),
      
      Skew = mean(.data$skewness),
      CIlow.skew = .data$Skew - sd(.data$skewness),
      CIhigh.skew = .data$Skew + sd(.data$skewness),
      
      Kurt = mean(.data$kurtosis),
      CIlow.kurt = .data$Kurt - sd(.data$kurtosis),
      CIhigh.Kurt = .data$Kurt + sd(.data$kurtosis)
    )
  
  return(sBootstrapMoments)
}



#' Summarise Imputed Traits
#' @param imputed_traits trait moments from trait_np_bootstrap
#' @description Shows at which level the data are coming from in each plot.
#' 
#' @return a ggplot
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise distinct bind_cols ungroup one_of
#' @importFrom rlang .data
#' @importFrom ggplot2 autoplot ggplot geom_col facet_wrap aes scale_y_continuous labs
#' @export


autoplot.imputed_trait <- function(imputed_traits){
  #get scale_hierarchy and concatenate to make an ID
  scale_hierarchy <- attr(imputed_traits, which = "scale_hierarchy")
  scale_hierarchy <- scale_hierarchy[scale_hierarchy != "global"]
  id <- imputed_traits %>% 
    ungroup() %>% 
    select(one_of(scale_hierarchy)) %>% 
    apply(1, paste, collapse = "_")
  
  imputed_traits %>% 
    bind_cols(.id = id) %>% 
    group_by(.data$.id, .data$level, .data$Taxon, add = TRUE) %>% 
    distinct(.data$abundance, .keep_all = TRUE) %>%  
    summarise(s = sum(.data$abundance)/.data$sum_abun) %>% 
    summarise(s = sum(.data$s)) %>% 
    ggplot(aes(x = .data$.id, 
               y = .data$s, fill = .data$level)) + 
    geom_col() + 
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Plot", y = "Proportion of cover", fill = "Data source") + 
    facet_wrap(~.data$Trait)
}
