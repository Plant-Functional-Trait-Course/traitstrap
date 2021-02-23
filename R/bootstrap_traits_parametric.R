
#'Bootstrap traits parametrically
#' 
#'trait_parametric_bootstrap uses parametric bootstrapping to infer community trait moments.
#' @description Parametric bootstrapping of imputed traits
#' @param imputed_traits imputed trait and community data in long format
#' @param nrep number of bootstrap replicates
#' @param sample_size number of points to sample for each distribution
#' @param distribution_type the type of statistical distribution to use. Character or named list. Currently accepts "normal","lognormal", and "beta".  
#' @notes Distribution type should be a single character (e.g. "normal") or a named list (e.g. list(height="normal",mass="lognormal"))
#' @return a tibble
#' @importFrom stats var
#' @importFrom e1071 skewness kurtosis
#' @importFrom magrittr %>%
#' @importFrom dplyr slice_sample group_by summarise
#' @importFrom purrr map_df
#' @importFrom fitdistrplus fitdist
#' @examples {
#' data("community")
#' data("trait")
#' 
#' trait$Trait <- log10(trait$trait)
#' imputed_traits <- trait_impute(comm = community, traits = trait,scale_hierarchy = c("Site","PlotID"),
#' taxon_col = "Taxon",global = T,trait_col = "Trait",value_col = "Value")
#' 
#' parametric_results <- trait_parametric_bootstrap(imputed_traits = imputed_traits, distribution_type = "normal") 
#' parametric_summary <- trait_summarise_boot_moments(bootstrap_moments = parametric_results)
#' 
#' } 
#' @export
trait_parametric_bootstrap <- function(imputed_traits, nrep = 100, sample_size = 200, distribution_type = "normal"){
  
  
  #Pull useful information from imputed traits object
    value_col <- attributes(imputed_traits)$attrib$value_col
    trait_col <- attributes(imputed_traits)$attrib$trait_col
    taxon_col <- attributes(imputed_traits)$attrib$taxon_col
    abundance_col <- attributes(imputed_traits)$attrib$abundance_col
    scale_hierarchy <- attributes(imputed_traits)$attrib$scale_hierarchy
    
  
  #If only a single type of distribution was supplied, apply it to all traits
    if(length(distribution_type)==1){
      distributions <- unique(as.data.frame(imputed_traits)[,trait_col])
      dist_list <- replicate(n = length(distributions),expr = distribution_type ) 
      names(dist_list) <- distributions
      distribution_type <- as.list(dist_list)
      rm(dist_list,distributions)
    }#if statement
  
#######################################################################################################################################

#Data checks
    
    #Distribution type
    if(any(!unique(distribution_type)%in%c("beta","normal","lognormal"))){
      
      stop(paste("Unsupported distribution type:", setdiff(y = c("beta","normal","lognormal"),x = unique(distribution_type))))
    }
    
    #Beta checks
    if("beta" %in% distribution_type){
     
      
      #If beta distribution is specified, make sure every species x hierarchy combination has at least 2 data points, and that the values make sense
      
          
      imputed_traits[which(imputed_traits[[trait_col]] %in% names(distribution_type)[distribution_type=="beta"]),] %>% 
        group_by_at(c(as.character(scale_hierarchy),trait_col,taxon_col)) %>% summarise(n=n(),.groups="keep") -> beta_counts

      if(any(beta_counts$n>=1)){stop("Fitting a beta distrbution requires at least 2 points per distribution. We suggest re-imputing traits with a minimum of (at least) 2 leaves.")}
      
    
      #check that 0 <= values <= 1
      
      if(any(imputed_traits[which(imputed_traits[,trait_col]==names(distribution_type)[distribution_type=="beta"]),value_col]>1 |
      imputed_traits[which(imputed_traits[,trait_col]==names(distribution_type)[distribution_type=="beta"]),value_col]<0 )){
        stop("For a Beta distribution values must be between 0 and 1.")
      }
    
      
    }    
    
    #lognormal checks
    if("lognormal" %in% distribution_type){
      
      
      #check that 0 <= values <= 1
      
      if(any( imputed_traits[which(imputed_traits[,trait_col]==names(distribution_type)[distribution_type=="lognormal"]),value_col]<0 )){
        stop("For a Log-Normal distribution values cannot be negative.")
      }
      
      
    }    
    
    
        
#######################################################################################################################################    
    
    #First, we need to fit distributions For each trait x species x level of the hierarchy
    
    #Create table to hold the distribution parameters
    distributions_to_fit <- unique(imputed_traits[c(as.character(scale_hierarchy),taxon_col,trait_col)])
    distribution_parameters <- as.data.frame(distributions_to_fit)
    
    
    #Add columns for different distribution parameters as needed
    
    #Normal
    if("normal" %in% names(distribution_type)){
      
      distribution_parameters$norm_mean<-NA
      distribution_parameters$norm_sd<-NA
      
    }
    
    #lognormal
    
    if("lognormal" %in% names(distribution_type) | "lognormal" %in% names(distribution_type)){
      
      distribution_parameters$lnorm_meanlog <- NA
      distribution_parameters$lnorm_sdlog <- NA
      
    }
    
    #beta
    
    if("beta" %in% names(distribution_type)){
      
      distribution_parameters$beta_shape1 <- NA
      distribution_parameters$beta_shape2 <- NA
      
    }
    
    
    
#####################################################################################################################################
    
    
    for(i in 1:nrow(distributions_to_fit)){
      
      #pull all data matching the criteria  
      data_i <- imputed_traits[eval(parse(text = paste("which(",paste(paste("imputed_traits$",names(distributions_to_fit)," == '",distributions_to_fit[i,],"'",sep = ""  ),collapse = " & " ),")"))),value_col]
      data_i <- na.omit(data_i)  
      
      if(nrow(data_i)==0){next}
      
      #fit the specified distribution
      
      #Figure out what type of distribution to fit
      dist_type_i <- distribution_type[as.character(distributions_to_fit[i,trait_col])]
      
      
      
      
      #Normal
      if(dist_type_i=="normal"){
        
        if(nrow(data_i)==1){
          
          distribution_parameters$norm_mean[i] <- data_i
          distribution_parameters$norm_sd[i] <- 0  
          
        }else{ 
          
          fit_i <- NULL
          try(fit_i <- fitdistrplus::fitdist(data = unlist(data_i),distr = "norm",method = "mme"),silent = T)
          if(is.null(fit_i)){fit_i <- fitdistrplus::fitdist(data = unlist(data_i),distr = "norm",method = "mle")}
          
          
          distribution_parameters$norm_mean[i] <- fit_i$estimate[1]
          distribution_parameters$norm_sd[i] <- fit_i$estimate[2]
        }
      }#normal fit
      
      #log-normal
      
      if(dist_type_i %in% c("lognormal","log-normal")){
        
        if(nrow(data_i)==1){
          
          distribution_parameters$lnorm_meanlog[i] <- log(data_i)
          distribution_parameters$lnorm_sdlog[i] <- 0  
          
        }else{ 
          
          
          fit_i<-NULL
          
          try(fit_i <- fitdistrplus::fitdist(data = unlist(data_i),distr = "lnorm",method = "mme"),silent = T)
          
          if(is.null(fit_i)){fit_i <- fitdistrplus::fitdist(data = unlist(data_i),distr = "lnorm",method = "mle")}
          
          
          distribution_parameters$lnorm_meanlog[i] <- fit_i$estimate[1]
          distribution_parameters$lnorm_sdlog[i] <- fit_i$estimate[2]
        }
        
      }#lognormal
      
      
      ###########
      #beta
      
      if(dist_type_i %in% c("beta") ){
        
          fit_i <- NULL
          
          try(fit_i <- fitdistrplus::fitdist(data = unlist(data_i),distr = "beta",method = "mme"),silent = T)
          
          if(is.null(fit_i)){fit_i <- fitdistrplus::fitdist(data = unlist(data_i),distr = "beta",method = "mle")}
          
          
          distribution_parameters$beta_shape1[i] <- fit_i$estimate[1]
          distribution_parameters$beta_shape2[i] <- fit_i$estimate[2]

      }#beta

      
      ###########
      
      
      
      
      
      
      
      
      
      
      
            
      
      
    }#i loop
    
#######################################################################################################################################    
  
    #Now, sample the distributions
        
    #Create output data.frame
      distribution_moments <- as.data.frame(unique(imputed_traits[c(as.character(scale_hierarchy),trait_col)]))
      
      distribution_moments$mean <- NA
      distribution_moments$variance <- NA
      distribution_moments$skewness <- NA
      distribution_moments$kurtosis <- NA
      
      
      #Duplicate output to hold all replicates
      distribution_moments <- distribution_moments[rep(row.names(distribution_moments),nrep),]
      distribution_moments$n <- sort(rep(x =1:nrep,((nrow(distribution_moments)/nrep))))
      rownames(distribution_moments) <- NULL
             
      #distribution_moments[which(colnames(distribution_moments)!="n")]
      
      # Iterate through each scale hierarchy x trait combo and sample nreps
  
      for( i in 1:nrow(unique(distribution_moments[which(colnames(distribution_moments)!="n")]))){
        
        # Get the focal distrbution we'll be looking at
        distr_i <- unique(distribution_moments)[i,which(colnames(distribution_moments)!="n")]  
        
        # Get the necessary abundance info  
        abundances_i <- unique(imputed_traits[eval(parse(text = paste("which(",paste(paste("imputed_traits$",as.character(scale_hierarchy)," == '",distr_i[as.character(scale_hierarchy)],"'",sep = ""  ),collapse = " & " ),")"))),c(taxon_col,abundance_col)])
        
        # draw the traits
        
        out_i <- replicate(n = nrep, 
                           expr = draw_distribution(abundances = abundances_i, 
                                                    trait = distr_i[[trait_col]], 
                                                    attribs =  attributes(imputed_traits)$attrib, 
                                                    samples_per_dist =  sample_size,
                                                    distribution_type_list = distribution_type,
                                                    distribution_parameters_table = distribution_parameters, 
                                                    distr = distr_i))
        
        
        distribution_moments[eval(parse(text = paste("which( ",paste(paste("distribution_moments$",as.character(scale_hierarchy)," == '",distr_i[as.character(scale_hierarchy)],"'",
                                                                           sep = ""  ),collapse = " & " )," & distribution_moments$",trait_col, "== '",distr_i[trait_col],"' )",sep = "")) ),
                             c("mean","variance","skewness","kurtosis")] <-
        

        data.frame(mean = apply(X = out_i,MARGIN = 2,FUN = mean),
                   var = apply(X = out_i,MARGIN = 2,FUN = var),
                   skew = apply(X = out_i,MARGIN = 2,FUN = moments::skewness),
                   kurt = apply(X = out_i,MARGIN = 2,FUN = moments::kurtosis))
        
        
      
      }#for I loop
      
      
      #Make output comparable to nonparametric output
      distribution_moments <- distribution_moments[c(as.character(scale_hierarchy),trait_col,"n","mean","variance","skewness","kurtosis")]
      
      
      attributes(distribution_moments)$attrib <- attributes(imputed_traits)$attrib
      
      
      
      return(distribution_moments)
      
  
}#end function



######################

#' Draw traits from a parametric distribution
#' @description Draws traits
#' @param abundances
#' @param trait
#' @param attribs
#' @param samples_per_dist
#' @param distrbution_type_list
#' @param distribution_parametrics_table
#' @param distr
#' @description
#'
#' @return a vector of traits drawn from the specified distributions
#'
#' @importFrom stats rnorm rlnorm rbeta
#' @keywords internal
draw_distribution <- function(abundances, 
                              trait, 
                              attribs, 
                              samples_per_dist, 
                              distribution_type_list,
                              distribution_parameters_table,
                              distr){
  
  species_draw <- sample(x = abundances[[attribs$taxon_col]],prob = as.data.frame(abundances)[[attribs$abundance_col]],size = samples_per_dist, replace = T)
  tab <- as.data.frame(table(species_draw))

  ###########
    #Normal

  if(as.character(distribution_type_list[trait])=="normal"){

      samp <- apply(X = tab,MARGIN = 1,FUN =   function(x){
    
    
      #sample from the corresponding distribution
      distr_t  <- distribution_parameters_table[eval(parse(text = paste("which( ",paste(paste("distribution_parameters_table$",as.character(attribs$scale_hierarchy)," == '",distr[as.character(attribs$scale_hierarchy)],"'",
                                                                                        sep = ""  ),collapse = " & " ),"& distribution_parameters_table$",attribs$taxon_col," == '",
                                                                  as.character(x[[1]]),"' & distribution_parameters_table$",attribs$trait_col, "== '",distr[attribs$trait_col],"' )",sep = ""))),]  
      
      return(rnorm(n = as.numeric(x[[2]]),mean = as.numeric(distr_t$norm_mean[1]), sd = as.numeric(distr_t$norm_sd[1] )))
      
      
    } #apply fx 
  )
  
  samp <- unlist(samp)
  return(samp)  
    
  }#if normal
  
  ###########
    #lognormal
  if(as.character(distribution_type_list[trait])=="lognormal"){
    
    samp <- apply(X = tab,MARGIN = 1,FUN =   function(x){
      
      #sample from the corresponding distribution
      distr_t  <- distribution_parameters_table[eval(parse(text = paste("which( ",paste(paste("distribution_parameters_table$",as.character(attribs$scale_hierarchy)," == '",distr[as.character(attribs$scale_hierarchy)],"'",
                                                                                        sep = ""  ),collapse = " & " ),"& distribution_parameters_table$",attribs$taxon_col," == '",
                                                                  as.character(x[[1]]),"' & distribution_parameters_table$",attribs$trait_col, "== '",distr[attribs$trait_col],"' )",sep = ""))),]  
      
      return(rlnorm(n = as.numeric(x[[2]]),
                    meanlog =  as.numeric(distr_t$lnorm_meanlog[1]), 
                    sdlog = as.numeric(distr_t$lnorm_sdlog[1]) ))
      
    } #apply fx 
    )
    
    samp <- unlist(samp)
    return(samp)  
    
    
    
    }#if lognormal  
  

  ###########
    #beta
  
  if(as.character(distribution_type_list[trait])=="beta"){
    
    samp <- apply(X = tab,MARGIN = 1,FUN =   function(x){
          
      #sample from the corresponding distribution
      distr_t  <- distribution_parameters_table[eval(parse(text = paste("which( ",paste(paste("distribution_parameters_table$",as.character(attribs$scale_hierarchy)," == '",distr_t[as.character(attribs$scale_hierarchy)],"'",
                                                                                        sep = ""  ),collapse = " & " ),"& distribution_parameters_table$",attribs$taxon_col," == '",
                                                                  as.character(x[[1]]),"' & distribution_parameters_table$",attribs$trait_col, "== '",distr_t[trait_col],"' )",sep = ""))),]  
      
      return(rbeta(n = as.numeric(x[[2]]),
                   shape1 =  as.numeric(distr_t$beta_shape1[[1]] ), 
                   shape2 = as.numeric(distr_t$beta_shape2[[1]] ) ))
      
      
    } #apply fx 
    )
    
    samp <- unlist(samp)
    return(samp)  
    
  }#if beta
  
  ###########
  
    
    
}#fx

######################################################################


