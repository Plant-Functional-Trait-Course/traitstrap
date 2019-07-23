#trait parametric bootstrap



#for each species x site x trait: fit a distribution

#for each iteration, draw traits from distribution proportional to abundance



bs <- trait_np_bootstrap(imputed_traits = timp,
                       nrep = 100,
                       sample_size = 10)


distribution_typetest<-c(a="normal")
distribution_type<-"normal"


#' Bootstrap traits
#' @description Bootstrap imputed traits using parametric bootstrapping 
#' @param imputed_traits imputed trait and community data in long format
#' @param nrep number of bootstrap replicates
#' @param samples_per_abundance number trait samples to draw per unit of abundance
#' @param distribution_type the type of statistical distribution to use. Character or named list. Currently accepts "normal" or "lognormal"  
#' @description 
#' 
#' @return 
#' 
#' @importFrom stats var 
#' @importFrom e1071 skewness kurtosis
#' @export
trait_parametric_bootstrap(imputed_traits,
                           nrep=100,
                           samples_per_abundance=1,
                           distribution_type){
  

  value_col <- attr(imputed_traits, "value_col")
  trait_col <- attr(imputed_traits, "trait_col")
  taxon_col <- attr(imputed_traits, "taxon_col")
  scale_hierarchy <- attr(imputed_traits, "scale_hierarchy")
  
  
  #If only a single type of distribution was supplied, apply it to all traits
    if(length(distribution_type)==1){
    
      distributions <- unique(imputed_traits$Trait)
      names(distributions)[1:length(distributions)] <- distribution_type
      distribution_type<-distributions
      rm(distributions)
    }
      
  #First, we need to fit distributions
    #For each trait x species x level of the hierarchy (note this only works if hierarchy is only geographic)
  
  distributions_to_fit <- unique(imputed_traits[c(scale_hierarchy,taxon_col,trait_col)])
  
  distribution_parameters<-distributions_to_fit
  
  for(i in 1:nrow(distributions_to_fit)){
  
    
    #pull all data matching the criteria  
      data_i <- imputed_traits[eval(parse(text = paste("which(",paste(paste("imputed_traits$",names(distributions_to_fit)," == '",distributions_to_fit[i,],"'",sep = ""  ),collapse = " & " ),")"))),value_col]
    
      
    #fit the specified distribution
      
      #Figure out what type of distribution to fit
      dist_type_i <- names(distribution_type)[which(distribution_type==as.character(distributions_to_fit[i,trait_col]))]
    
      #Normal
      if(dist_type_i=="normal"){
      
      fit_i <- fitdistrplus::fitdist(data = unlist(data_i),distr = "norm",method = "mle")
        
        distribution_parameters$norm_mean[i] <- fit_i$estimate[1]
        distribution_parameters$norm_sd[i] <- fit_i$estimate[2]
      
      }
      
      #log-normal
      
      if(dist_type_i %in% c("lognormal","log-normal")){
        
        fit_i <- fitdistrplus::fitdist(data = unlist(data_i),distr = "lnorm")
        
        
        distribution_parameters$lnorm_meanlog[i] <- fit_i$estimate[1]
        distribution_parameters$lnorm_sdlog[i] <- fit_i$estimate[2]
        
      }
      
      
    
    
    
  }#i loop
  
  
  

}#function











