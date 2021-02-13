#' Bootstrap traits parametrically
#' @description Bootstrap imputed traits using parametric bootstrapping 
#' @param imputed_traits imputed trait and community data in long format
#' @param distribution_type the type of statistical distribution to use. Character or named list. Currently accepts "normal","lognormal", and "beta".  
#' @param nrep number of bootstrap replicates
#' @param samples_per_abundance number trait samples to draw per unit of abundance
#' @notes To use multiple types of distributions, supply a named list: list( normal = "leaf_area_mm2", 
#' lognormal = "dry_mass_mg", normal = "SLA_m2_kg", lognormal = "height", normal = "biomass_per_ind")
#' @notes No, I don't remember why I specified the list this way rather than more intuitively.  Meh.
#' 
#' @return 
#' 
#' @importFrom stats var 
#' @importFrom e1071 skewness kurtosis
#' @export
trait_parametric_bootstrap <- function( imputed_traits, 
                                       distribution_type, 
                                       nrep = 100, 
                                       samples_per_abundance = 1 ){
  

  value_col <- attributes(imputed_traits)$attrib$value_col
  trait_col <- attributes(imputed_traits)$attrib$trait_col
  taxon_col <- attributes(imputed_traits)$attrib$taxon_col
  scale_hierarchy <- attributes(imputed_traits)$attrib$scale_hierarchy
  
  
  #If only a single type of distribution was supplied, apply it to all traits
    if(length(distribution_type)==1){
    
      #distributions <- unique(imputed_traits$Trait)
      distributions <- unique(as.data.frame(imputed_traits)[,trait_col])
      
      names(distributions)[1:length(distributions)] <- distribution_type
      distribution_type<-distributions
      rm(distributions)
    }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  #First, we need to fit distributions
    #For each trait x species x level of the hierarchy (note this only works if hierarchy is only geographic)
  
  distributions_to_fit <- unique(imputed_traits[c(as.character(scale_hierarchy),taxon_col,trait_col)])
  
  distribution_parameters <- as.data.frame(distributions_to_fit)

  
  #Add columns for different distribution parameters as needed
  
  #Normal
    if("normal" %in% names(distribution_type)){
      
      distribution_parameters$norm_mean<-NA
      distribution_parameters$norm_sd<-NA
      
    }
  
  #lognormal
  
    if("lognormal" %in% names(distribution_type) | "log-normal" %in% names(distribution_type)){
      
      distribution_parameters$lnorm_meanlog <- NA
      distribution_parameters$lnorm_sdlog <- NA
      
    }
    
  #beta
  
  if("beta" %in% names(distribution_type)){
    
    distribution_parameters$beta_shape1 <- NA
    distribution_parameters$beta_shape2 <- NA
    
  }
  
  
  
  
  for(i in 1:nrow(distributions_to_fit)){
  
    
    #pull all data matching the criteria  
      data_i <- imputed_traits[eval(parse(text = paste("which(",paste(paste("imputed_traits$",names(distributions_to_fit)," == '",distributions_to_fit[i,],"'",sep = ""  ),collapse = " & " ),")"))),value_col]
      data_i <- na.omit(data_i)  
      
      if(nrow(data_i)==0){next}
      
    #fit the specified distribution
      
      #Figure out what type of distribution to fit
      dist_type_i <- names(distribution_type)[which(distribution_type==as.character(distributions_to_fit[i,trait_col]))]
    
      #Normal
      if(dist_type_i=="normal"){
      
        if(nrow(data_i)==1){
        
          distribution_parameters$norm_mean[i] <- data_i
          distribution_parameters$norm_sd[i] <- 0  
            
        }else{ 
      
              fit_i<-NULL
              
              try(fit_i <- fitdistrplus::fitdist(data = unlist(data_i),distr = "norm",method = "mle"),silent = T)
              
              if(is.null(fit_i)){fit_i <- fitdistrplus::fitdist(data = unlist(data_i),distr = "norm",method = "mme")}
              
              
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
          
          try(fit_i <- fitdistrplus::fitdist(data = unlist(data_i),distr = "lnorm",method = "mle"),silent = T)
          
          if(is.null(fit_i)){fit_i <- fitdistrplus::fitdist(data = unlist(data_i),distr = "lnorm",method = "mme")}
    
          
        distribution_parameters$lnorm_meanlog[i] <- fit_i$estimate[1]
        distribution_parameters$lnorm_sdlog[i] <- fit_i$estimate[2]
        }
        
      }#lognormal
      
      #beta
      
      if(dist_type_i %in% c("beta") ){
        
        if(nrow(data_i)==1){
          
          if(!"unif_min" %in% colnames(distribution_parameters)){
          distribution_parameters$unif_min <- NA
          distribution_parameters$unif_max <- NA
            
          }
          
          
          distribution_parameters$unif_min[i] <- data_i
          distribution_parameters$unif_max[i] <- data_i
          
          
        }else{ 
          
          
          fit_i <- NULL
          
          try(fit_i <- fitdistrplus::fitdist(data = unlist(data_i),distr = "beta",method = "mle"),silent = T)
          
          if(is.null(fit_i)){fit_i <- fitdistrplus::fitdist(data = unlist(data_i),distr = "beta",method = "mme")}
          
          
          distribution_parameters$beta_shape1[i] <- fit_i$estimate[1]
          distribution_parameters$beta_shape2[i] <- fit_i$estimate[2]
        }
        
      }#lognormal
      
    

  }#i loop
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  
  #Next, we need to use the fitted distributions to estimate moments
  
  distribution_moments <- as.data.frame(unique(imputed_traits[c(as.character(scale_hierarchy),trait_col)]))
  
  distribution_moments$mean <- NA
  distribution_moments$variance <- NA
  distribution_moments$skewness <- NA
  distribution_moments$kurtosis <- NA
  
  
  #Duplicate output to hold all replicates
  distribution_moments <- distribution_moments[rep(row.names(distribution_moments),nrep),]
  distribution_moments$n <- NA
  
  #Now, iterate through each scale hierarchy x trait combo and sample nreps
  
  
  for( i in 1:nrow(unique(distribution_moments))){
  #print(i)
  distr_i <- unique(distribution_moments)[i,]  
    
  #abundances_i <- unique(imputed_traits[eval(parse(text = paste("which(",paste(paste("imputed_traits$",scale_hierarchy," == '",distr_i[scale_hierarchy],"'",sep = ""  ),collapse = " & " ),")"))),c(taxon_col,"abundance")])
  
  abundances_i <- unique(imputed_traits[eval(parse(text = paste("which(",paste(paste("imputed_traits$",as.character(scale_hierarchy)," == '",distr_i[as.character(scale_hierarchy)],"'",sep = ""  ),collapse = " & " ),")"))),c(taxon_col,"abundance")])
  
  out_i <- matrix(nrow = nrep, ncol=sum(floor(abundances_i$abundance*samples_per_abundance)))
  
  
  #pull all data matching the criteria  
  
  for(t in 1:nrow(abundances_i)){
  
  
    

  #distr_t  <- distribution_parameters[eval(parse(text = paste("which( ",paste(paste("distribution_parameters$",scale_hierarchy," == '",distr_i[scale_hierarchy],"'",
  #                                                                        sep = ""  ),collapse = " & " ),"& distribution_parameters$",taxon_col," == '",
  #                                                  abundances_i[t,taxon_col],"' & distribution_parameters$",trait_col, "== '",distr_i[trait_col],"' )",sep = ""))),]  
  
  distr_t  <- distribution_parameters[eval(parse(text = paste("which( ",paste(paste("distribution_parameters$",as.character(scale_hierarchy)," == '",distr_i[as.character(scale_hierarchy)],"'",
                                                                                    sep = ""  ),collapse = " & " ),"& distribution_parameters$",taxon_col," == '",
                                                              abundances_i[t,taxon_col],"' & distribution_parameters$",trait_col, "== '",distr_i[trait_col],"' )",sep = ""))),]  
  
  
  

  #If the specified trait was not fit (missing data), supply NAs
  if(nrow(distr_t)==0){
    

    draw_t <- rep(NA,(floor(abundances_i$abundance[t]*samples_per_abundance)*nrep))
    
    #Record output (if any)
    
    if(length(draw_t)>0){
      
      
      #This bit fills the first empty columns in the output file
      out_i[,which(apply(X = out_i,MARGIN = 2,FUN = function(x){all(is.na(x))}))[1]:(which(apply(X = out_i,MARGIN = 2,FUN = function(x){all(is.na(x))}))[1]+(floor(abundances_i$abundance[t]*samples_per_abundance)-1))] <- draw_t
      
    }#if there is output
    
    #break out of iteration t with next, otherwise errors occur with following if statements
    next 
  }
  

  #If the specified trait has a normal distribution, draw from that
  if("norm_mean" %in% colnames(distr_t)){
  if(!is.na(distr_t$norm_mean)){
  
    draw_t <- rnorm(n = (floor(abundances_i$abundance[t]*samples_per_abundance)*nrep),
                    mean = distr_t$norm_mean[[1]],
                    sd =  distr_t$norm_sd[[1]])


        #Record output (if any)
    
        if(length(draw_t)>0){
          
          
        #This bit fills the first empty columns in the output file
        out_i[,which(apply(X = out_i,MARGIN = 2,FUN = function(x){all(is.na(x))}))[1]:(which(apply(X = out_i,MARGIN = 2,FUN = function(x){all(is.na(x))}))[1]+(floor(abundances_i$abundance[t]*samples_per_abundance)-1))] <- draw_t

        }#if there is output
  
    
  }}
  
  
  
  #Log normal fitting if needed

  if("lnorm_meanlog" %in% colnames(distr_t)){
  
    if(!is.na(distr_t$lnorm_meanlog)){
      
      draw_t <- rlnorm(n = (floor(abundances_i$abundance[t]*samples_per_abundance)*nrep),
                      meanlog = distr_t$lnorm_meanlog[[1]],
                      sdlog = distr_t$lnorm_sdlog[[1]])

      #Record output (if any)
      
      if(length(draw_t)>0){
        
        
        #This bit fills the first empty columns in the output file
        out_i[,which(apply(X = out_i,MARGIN = 2,FUN = function(x){all(is.na(x))}))[1]:(which(apply(X = out_i,MARGIN = 2,FUN = function(x){all(is.na(x))}))[1]+(floor(abundances_i$abundance[t]*samples_per_abundance)-1))] <- draw_t
        
      }#if there is output
      
      
      
    }  
    
      
  }
  
  #beta
  
  if("beta_shape1" %in% colnames(distr_t)){
    
    if(!is.na(distr_t$beta_shape1)){
      
      draw_t <- rbeta(n = (floor(abundances_i$abundance[t]*samples_per_abundance)*nrep),
                      shape1 = distr_t$beta_shape1[[1]],
                      shape2 = distr_t$beta_shape2[[1]]
                     )
      
      #Record output (if any)
      
      if(length(draw_t)>0){
        
        
        #This bit fills the first empty columns in the output file
        out_i[,which(apply(X = out_i,MARGIN = 2,FUN = function(x){all(is.na(x))}))[1]:(which(apply(X = out_i,MARGIN = 2,FUN = function(x){all(is.na(x))}))[1]+(floor(abundances_i$abundance[t]*samples_per_abundance)-1))] <- draw_t
        
      }#if there is output
      
      
      
    }  
    
    
  }
  
  
  #unif
  
  if("unif_min" %in% colnames(distr_t)){
    
    if(!is.na(distr_t$unif_min)){
      
      draw_t <- runif(n = (floor(abundances_i$abundance[t]*samples_per_abundance)*nrep),
                      min = distr_t$unif_min[[1]],
                      max = distr_t$unif_max[[1]]
      )
      
      #Record output (if any)
      
      if(length(draw_t)>0){
        
        
        #This bit fills the first empty columns in the output file
        out_i[,which(apply(X = out_i,MARGIN = 2,FUN = function(x){all(is.na(x))}))[1]:(which(apply(X = out_i,MARGIN = 2,FUN = function(x){all(is.na(x))}))[1]+(floor(abundances_i$abundance[t]*samples_per_abundance)-1))] <- draw_t
        
      }#if there is output
      
      
      
    }  
    
    
  }
  
    
  }#t loop, drawing traits for each species
  
  
  
  #Calculate the moments of each replicate and populate the distribution_moments object accordingly
  
  distribution_moments[eval(parse(text = paste("which( ",paste(paste("distribution_moments$",as.character(scale_hierarchy)," == '",distr_i[as.character(scale_hierarchy)],"'",
                                           sep = ""  ),collapse = " & " )," & distribution_moments$",trait_col, "== '",distr_i[trait_col],"' )",sep = "")) ),
                       c("mean","variance","skewness","kurtosis")] <-
  
  
    cbind(
      rowMeans(out_i),
      apply(X = out_i,MARGIN = 1,FUN = var),
      apply(X = out_i,MARGIN = 1,FUN = skewness),
      apply(X = out_i,MARGIN = 1,FUN = kurtosis)
    )
  
  
  #Assign rep n
  distribution_moments[eval(parse(text = paste("which( ",paste(paste("distribution_moments$",as.character(scale_hierarchy)," == '",distr_i[as.character(scale_hierarchy)],"'",
                                                                     sep = ""  ),collapse = " & " )," & distribution_moments$",trait_col, "== '",distr_i[trait_col],"' )",sep = "")) ),"n"] <- 1:nrep
    
    
    
  
    
  }#i loop sampling para distr.
  

  #Make output comparable to nonparametric output
  distribution_moments <- distribution_moments[c(as.character(scale_hierarchy),trait_col,"n","mean","variance","skewness","kurtosis")]
  
  
  attributes(distribution_moments)$attrib <- attributes(imputed_traits)$attrib
  
  
  
  return(distribution_moments)
  
  
}#function




