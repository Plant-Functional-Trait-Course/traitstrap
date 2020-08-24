### CALCULATE COMMUNITY WEIGHTED MEANS ###

# NORMAL WAY WITHOUT BOOTSTRAPPING
CommunityW_GlobalAndLocalMeans <- function(dat){
  meanTraits <- dat$trait %>% 
    # Global means
    group_by(Taxon, Trait) %>% 
    mutate(TraitMean_global = mean(Value, na.rm = TRUE)) %>% 
    
    # Site means (site level)
    group_by(Site, Taxon, Trait) %>%
    mutate(TraitMean_site = mean(Value, na.rm = TRUE)) %>% 
    
    # Plot means
    group_by(PlotID, BlockID, Site, Taxon, Trait) %>%
    mutate(TraitMean_plot = mean(Value, na.rm = TRUE)) %>% 
    select(-Year, -Value) %>% 
    ungroup() %>% 
    distinct()
  
  dat2 <- dat$community 
  
  if(!dat2$Country[1] %in% c("NO", "CO")) {
    dat2 <- dat2 %>% 
      # join plot level means
      left_join(meanTraits %>% select(-TraitMean_global, -TraitMean_site))
  }
  
  dat2 <- dat2 %>% 
    # join site level means
    left_join(meanTraits %>% select(-TraitMean_global, -TraitMean_plot, -BlockID, -PlotID) %>% distinct()) %>% 
    
    # join global means
    left_join(meanTraits %>% select(-TraitMean_plot, -TraitMean_site, -Site, -BlockID, -PlotID) %>% distinct())
  
  if(dat2$Country[1] %in% c("NO", "CO")) {
    dat2 <- dat2 %>% 
      mutate(TraitMean = coalesce(TraitMean_site, TraitMean_global))
  } else{
    dat2 <- dat2 %>% 
      mutate(TraitMean = coalesce(TraitMean_plot, TraitMean_site, TraitMean_global))
  }
  
  ### Calculate Community weighted means
  dat2 <- dat2 %>% 
    group_by(Trait, Site, PlotID) %>% 
    mutate(CWTraitMean = weighted.mean(TraitMean, Cover, na.rm=TRUE)) %>% 
    ungroup()
  
  return(dat2)
}

# Reduce to only CWMeans
CommunityW_Means <- function(TraitMeans_All){
  CWTraitMeans <- TraitMeans_All %>% 
    filter(!is.na(Trait)) %>% 
    select(-Taxon, -Cover, -TraitMean_plot, -TraitMean_site, -TraitMean_global, -TraitMean) %>% 
    distinct()
  return(CWTraitMeans)
}


# TRANSFORMING THE TRAITS
LogTranformation <- function(dat){
  dat$trait_trans <- dat$trait %>% 
    
    # log transform
    mutate(Value = ifelse(Trait %in% c("Plant_Height_cm", "Wet_Mass_g", "Dry_Mass_g", "Leaf_Area_cm2"), log(Value), Value),
           Trait = recode(Trait, "Plant_Height_cm" = "Plant_Height_cm_log", "Wet_Mass_g" = "Wet_Mass_g_log", "Dry_Mass_g" = "Dry_Mass_g_log", "Leaf_Area_cm2" = "Leaf_Area_cm2_log"))
  return(dat)
}


# WITH BOOTSTRAPPING
CWM_Bootstrapping <- function(dat, nrep = 100, samplesize = 200){
  comm <- dat$community %>% 
    filter(!Cover == 0) %>% 
    group_by(Country, Year, Site, Gradient, BlockID, PlotID) %>% 
    mutate(sumCover = sum(Cover))
  
  trait <- dat$trait_trans %>% filter(!is.na(Value))
  
  TraitWeights_plot <- comm %>% 
    left_join(trait %>% select(-Year), by = c("Country", "Site", "Gradient", "BlockID", "PlotID", "Taxon")) %>% 
    group_by(Country, Year, Site, Gradient, BlockID, PlotID, Taxon, Trait) %>% 
    mutate(weight = Cover/n()) %>% 
    group_by(Country, Year, Site, Gradient, BlockID, PlotID, Trait)
  
  # Site level weights and traits
  TraitWeights_site <- comm %>%
    left_join(trait %>% select(-Year, -BlockID, -PlotID), by = c("Country", "Site", "Gradient", "Taxon")) %>% 
    group_by(Country, Year, Site, Gradient, Taxon, Trait) %>% 
    mutate(weight = Cover/n()) %>% 
    group_by(Country, Year, Site, Gradient, Trait) 
  
  
  # Global level weights and traits
  TraitWeights_global <- comm %>% 
    left_join(trait %>% select(-Year, -BlockID, -PlotID, -Site), by = c("Country", "Gradient", "Taxon")) %>% 
    group_by(Country, Year, Gradient, Taxon, Trait) %>% 
    mutate(weight = Cover/n()) %>% 
    group_by(Country, Year, Gradient, Trait) 
  
  
  TraitWeights_all <- bind_rows(plot = TraitWeights_plot, site = TraitWeights_site, global = TraitWeights_global, .id = "level") %>% 
    mutate(level = factor(level, levels = c("plot", "site", "global"), ordered = TRUE)) %>%
    filter(!is.na(Value)) %>% 
    group_by(Country, Year, Site, Gradient, BlockID, PlotID, Trait, Taxon) %>% 
    filter(level == min(level)) %>% 
    group_by(Country, Year, Site, Gradient, BlockID, PlotID, Trait)
  
  
  BootstrapMoments_All <- rerun(.n = nrep, 
                                slice_sample(TraitWeights_all, 
                                         n = samplesize,  
                                         replace = TRUE, 
                                         weight_by = TraitWeights_all$weight)
                                ) %>%
    bind_rows(.id = "n") %>% 
    group_by(n, add = TRUE) %>% 
    # get all the happy moments
    summarise(Mean = mean(Value), Variance = var(Value), Skewness = skewness(Value), Kurtosis = kurtosis(Value))
  
  return(BootstrapMoments_All)
}

SummarizeBootMoments <- function(BootstrapMoments_All){
  # calculate means and 
  BootstrapMoments <- BootstrapMoments_All %>% 
    group_by(Country, Year, Site, Gradient, BlockID, PlotID, Trait) %>% 
    summarise(n = n(),
              meanMean = mean(Mean), CIlow.Mean = meanMean - sd(Mean), CIhigh.Mean = meanMean + sd(Mean),
              meanVar = mean(Variance), CIlow.Var = meanVar - sd(Variance), CIhigh.Var = meanVar + sd(Variance),
              meanSkew = mean(Skewness), CIlow.Skew = meanSkew - sd(Skewness), CIhigh.Skew = meanSkew + sd(Skewness),
              meanKurt = mean(Kurtosis), CIlow.Kurt = meanKurt - sd(Kurtosis), CIhigh.Kurt = meanKurt + sd(Kurtosis)) 
  
  return(BootstrapMoments)
}


#NEEDS TO BE INCORPORATED IN BOOTSTRAPPING FUNCTION !!!!!
#TraitWeights_all %>% 
#group_by(Country, Year, Site, Gradient, BlockID, PlotID, Trait, Taxon) %>% 
#slice(1) %>% 
#group_by(Country, Year, Site, Gradient, BlockID, PlotID, Trait) %>% summarise(covered = sum(Cover)) %>% mutate(percent_cover = covered/sumCover * 100) %>% arrange(percent_cover) %>% pn


#### Filtering out turfs with less than 70% of the community present ###

#check_community_df <- wcommunity %>%
#group_by(Site, Species, turfID)%>%
#select(Site, turfID, Species, cover, SLA_mean, Lth_mean, Height_mean, LDMC_mean, LA_mean, CN_ratio_mean, sum_cover)%>%
#unique()%>%
#ungroup()%>%
#group_by(turfID)%>%
#mutate(cover_traits = (sum(cover)))%>%
#filter(!is.na(SLA_mean))%>%
#mutate(community_covered_trait=cover_traits/sum_cover*100)

#complete_turf <- check_community_df%>%
#filter(community_covered_trait>80)%>%
#distinct(turfID, .keep_all=TRUE)

#Complete_turfs<-as.vector(complete_turf$turfID)

#wcommunity_df <- filter(wcommunity, turfID %in% Complete_turfs)




