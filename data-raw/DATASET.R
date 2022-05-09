library("tidyverse")
library("dataDownloader")

#### community ####
load("data-raw/community.RData")
community <- community %>% 
  mutate(Cover = as.numeric(Cover)) %>% 
  filter(!is.na(Cover))
usethis::use_data(community)

#### trait ####
dataDownloader::get_file(node = "7mzjk", remote_path = "Svalbard", file = "traitsGradients_SV_2018.Rdata", path = "data-raw/")

load("data-raw/traitsGradients_SV_2018.Rdata")

trait <- traitsGradients_SV_2018 |> 
  filter(Project == "T", Site %in% 1:2, Gradient == "C") |> 
  select(Site, PlotID, Taxon, ID, Plant_Height_cm, Wet_Mass_g, Leaf_Thickness_Ave_mm) |> 
  pivot_longer(c(Plant_Height_cm, Wet_Mass_g, Leaf_Thickness_Ave_mm), 
               names_to = "Trait", 
               values_to = "Value", 
               values_drop_na = TRUE)

usethis::use_data(trait, overwrite = TRUE)
