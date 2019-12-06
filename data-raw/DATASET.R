library("tidyverse")
#### community ####
load("data-raw/community.RData")
community <- community %>% 
  mutate(Cover = as.numeric(Cover)) %>% 
  filter(!is.na(Cover))
usethis::use_data(community)

#### trait ####
load("data-raw/trait.RData")
trait <- trait %>% filter(!is.na(Value))
usethis::use_data(trait)

