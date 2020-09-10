context("trait_impute")

test_that("trait_imputation works", {
library(tidyr)
    mini_trait <- crossing(
      taxon = c("sp1", "sp2"), 
      site = c("A", "B"), 
      plot = 1:2, 
      trait = "trait"
    ) %>% 
    mutate(value = 1:8)
  
  mini_comm <- crossing(
    taxon =  c("sp1", "sp2"),
    site = c("A", "B"),
    plot = 1:2
  ) %>%
    mutate(cover = 5)
  
# impute site A plot 2 - should get value from A1
  cond <- with(mini_trait, !(taxon == "sp1" & site == "A" & plot == 2))
  mini_trait2 <- mini_trait %>% 
    filter(cond)
  
  ti2 <- trait_impute(
    comm = mini_comm, 
    traits = mini_trait2,
    scale_hierarchy = c("site", "plot"), 
    taxon_col = "taxon",
    trait_col = "trait",
    value_col = "value",
    abundance_col = "cover"
  ) 
  
  #check expected value of trait imputed (A1 sp1)
  expect_equal(
    ti2[!cond, "value"], 
    mini_trait[(mini_trait$taxon == "sp1" & mini_trait$site == "A" & mini_trait$plot == 1) , "value"]
  )
})


