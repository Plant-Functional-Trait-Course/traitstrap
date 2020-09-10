context("trait_impute")

test_that("trait_imputation works", {

    mini_trait <- tidyr::crossing(
      taxon = c("sp1", "sp2"), 
      site = c("A", "B"), 
      plot = 1:2, 
      trait = "trait"
    ) %>% 
    mutate(value = 1:8)
  
  mini_comm <- tidyr::crossing(
    taxon =  c("sp1", "sp2"),
    site = c("A", "B"),
    plot = 1:2
  ) %>%
    mutate(cover = 5)
  
# impute site A plot 2 - should get value from A1
  mini_trait2 <- mini_trait %>% 
    filter(!(taxon == "sp1" & site == "A" & plot == 2))
  
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
  cond <- with(ti2, (taxon == "sp1" & site == "A" & plot == 2))
  expect_equal(
    ti2[cond, "value"], 
    mini_trait[(mini_trait$taxon == "sp1" & mini_trait$site == "A" & mini_trait$plot == 1) , "value"]
  )
  #check imputation from correct level (site)
  expect_equal(
    as.vector(ti2[cond, "level", drop = TRUE]), 
    "site"
  )
  #check imputation from correct level (site)
  expect_equal(
    ti2[cond, "weight", drop = TRUE], 
    mini_comm[(mini_trait$taxon == "sp1" & mini_trait$site == "A" & mini_trait$plot == 1) , "cover", drop = TRUE]
  )
  
})


