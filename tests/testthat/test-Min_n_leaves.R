context("trait_impute minimum number of leaves")

test_that("trait_impute minimum number of leaves", {
  #### set-up ####
  mini_trait <- tidyr::crossing(
    genus = c("G1", "G2"),
    taxon = c("sp1", "sp2"),
    site = c("A", "B"),
    plot = 1:2,
    trait = "trait",
    leaf_n = 1:5 # 5 leaves per species-trait
  ) %>%
    mutate(taxon = paste(genus, taxon)) %>%
    mutate(value = 1:80)

  mini_comm <- tidyr::crossing(
    genus = c("G1", "G2"),
    taxon =  c("sp1", "sp2"),
    site = c("A", "B"),
    plot = 1:2
  ) %>%
    mutate(taxon = paste(genus, taxon)) %>%
    mutate(cover = 5)

  #### test set 1 ####
  # site A plot 2 G1 sp1 has only 2 leaves
  # impute should include leaves from from plot 2

  mini_trait1 <- mini_trait %>%
    filter(!(taxon == "G1 sp1" & site == "A" & plot == 2 & leaf_n > 2))

  ti_1 <- trait_impute(
    comm = mini_comm,
    traits = mini_trait1,
    scale_hierarchy = c("site", "plot"),
    taxon_col = "taxon",
    trait_col = "trait",
    value_col = "value",
    abundance_col = "cover",
    min_n_leaves = 5
  )

  #check expected number of leaves (two from plot 2, five from plot 2)
  cond <- with(ti_1, (taxon == "G1 sp1" & site == "A" & plot == 2))
  expect_equal(
    unique(ti_1[cond, "n_leaves", drop = TRUE]),
    5 + 2
  )

  #check imputation from correct level (site)
  expect_equal(
    as.character(unique(ti_1[cond, "level", drop = TRUE])),
    "site"
  )

  #### test set 2 ####
  # site A G1 sp1 has only 2 leaves in each plot = 4 overall
  # impute should be at global level

  mini_trait2 <- mini_trait %>%
    filter(!(taxon == "G1 sp1" & site == "A" & leaf_n > 2))

  ti_2 <- trait_impute(
    comm = mini_comm,
    traits = mini_trait2,
    scale_hierarchy = c("site", "plot"),
    taxon_col = "taxon",
    trait_col = "trait",
    value_col = "value",
    abundance_col = "cover",
    min_n_leaves = 5
  )

  #check expected number leaves (2 each from A plots, 5 each from B plots) = 14
  cond <- with(ti_2, (taxon == "G1 sp1" & site == "A" & plot == 2))
  expect_equal(
    unique(ti_2[cond, "n_leaves", drop = TRUE]),
    (5 + 2) * 2
  )

  #check imputation from correct level (global)
  expect_equal(
    as.character(unique(ti_2[cond, "level", drop = TRUE])),
    "global"
  )

  #### test set 3 ####
  # G1 sp1 has only 1 leaf in each plot = 4 overall !< n_min_leaves
  # impute should be at global level

  mini_trait3 <- mini_trait %>%
    filter(!(taxon == "G1 sp1" & leaf_n > 1))

  ti_3 <- trait_impute(
    comm = mini_comm,
    traits = mini_trait3,
    scale_hierarchy = c("site", "plot"),
    taxon_col = "taxon",
    trait_col = "trait",
    value_col = "value",
    abundance_col = "cover",
    min_n_leaves = 5
  )

  #check expected number of leaves (one from each plot) == 4
  cond <- with(ti_3, (taxon == "G1 sp1" & site == "A" & plot == 2))
  expect_equal(
    unique(ti_3[cond, "n_leaves", drop = TRUE]),
    1 * 4
  )

  #check imputation from correct level (global)
  expect_equal(
    as.character(unique(ti_3[cond, "level", drop = TRUE])),
    "global"
  )
})
