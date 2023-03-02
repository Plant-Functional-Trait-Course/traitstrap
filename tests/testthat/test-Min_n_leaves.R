context("trait_fill minimum number in sample")

test_that("trait_fill minimum number in sample", {
  #### set-up ####
  mini_trait <- tidyr::crossing(
    genus = c("G1", "G2"),
    taxon = c("sp1", "sp2"),
    site = c("A", "B"),
    plot = 1:2,
    trait = "trait",
    n_sample = 1:5 # 5 samples per species-trait
  ) |>
    mutate(taxon = paste(genus, taxon)) |>
    mutate(value = 1:80)

  mini_comm <- tidyr::crossing(
    genus = c("G1", "G2"),
    taxon = c("sp1", "sp2"),
    site = c("A", "B"),
    plot = 1:2
  ) |>
    mutate(taxon = paste(genus, taxon)) |>
    mutate(cover = 5)

  #### test set 1 ####
  # site A plot 2 G1 sp1 has only 2 samples
  # fill should include samples from plot 2

  mini_trait1 <- mini_trait |>
    filter(!(taxon == "G1 sp1" & site == "A" & plot == 2 & n_sample > 2))

  ti_1 <- trait_fill(
    comm = mini_comm,
    traits = mini_trait1,
    scale_hierarchy = c("site", "plot"),
    taxon_col = "taxon",
    trait_col = "trait",
    value_col = "value",
    abundance_col = "cover",
    min_n_in_sample = 5
  )

  # check expected number of samples (two from plot 2, five from plot 2)
  cond <- with(ti_1, (taxon == "G1 sp1" & site == "A" & plot == 2))
  expect_equal(
    unique(ti_1[cond, "n_sample", drop = TRUE]),
    5 + 2
  )

  # check selection from correct level (site)
  expect_equal(
    as.character(unique(ti_1[cond, "level", drop = TRUE])),
    "site"
  )

  #### test set 2 ####
  # site A G1 sp1 has only 2 samples in each plot = 4 overall
  # fill should be at global level

  mini_trait2 <- mini_trait |>
    filter(!(taxon == "G1 sp1" & site == "A" & n_sample > 2))

  ti_2 <- trait_fill(
    comm = mini_comm,
    traits = mini_trait2,
    scale_hierarchy = c("site", "plot"),
    taxon_col = "taxon",
    trait_col = "trait",
    value_col = "value",
    abundance_col = "cover",
    min_n_in_sample = 5
  )

  # check expected number in sample
  # (2 each from A plots, 5 each from B plots) = 14
  cond <- with(ti_2, (taxon == "G1 sp1" & site == "A" & plot == 2))
  expect_equal(
    unique(ti_2[cond, "n_sample", drop = TRUE]),
    (5 + 2) * 2
  )

  # check selection from correct level (global)
  expect_equal(
    as.character(unique(ti_2[cond, "level", drop = TRUE])),
    "global"
  )

  #### test set 3 ####
  # G1 sp1 has only 1 sample in each plot = 4 overall !< min_n_in_sample
  # fill should be at global level

  mini_trait3 <- mini_trait |>
    filter(!(taxon == "G1 sp1" & n_sample > 1))

  ti_3 <- trait_fill(
    comm = mini_comm,
    traits = mini_trait3,
    scale_hierarchy = c("site", "plot"),
    taxon_col = "taxon",
    trait_col = "trait",
    value_col = "value",
    abundance_col = "cover",
    min_n_in_sample = 5
  )

  # check expected number of leaves (one from each plot) == 4
  cond <- with(ti_3, (taxon == "G1 sp1" & site == "A" & plot == 2))
  expect_equal(
    unique(ti_3[cond, "n_sample", drop = TRUE]),
    1 * 4
  )

  # check selection from correct level (global)
  expect_equal(
    as.character(unique(ti_3[cond, "level", drop = TRUE])),
    "global"
  )
})
