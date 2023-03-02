context("trait_fill taxonomic filling")

test_that("trait_fill taxonomic filling", {
  #### set-up ####
  mini_trait <- tidyr::crossing(
    genus = c("G1", "G2"),
    taxon = c("sp1", "sp2"),
    site = c("A", "B"),
    plot = 1:2,
    trait = "trait"
  ) |>
    mutate(taxon = paste(genus, taxon)) |>
    mutate(value = 1:16)

  mini_comm <- tidyr::crossing(
    genus = c("G1", "G2"),
    taxon = c("sp1", "sp2"),
    site = c("A", "B"),
    plot = 1:2
  ) |>
    mutate(taxon = paste(genus, taxon)) |>
    mutate(cover = 5)

  #### test set 1 ####
  # fill site A plot 2 - should get value from A1.
  # Should be no genus level filling
  mini_trait1 <- mini_trait |>
    filter(!(taxon == "G1 sp1" & site == "A" & plot == 2))

  ti_1 <- trait_fill(
    comm = mini_comm,
    traits = mini_trait1,
    scale_hierarchy = c("site", "plot"),
    taxon_col = c("taxon", "genus"),
    trait_col = "trait",
    value_col = "value",
    abundance_col = "cover",
    min_n_in_sample = 1
  )

  # check expected value of trait filled (A1 sp1)
  cond <- with(ti_1, (taxon == "G1 sp1" & site == "A" & plot == 2))
  expect_equal(
    ti_1[cond, "value"],
    mini_trait[(mini_trait$taxon == "G1 sp1" &
      mini_trait$site == "A" &
      mini_trait$plot == 1), "value"]
  )

  # check filling from correct level (site)
  expect_equal(
    as.vector(ti_1[cond, "level", drop = TRUE]),
    "site"
  )
  # check filling from correct level (site)
  expect_equal(
    ti_1[cond, "weight", drop = TRUE],
    mini_comm[(mini_trait$taxon == "G1 sp1" &
      mini_trait$site == "A" &
      mini_trait$plot == 1), "cover", drop = TRUE]
  )
  # check all filling at correct taxon level (taxon)
  expect_true(
    all(ti_1[cond, "taxon_level", drop = TRUE] == "taxon")
  )

  #### test set 1b ####
  # exactly the same as test set 1a but with pipes
  ti_1 <- trait_fill(
    comm = mini_comm,
    traits = mini_trait1,
    scale_hierarchy = c("site", "plot"),
    taxon_col = c("taxon", "genus"),
    trait_col = "trait",
    value_col = "value",
    abundance_col = "cover",
    min_n_in_sample = 1
  )

  # check expected value of trait filled (A1 sp1)
  cond <- with(ti_1, (taxon == "G1 sp1" & site == "A" & plot == 2))
  expect_equal(
    ti_1[cond, "value"],
    mini_trait[(mini_trait$taxon == "G1 sp1" &
      mini_trait$site == "A" &
      mini_trait$plot == 1), "value"]
  )

  # check filling from correct level (site)
  expect_equal(
    as.vector(ti_1[cond, "level", drop = TRUE]),
    "site"
  )
  # check filling from correct level (site)
  expect_equal(
    ti_1[cond, "weight", drop = TRUE],
    mini_comm[(mini_trait$taxon == "G1 sp1" &
      mini_trait$site == "A" &
      mini_trait$plot == 1), "cover", drop = TRUE]
  )
  # check all filling at correct taxon level (taxon)
  expect_true(
    all(ti_1[cond, "taxon_level", drop = TRUE] == "taxon")
  )


  #### test set 2 ####
  # fill site A G1 sp1 - should get value from site A G1 sp2
  mini_trait2 <- mini_trait |>
    filter(taxon != "G1 sp1")

  ti_2 <- trait_fill(
    comm = mini_comm,
    traits = mini_trait2,
    scale_hierarchy = c("site", "plot"),
    taxon_col = c("taxon", "genus"),
    trait_col = "trait",
    value_col = "value",
    abundance_col = "cover",
    min_n_in_sample = 1
  )

  # for site A, should fill from G1 sp2 at site A
  got <- ti_2 |> filter(site == "A", taxon == "G1 sp1")
  target <- mini_trait |> filter(site == "A", taxon == "G1 sp2")
  expect_equal(
    sort(got$value),
    sort(target$value)
  )
})
