context("trait_fill")

test_that("trait_fill works", {
  #### setup ####
  mini_trait <- tidyr::crossing(
    taxon = c("sp1", "sp2"),
    site = c("A", "B"),
    plot = 1:2,
    trait = "trait"
  ) |>
    mutate(value = 1:8)

  mini_comm <- tidyr::crossing(
    taxon =  c("sp1", "sp2"),
    site = c("A", "B"),
    plot = 1:2
  ) |>
    mutate(cover = 5)

  #### test set 1 ####
  # select site A plot 2 - should get value from A1
  mini_trait1 <- mini_trait |>
    filter(!(taxon == "sp1" & site == "A" & plot == 2))

  ti_1 <- trait_fill(
    comm = mini_comm,
    traits = mini_trait1,
    scale_hierarchy = c("site", "plot"),
    taxon_col = "taxon",
    trait_col = "trait",
    value_col = "value",
    abundance_col = "cover",
    min_n_in_sample = 1
  )

  # check expected value of trait selected (A1 sp1)
  cond <- with(ti_1, (taxon == "sp1" & site == "A" & plot == 2))
  expect_equal(
    ti_1[cond, "value"],
    mini_trait[
      (mini_trait$taxon == "sp1" &
        mini_trait$site == "A" &
        mini_trait$plot == 1),
      "value"
    ]
  )
  # check selection from correct level (site)
  expect_equal(
    as.vector(ti_1[cond, "level", drop = TRUE]),
    "site"
  )
  # check selection from correct level (site)
  expect_equal(
    ti_1[cond, "weight", drop = TRUE],
    mini_comm[(mini_trait$taxon == "sp1" &
      mini_trait$site == "A" &
      mini_trait$plot == 1), "cover", drop = TRUE]
  )

  #### test set 2 ####
  # select site A sp1 - should get value from site B
  mini_trait2 <- mini_trait |>
    filter(!(taxon == "sp1" & site == "A"))

  ti_2 <- trait_fill(
    comm = mini_comm,
    traits = mini_trait2,
    scale_hierarchy = c("site", "plot"),
    taxon_col = "taxon",
    trait_col = "trait",
    value_col = "value",
    abundance_col = "cover",
    min_n_in_sample = 1
  )

  # check 2 entries for each sp1 in each plot in A (leave from each plot in B)
  expect_equal(
    nrow(with(ti_2, ti_2[site == "A" & taxon == "sp1" & plot == 1, ])),
    2
  )

  # check expected value of trait selected (B1B2 sp1)
  cond <- with(ti_2, (taxon == "sp1" & site == "A"))
  expect_equal(
    ti_2[cond, "value"][1:2, ],
    mini_trait[(mini_trait$taxon == "sp1" & mini_trait$site == "B"), "value"]
  )
  # check selection from correct level (global)
  expect_equal(
    as.vector(ti_2[cond, "level", drop = TRUE])[1],
    "global"
  )
  # check selection from correct level (site)
  expect_equal(
    ti_2[cond, "weight", drop = TRUE][1],
    mini_comm[(mini_trait$taxon == "sp1" &
      mini_trait$site == "A" &
      mini_trait$plot == 1), "cover", drop = TRUE] / 2
  )
})
