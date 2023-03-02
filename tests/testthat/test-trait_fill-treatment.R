context("trait_fill with treatment")

test_that("trait_fill with treatment", {
  #### set-up ####
  mini_trait <- tidyr::crossing(
    taxon = c("sp1", "sp2"),
    site = c("A", "B"),
    plot = 1:3,
    trait = "trait",
  ) |>
    mutate(
      treatment = dplyr::recode(plot, `1` = "a", `2` = "b", `3` = "c"),
      treatment = factor(treatment, levels = c("c", "a", "b")),
      value = 1:12
    )

  mini_comm <- tidyr::crossing(
    taxon = c("sp1", "sp2"),
    site = c("A", "B"),
    plot = 1:3,
    cover = 5
  ) |>
    mutate(
      treatment = dplyr::recode(plot, `1` = "a", `2` = "b", `3` = "c"),
      treatment = factor(treatment, levels = c("c", "a", "b"))
    )

  #### test set 1 ####
  # fill site A plot 2 (treat b) - should get value from A2b and then A3c.
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
    treatment_col = "treatment",
    treatment_level = "site",
    min_n_in_sample = 1
  )

  # check expected value of trait filled (A1 sp1)
  got <- ti_1 |> filter(site == "A", taxon == "sp1", plot == 2)
  target <- mini_trait |> filter(site == "A", taxon == "sp1", treatment == "c")
  expect_equal(
    sort(got$value),
    sort(target$value)
  )
})
