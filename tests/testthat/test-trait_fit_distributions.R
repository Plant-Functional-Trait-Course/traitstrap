context("trait_fit_distributions")

test_that("returns parametric distribution tibbles of proper size", {
  # Get test data

  data(community)

  data(trait)

  filled_traits <- trait_fill(
    comm = community,
    traits = trait,
    scale_hierarchy = c("Site", "PlotID"),
    global = TRUE,
    taxon_col = "Taxon",
    trait_col = "Trait",
    min_n_in_sample = 3
  )



  expect_true(all(c("parametric_distributions", "tbl") %in%

    class(trait_fit_distributions(
      filled_traits = filled_traits,
      distribution_type = "lognormal"
    ))))

  filled_traits$Value <- log10(filled_traits$Value)

  expect_true(all(c("parametric_distributions", "tbl") %in%
    class(trait_fit_distributions(
      filled_traits = filled_traits,
      distribution_type = "normal"
    ))))

  filled_traits$Value <- rbeta(
    n = nrow(filled_traits),
    shape1 = .5, shape2 = .5
  )

  expect_true(all(c("parametric_distributions", "tbl") %in%
    class(trait_fit_distributions(
      filled_traits = filled_traits,
      distribution_type = "beta"
    ))))


  scale_hierarchy <- attr(filled_traits, "attr")$scale_hierarchy
  taxon_col <- attr(filled_traits, "attr")$taxon_col
  trait_col <- attr(filled_traits, "attr")$trait_col

  expect_equal(
    object = nrow(trait_fit_distributions(
      filled_traits =
        filled_traits,
      distribution_type = "beta"
    )),
    expected = filled_traits |>
      group_by_at(c(
        as.character(scale_hierarchy),
        taxon_col, trait_col
      )) |>
      dplyr::n_groups()
  )
})


test_that("bad inputs return errors", {
  # Get test data

  data(community)

  data(trait)

  filled_traits <- trait_fill(
    comm = community,
    traits = trait,
    scale_hierarchy = c("Site", "PlotID"),
    global = TRUE,
    taxon_col = "Taxon",
    trait_col = "Trait",
    min_n_in_sample = 3
  )

  expect_error(object = trait_fit_distributions(
    filled_traits = "a", distribution_type = "Megatron"
  ))

  expect_error(object = trait_fit_distributions(
    filled_traits = filled_traits, distribution_type = "Soundwave"
  ))
})
