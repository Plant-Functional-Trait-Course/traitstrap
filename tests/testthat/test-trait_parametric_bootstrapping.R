context("trait_parametric_bootstrap")

test_that("output is properly formatted", {
  # Get test data
  data(community)
  data(trait)
  nrep_pbs <- 10
  filled_traits <- trait_fill(
    comm = community,
    traits = trait,
    scale_hierarchy = c("Site", "PlotID"),
    global = TRUE,
    taxon_col = "Taxon",
    trait_col = "Trait",
    min_n_in_sample = 3
  )

  fitted_distributions <- trait_fit_distributions(
    filled_traits = filled_traits, distribution_type = "lognormal"
  )

  pbs_out <- trait_parametric_bootstrap(
    fitted_distributions = fitted_distributions, nrep = nrep_pbs
  )

  expect_true("tbl" %in% class(pbs_out))

  expect_true(all(c("mean", "variance", "skewness", "kurtosis")
  %in% colnames(pbs_out)))

  expect_true("tbl" %in%
    class(trait_summarise_boot_moments(bootstrap_moments = pbs_out)))

  scale_hierarchy <- attr(fitted_distributions, "attr")$scale_hierarchy
  taxon_col <- attr(fitted_distributions, "attr")$taxon_col
  trait_col <- attr(fitted_distributions, "attr")$trait_col


  expect_equal(
    object = nrow(pbs_out),
    expected = filled_traits |>
      group_by_at(c(
        as.character(scale_hierarchy),
        trait_col
      )) |>
      dplyr::n_groups() * nrep_pbs
  )
})


test_that("Bad inputs return errors", {
  # Get test data
  data(community)
  data(trait)
  nrep_pbs <- 10
  filled_traits <- trait_fill(
    comm = community,
    traits = trait,
    scale_hierarchy = c("Site", "PlotID"),
    global = TRUE,
    taxon_col = "Taxon",
    trait_col = "Trait",
    min_n_in_sample = 3
  )

  fitted_distributions <- trait_fit_distributions(
    filled_traits = filled_traits, distribution_type = "lognormal"
  )

  expect_error(object = trait_parametric_bootstrap(
    fitted_distributions = "Wreck-Gar", nrep = 10
  ))

  expect_error(object = trait_parametric_bootstrap(
    fitted_distributions = fitted_distributions, nrep = "Arcee"
  ))

  expect_error(object = trait_parametric_bootstrap(
    fitted_distributions = fitted_distributions,
    nrep = 10, sample_size = "Kup"
  ))
})
