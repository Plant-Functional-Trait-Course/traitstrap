context("trait_fit_distributions")

test_that("returns parametric distribution tibbles of proper size", {

  #Get test data
  data(community)
  data(trait)
  selected_traits <- trait_select(comm = community,
                                 traits = trait,
                                 scale_hierarchy = c("Site", "PlotID"),
                                 global = TRUE,
                                 taxon_col = "Taxon",
                                 trait_col = "Trait",
                                 min_n_in_sample = 3)

  expect_true(all(c("parametric_distributions", "tbl") %in%
                    class(trait_fit_distributions(
                      selected_traits = selected_traits,
                      distribution_type = "lognormal"))))

  selected_traits$Value <- log10(selected_traits$Value)
  expect_true(all(c("parametric_distributions", "tbl") %in%
                    class(trait_fit_distributions(
                      selected_traits = selected_traits,
                      distribution_type = "normal"))))

  selected_traits$Value <- rbeta(n = nrow(selected_traits),
                                shape1 = .5, shape2 = .5)
  expect_true(all(c("parametric_distributions", "tbl") %in%
                    class(trait_fit_distributions(
                      selected_traits = selected_traits,
                      distribution_type = "beta"))))


  scale_hierarchy <- attr(selected_traits, "attr")$scale_hierarchy
  taxon_col <- attr(selected_traits, "attr")$taxon_col
  trait_col <- attr(selected_traits, "attr")$trait_col


  expect_equal(object = nrow(trait_fit_distributions(selected_traits =
                                              selected_traits,
                                            distribution_type = "beta")),
               expected = selected_traits %>%
                 group_by_at(c(as.character(scale_hierarchy),
                               taxon_col, trait_col)) %>%
                 dplyr::n_groups())

})


test_that("bad inputs return errors", {

  #Get test data
  data(community)
  data(trait)
  selected_traits <- trait_select(comm = community,
                                 traits = trait,
                                 scale_hierarchy = c("Site", "PlotID"),
                                 global = TRUE,
                                 taxon_col = "Taxon",
                                 trait_col = "Trait",
                                 min_n_in_sample = 3)

  expect_error(object = trait_fit_distributions(
    selected_traits = "a", distribution_type = "Megatron"))

  expect_error(object = trait_fit_distributions(
    selected_traits = selected_traits, distribution_type = "Soundwave"))


})
