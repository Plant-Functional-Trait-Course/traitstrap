context("get_dist_parms")

test_that("returns dataframes", {
  test_vector <- rnorm(n = 100)
  expect_equal(
    object = class(get_dist_parms(
      data = test_vector, distribution_type = "normal"
    )),
    expected = "data.frame"
  )

  test_vector <- rlnorm(n = 100)
  expect_equal(
    object = class(get_dist_parms(
      data = test_vector, distribution_type = "lognormal"
    )),
    expected = "data.frame"
  )

  test_vector <- rbeta(n = 100, shape1 = .5, shape2 = .5)
  expect_equal(
    object = class(get_dist_parms(
      data = test_vector, distribution_type = "beta"
    )),
    expected = "data.frame"
  )
})

test_that("nonsuported distributions return Error", {
  test_vector <- rnorm(n = 100)
  expect_error(object = get_dist_parms(
    data = test_vector, distribution_type = "OptimusPrime"
  ))
})
