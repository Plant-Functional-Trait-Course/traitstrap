context("distribution_handler")

test_that("returns numeric", {
  expect_equal(
    object = class(distribution_handler(
      parm1 = 0, parm2 = 1, n = 1, type = "normal"
    )),
    expected = "numeric"
  )

  expect_equal(
    object = class(distribution_handler(
      parm1 = 0, parm2 = 1, n = 10, type = "lognormal"
    )),
    expected = "numeric"
  )

  expect_equal(
    object = class(distribution_handler(
      parm1 = .5, parm2 = .5, n = 10, type = "beta"
    )),
    expected = "numeric"
  )
})

test_that("nonsupported distributions return NULL", {
  expect_null(object = distribution_handler(
    parm1 = 0, parm2 = 1, n = 1, type = "Unicron"
  ))
})
