context("get_ci")

test_that("get_ci ranks are appropriate", {
  test_vector <- runif(n = 100, min = 0, max = 100)

  expect_gt(
    object = get_ci(data = test_vector, which = "high"),
    expected = get_ci(data = test_vector, which = "low")
  )

  expect_gt(
    object = get_ci(data = test_vector, which = "high", sd_mult = 2),
    expected = get_ci(data = test_vector, which = "high")
  )

  expect_gt(
    object = get_ci(data = test_vector, which = "high", parametric = FALSE),
    expected = get_ci(data = test_vector, which = "low", parametric = FALSE)
  )

  expect_gt(
    object = get_ci(
      data = test_vector, which = "low",
      parametric = FALSE, ci = 0.45
    ),
    expected = get_ci(
      data = test_vector, which = "low",
      parametric = FALSE, ci = 0.99
    )
  )
})

test_that("get_ci methods return same class and type", {
  test_vector <- runif(n = 100, min = 0, max = 100)
  expect_equal(
    object = typeof(get_ci(data = test_vector, which = "high")),
    expected = typeof(get_ci(data = test_vector, which = "high"))
  )
})

test_that("get_ci data order doesn't matter", {
  test_vector <- runif(n = 100, min = 0, max = 100)

  expect_equal(
    object = get_ci(data = test_vector[100:1], which = "high"),
    expected = get_ci(data = test_vector[1:100], which = "high")
  )

  expect_equal(
    object = get_ci(
      data = test_vector[100:1],
      which = "high", parametric = FALSE
    ),
    expected = get_ci(
      data = test_vector[1:100],
      which = "high", parametric = FALSE
    )
  )
})
