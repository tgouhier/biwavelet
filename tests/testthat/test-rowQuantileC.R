context("rowQuantileC implementation in C++")

test_that("Quantile out of range 0 to 1 should fail", {
  data <- matrix(1:10, nrow = 2)

  expect_error(rowQuantileC(data, 2), regexp = "out of the allowed range")
  expect_error(rowQuantileC(data, -1), regexp = "out of the allowed range")
  expect_error(rowQuantileC(data, 1.1), regexp = "out of the allowed range")

  expect_error(
    .C("rowQuantileC", data = as.double(data),
       nrow = nrow(data), ncol = ncol(data),
       p = 2, quantiles = rep(0, nrow(data))),
    regexp = "out of range"
  )
})

test_that("Multiple quantiles as parameter should also fail", {
  data <- matrix(1:10, nrow = 2)
  expect_error(rowQuantileC(data, c(.5, .75)))
  expect_error(rowQuantileC(data, 1:5))
})
