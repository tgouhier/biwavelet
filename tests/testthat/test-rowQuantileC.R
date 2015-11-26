context("rowQuantileC implementation in C++")

test_that("Quantile out of range should fail", {
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

test_that("R and C++ versions of quantile using 1xN matrix shold be equal", {
  data <- matrix(rnorm(100), nrow = 1)
  for (p in seq(0, 1, length.out = 20)) {
    expect_equal(rowQuantileC(data, p),
                 as.double(quantile(data, p)),
                 tolerance = 1e-10)
  }
})

test_that("rowQuantileC should ignore NAs without errors", {
  data <- matrix(c(1,1,NA,2,3,3,3,5,NA,NA), nrow = 1)
  for (p in seq(0, 1, length.out = 20)) {
    expect_equal(rowQuantileC(data, 0.9),
                 as.double(quantile(data, 0.9, na.rm = TRUE)),
                 tolerance = 1e-10)
  }
})

test_that("Multiple quantiles as parameter should fail", {
  data <- matrix(1:10, nrow = 2)
  expect_error(rowQuantileC(data, c(.5, .75)))
  expect_error(rowQuantileC(data, 1:5))
})

test_that("Multiple quantiles as parameter are not supported", {
  data <- matrix(1:10, nrow = 2)
  expect_error(rowQuantileC(data, c(.5, .75)), regexp = "one quantile")
  expect_error(rowQuantileC(data, 1:5), regexp = "one quantile")
})



test_that("", {
   data <- as.double(c(1,1,1,1,2,3,4))
   .C("pivot", v = data, len = length(data), target = 3)
 })
