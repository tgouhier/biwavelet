context("rcpp_row_quantile implementation in C++")

test_that("Quantile out of range should fail", {
  data <- matrix(1:10, nrow = 2)
  expect_error(rcpp_row_quantile(data, 2), regexp = "out of range")
  expect_error(rcpp_row_quantile(data, -1), regexp = "out of range")
  expect_error(rcpp_row_quantile(data, 1.1), regexp = "out of range")
})

test_that("R and C++ versions of quantile using 1xN matrix shold be equal", {
  data <- matrix(rnorm(100), nrow = 1)
  for (p in seq(0, 1, length.out = 20)) {
    expect_equal(rcpp_row_quantile(data, p),
                 as.double(quantile(data, p)),
                 tolerance = 1e-10)
  }
})

test_that("R and C++ version work the same on small inputs", {
  for (data in list(
    matrix(1),
    matrix(1:2, nrow = 1),
    matrix()
  )) {
    expect_equal(rcpp_row_quantile(data, .2),
                 as.double(quantile(data, .2, na.rm = TRUE)),
                 tolerance = 1e-10)
  }
})

test_that("rcpp_row_quantile should ignore NAs without errors", {
  data <- matrix(c(1,1,NA,2,3,3,3,5,NA,NA), nrow = 1)
  for (p in seq(0, 1, length.out = 20)) {
    expect_equal(rcpp_row_quantile(data, 0.9),
                 as.double(quantile(data, 0.9, na.rm = TRUE)),
                 tolerance = 1e-10)
  }
})

test_that("Matrix without columns should produce NAs", {
  data <- matrix(nrow = 10, ncol = 0)
  expect_equal(rcpp_row_quantile(data, .7), rep(NA_real_, 10))
})

test_that("Multiple quantiles as parameter are not supported", {
  data <- matrix(1:10, nrow = 2)
  expect_error(rcpp_row_quantile(data, c(.5, .75)), regexp = "single value")
  expect_error(rcpp_row_quantile(data, 1:5), regexp = "single value")
})
