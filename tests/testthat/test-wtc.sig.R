# Setup variables =====

size <- 64
d1 <- cbind(1:size, rnorm(size))
d2 <- cbind(1:size, rnorm(size))

d1.ar1 <- arima(d1[,2], order = c(1, 0, 0))$coef[1]
d2.ar1 <- arima(d2[,2], order = c(1, 0, 0))$coef[1]

checked <- check.data(y = d1, x1 = d2)
dt <- checked$y$dt
n <- checked$y$n.obs
s0 <- 2 * dt

# Tests ==========
context("wtc.sig")

test_that("Quiet mode without progress bar should not throw errors", {
  out <- wtc.sig(quiet = TRUE, nrands = 2, lag1 = c(d1.ar1, d2.ar1),
                 dt = dt, ntimesteps = n, s0 = s0, J1 = NULL)
  expect_true( is.matrix(out) )
})

test_that("Error message for unsupported mother wavelet", {
  expect_error(
    wtc.sig( mother = "dummy",
             lag1 = c(d1.ar1, d2.ar1),
             dt = dt, ntimesteps = n, s0 = s0, J1 = NULL),
    regexp = "should be one of" )
})

test_that("Testing whether all mother wavelets work for wtc.sig", {
  for (M in MOTHERS) {
    out <- wtc.sig( quiet = TRUE, nrands = 2, mother = M,
                    lag1 = c(d1.ar1, d2.ar1),
                    dt = dt, ntimesteps = n, s0 = s0, J1 = NULL)
    expect_true( is.matrix(out) )
  }
})
