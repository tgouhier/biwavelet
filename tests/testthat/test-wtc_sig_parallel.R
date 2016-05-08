context("Parallelized significance of wavelet coherence")

# Setup variables =====

foreach::registerDoSEQ() # to avoid warning when calling %dopar%

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

test_that("Quiet mode without progress bar should not throw errors", {
  out <- wtc_sig_parallel(quiet = TRUE, nrands = 10, lag1 = c(d1.ar1, d2.ar1),
                 dt = dt, ntimesteps = n, s0 = s0, J1 = NULL)
  expect_true( is.matrix(out) )
})

test_that("Error message for unsupported mother wavelet", {
  expect_error(
    wtc_sig_parallel( mother = "dummy",
             lag1 = c(d1.ar1, d2.ar1),
             dt = dt, ntimesteps = n, s0 = s0, J1 = NULL),
    regexp = "must be one of" )
})

test_that("Testing whether all mother wavelets are working", {
  for (M in MOTHERS) {
    out <- wtc_sig_parallel(quiet = TRUE, nrands = 2, mother = M,
                            lag1 = c(d1.ar1, d2.ar1),
                            dt = dt, ntimesteps = n, s0 = s0, J1 = NULL)
    expect_true( is.matrix(out) )
    expect_equal( dim(out), c(42, 64) )
  }
})

test_that("nrands<0 should behave nicely", {
  expect_equal(wtc_sig_parallel(
    nrands = 0,
    lag1 = 0, dt = 0, ntimesteps = 0, s0 = 0, J1 = 0), NA)
})

test_that("Progressbar should not cause errors", {
  expect_output(
    out <- wtc_sig_parallel(quiet = FALSE, nrands = 1, lag1 = c(d1.ar1, d2.ar1),
                   dt = dt, ntimesteps = n, s0 = s0, J1 = NULL),
    regexp = "\\|=+=\\| 100%"
  )
})
