context("Significance of wavelet transform (wt.sig)")

# Setup variables =====

data <- subset(enviro.data, select = c("date", "mei"))

# Tests ==========

test_that("Testing whether all mother wavelets work for wt.sig", {
  for (M in MOTHERS) {
    out <- wt.sig(data, dt = 1, scale = 1, sig.test = 0, mother = M)
    expect_true( is.list(out) )
  }
})

test_that("Testing whether all significance tests work with morlet", {
  expect_equal(as.numeric(wt.sig(data, dt = 1, scale = 1, sig.test = 0,
                                 arima.method = "CSS-ML")$signif),
               6.301, tolerance = 0.001) # default till v0.20.14

  expect_equal(as.numeric(wt.sig(data, dt = 1, scale = 1, sig.test = 0,
                                 arima.method = "CSS")$signif),
               6.357, tolerance = 0.001) # default since v0.20.15

  expect_equal(as.numeric(wt.sig(data, dt = 1, scale = 1, sig.test = 1,
                                 arima.method = "CSS-ML")$signif),
               5.721, tolerance = 0.001) # default till v0.20.14

  expect_equal(as.numeric(wt.sig(data, dt = 1, scale = 1, sig.test = 1,
                                 arima.method = "CSS")$signif),
               5.77, tolerance = 0.001) # default since v0.20.15
})

test_that("sig.test=2 failing on other params", {
  expect_error( wt.sig(data, dt = 1, scale = 1, sig.test = 2, dof = 1),
                regexp = "DOF must be set to")

  expect_error( wt.sig(data, dt = 1, scale = 1, sig.test = 2, dof = 1:3),
                regexp = "DOF must be set to")
})
