context("Significance of wavelet transform (wt.sig)")

# Setup variables =====

data <- subset(enviro.data, select = c("date", "mei"))

# Tests ==========

test_that("Testing whether all mother wavelets work for wt.sig", {
  for (M in MOTHERS) {
    out <- wt.sig(data, dt = 1, scale = 1, sig.test = 0, mother = M )
    expect_true( is.list(out) )
  }
})

test_that("Testing whether all significance tests work with morlet", {
  expect_equal(as.numeric(wt.sig(data, dt = 1, scale = 1, sig.test = 0)$signif),
               6.3, tolerance = 0.01)

  expect_equal(as.numeric(wt.sig(data, dt = 1, scale = 1, sig.test = 1)$signif),
               5.7, tolerance = 0.01)

  # TODO significance test for sig.test=2
})

test_that("sig.test=2 failing on other params", {
  expect_error( wt.sig(data, dt = 1, scale = 1, sig.test = 2, dof = 1),
                regexp = "DOF must be set to")

  expect_error( wt.sig(data, dt = 1, scale = 1, sig.test = 2, dof = 1:3),
                regexp = "DOF must be set to")
})
