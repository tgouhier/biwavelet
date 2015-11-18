context("Compute wavelet (wt.bases)")

test_that("List of supported mother wavelets as global promise", {

  b <- wt.bases("morlet", c(1), 1)

  expect_true(is.list(b))
  expect_equal(b$daughter, NA_real_)
  expect_equal(round(b$coi, 2), .73)
})

test_that("Unsupported mother wavelet should produce an error", {
  expect_error(wt.bases("dummy", c(1), 1))
})

test_that("Supported mothers should work without errors", {
  for (m in MOTHERS) {
    expect_true(is.list(wt.bases(m, c(1), 1)))
  }
})

test_that("When length(k) < 2 then daughter should be NA", {

  # dispatched
  expect_equal(wt.bases("morlet", c(1), 1)$daughter, NA_real_)
  expect_equal(wt.bases("paul", c(1), 1)$daughter, NA_real_)
  expect_equal(wt.bases("dog", c(1), 1)$daughter, NA_complex_)

  # original R implementation
  expect_equal(wt.bases.morlet(c(1), 1)$daughter, NA_real_)
  expect_equal(wt.bases.paul(c(1), 1)$daughter, NA_real_)
  expect_equal(wt.bases.dog(c(1), 1)$daughter, NA_complex_)

  # improved C++ implementation
  expect_equal(rcpp_wt_bases_morlet(c(1), 1)$daughter, NA_real_)
  expect_equal(rcpp_wt_bases_paul(c(1), 1)$daughter, NA_real_)
  expect_equal(rcpp_wt_bases_dog(c(1), 1)$daughter, NA_complex_)
})
