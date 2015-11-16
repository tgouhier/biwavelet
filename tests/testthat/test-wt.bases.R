context("wt.bases")

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
