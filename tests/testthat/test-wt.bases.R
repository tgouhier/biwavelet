context("wt.bases")
test_that("List of supported mother wavelets as global promise", {

  b <- wt.bases("morlet", c(1), 1)

  expect_true(is.list(b))
  expect_equal(b$daughter, NA_real_)
  expect_equal(round(b$coi, 2), .73)
})
