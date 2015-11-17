# Tests ==========
context("Cross-wavelet (xwt)")
pdf(NULL)

test_that("Plotting simple cross-wavelet should work", {
  t1 <- cbind(1:100, rnorm(100))
  t2 <- cbind(1:100, rnorm(100))
  xwt.t1t2 <- xwt(t1, t2)
  expect_null(plot(xwt.t1t2))
})

test_that("xwt should work with sample enviro.data", {

  x <- xwt(subset(enviro.data, select = c("date", "mei")),
           subset(enviro.data, select = c("date", "npgo")))

  expect_equal(class(x), "biwavelet")
  expect_equal(x$type, "xwt")
  expect_null(plot(x))
})
