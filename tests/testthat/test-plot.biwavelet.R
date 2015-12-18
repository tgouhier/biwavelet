context("plot.biwavelet")

# suppress generating any PDFs
pdf(NULL)

some_wt <- wt(cbind(1:720, enviro.data$mei))

test_that("Simple plotting without additional parameters should work", {
  expect_null( plot.biwavelet(some_wt) )
})

test_that("Plotting in black and white should work", {
  expect_null( plot.biwavelet(some_wt, bw = TRUE) )
})

test_that("Plotting wavelet cogerence (wtc) should work", {
  t1 <- cbind(1:100, rnorm(100))
  t2 <- cbind(1:100, rnorm(100))
  wtc.t1t2 <- wtc(t1, t2, nrands = 10, quiet = TRUE)
  expect_null( plot.biwavelet(wtc.t1t2) )
})

test_that("Setting ylim value", {
  expect_null( plot.biwavelet(some_wt, ylim = range(some_wt$period)) )
})

test_that("Plotting color bar should work", {
  expect_null( plot.biwavelet(some_wt, plot.cb = TRUE) )
})

test_that("Plotting phases should work", {
  expect_null( plot.biwavelet(some_wt, plot.phase = TRUE) )
})

test_that("Error message for unsupported plot type", {
  expect_error(
    plot.biwavelet(some_wt, type = "dummy"),
    regexp = "should be one of" )
})

test_that("Plotting all supported types should also work", {
  types <- c("power.corr.norm", "power.corr", "power.norm",
             "power", "wavelet", "phase")
  for (t in types) {
    expect_null( plot.biwavelet(some_wt, type = t) )
  }
})
