context("plot.biwavelet")

# suppress generating any PDFs
pdf(NULL)

some_wt <- wt(cbind(1:720, enviro.data$mei))

test_that("Simple plotting without additional parameters should work", {
  expect_null( plot.biwavelet(some_wt) )
})

test_that("Error message for unsupported plot type", {
  expect_error(
    plot.biwavelet(some_wt, type = "dummy"),
    regexp = "should be one of" )
})

test_that("Plotting all supported types should also work", {
  types <- c("power.corr.norm", "power.corr", "power.norm", "power", "wavelet", "phase")
  for (t in types) {
    expect_null( plot.biwavelet(some_wt, type = t) )
  }
})
