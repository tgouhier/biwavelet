library(testthat)
library(biwavelet)

# suppress generating any PDFs
pdf(NULL)

test_that("plot.biwavelet", {
  plot.biwavelet(wt(cbind(1:720, enviro.data$mei)))
})
