library(testthat)
library(biwavelet)
library(dplyr)

test_that("plot.biwavelet", {
  cbind(1:720, enviro.data$mei) %>% wt %>% plot.biwavelet
})