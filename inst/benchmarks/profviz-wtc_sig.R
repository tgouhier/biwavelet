if (!require("profvis", character.only = TRUE)) {
  devtools::install_github("rstudio/profvis")
  library(profvis)
}
library(biwavelet)

d1 <- cbind(1:100, rnorm(100))
d2 <- cbind(1:100, rnorm(100))

d1.ar1 <- arima(d1[,2], order = c(1, 0, 0))$coef[1]
d2.ar1 <- arima(d2[,2], order = c(1, 0, 0))$coef[1]

checked <- check.data(y = d1, x1 = d2)
dt <- checked$y$dt
n <- checked$y$n.obs
s0 <- 2 * dt

profvis({
  wtcsig <- wtc.sig(nrands = 100, lag1 = c(d1.ar1, d2.ar1), dt = dt,
                    ntimesteps = n,
                    pad = TRUE, dj = 1 / 12, s0 = s0,
                    J1 = NULL, max.scale = NULL, quiet = FALSE)
})

library(microbenchmark)

load("inst/benchmarks/data1.Rd")
out2 <- matrix(nrow = 10, ncol = 10)
out3 <- matrix(nrow = 10, ncol = 10)
microbenchmark(
  function() {
    out1 <- apply(ts, MARGIN = c(1,2), quantile, 0.95, na.rm = TRUE)
  },
  function() {
    for (i in 1:10) {
      out2[i,] <- row_quantile(ts[i,,], 0.95)
    }
  },
  function() {
    # this version is the best
    for (i in 1:10) {
      out3[,i] <- row_quantile(ts[,i,], 0.95)
    }
  },
  function() {
    out4 <- apply(ts, 2, row_quantile, 0.95)
  },
  times = 1000
)
