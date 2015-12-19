# This file contains some ad-hoc code snippets
# There are some attempts to perform code profiling and benchmarking
# It would be nice to split this into multiple files and rewrite
# it in a nicer reproducible manner, including some plots and nice
# reports - "Reproducible Research"

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
  wtcsig <- wtc.sig(nrands = 400, lag1 = c(d1.ar1, d2.ar1), dt = dt,
                    ntimesteps = n,
                    pad = TRUE, dj = 1 / 12, s0 = s0,
                    J1 = NULL, max.scale = NULL, quiet = FALSE)
})

library(foreach)
library(parallel)
cl <- makeCluster(4, outfile = "") # number of cores. Notice 'outfile'
registerDoParallel(cl)

wtcsig <- wtc_sig_parallel(nrands = 4, lag1 = c(d1.ar1, d2.ar1), dt = dt,
                           ntimesteps = n,
                           pad = TRUE, dj = 1 / 12, s0 = s0,
                           J1 = NULL, max.scale = NULL)

wtcsig <- wtc_sig_parallel(nrands = 400, lag1 = c(d1.ar1, d2.ar1), dt = dt,
                  ntimesteps = n,
                  pad = TRUE, dj = 1 / 12, s0 = s0,
                  J1 = NULL, max.scale = NULL)

nrands <- 4
foreach(r = seq_len(nrands),
        .init = array(dim = c(3,3,nrands)),
        .combine = function(y, x){
          y[,,x$r] <- x$a
          return(y)
        }) %dopar% {
          list(r = r, a = array(dim = c(3,3), (1:9) * r))
        }

stopCluster(cl)

image(wtcsig)

Rprof(interval = 0.001)
wtcsig <- wtc.sig(nrands = 100, lag1 = c(d1.ar1, d2.ar1), dt = dt,
                  ntimesteps = n,
                  pad = TRUE, dj = 1 / 12, s0 = s0,
                  J1 = NULL, max.scale = NULL, quiet = FALSE)
Rprof(NULL)
summaryRprof()$by.self

library(microbenchmark)
microbenchmark(
  {seq_len(1000000)},
  {1:1000000},
  times = 1000
)



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

ntimesteps <- 100
lag1 <- c(-0.0693, 0.1581)

mr1 <- get_minroots(lag1[1])
mr2 <- get_minroots(lag1[2])
ntseq <- seq_len(ntimesteps)

microbenchmark(
  function() {
    d1  <- cbind(1:ntimesteps,
                 arima.sim(model = list(ar = lag1[1], ma = 0), n = ntimesteps))
    d2  <- cbind(1:ntimesteps,
                 arima.sim(model = list(ar = lag1[2], ma = 0), n = ntimesteps))
  },
  function() {
    d1 <- cbind(ntseq, ar1_ma0_sim(mr1, lag1[1], ntimesteps))
    d2 <- cbind(ntseq, ar1_ma0_sim(mr2, lag1[2], ntimesteps))
  },
  times = 100000
)

ts1 <- arima.sim(model = list(ar = lag1[1], ma = 0), n = ntimesteps)
ar <- lag1[1]
minroots <- get_minroots(ar)
ts2 <- ar1_ma0_sim(minroots, ar = ar, ntimesteps)
n <- ntimesteps

get_minroots <- function(ar) {
  min(Mod(polyroot(c(1, -ar))))
}

ar1_ma0_sim <- function(minroots, ar, n) {

  if (minroots <= 1) {
    stop("'ar' part of model is not stationary")
  }
  n.start <- 2 + ceiling(6 / log(minroots))

  x <- ts(data = rnorm(n + n.start), start = 1 - n.start)
  x <- filter(x, ar, method = "recursive")
  x[-seq_len(n.start)] # maybe as.ts(x)
}
ts2 <- ar1_ma0_sim(minroots, ar = ar, 10)

twave <- matrix(nrow = 100, ncol = 10, 1)
ker <- 1:3
