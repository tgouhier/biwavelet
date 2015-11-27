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

ntimesteps <- 100L
lag1 <- c(-0.0693, 0.1581)

microbenchmark(
  function() {
    d1  <- cbind(1:ntimesteps,
                 arima.sim(model = list(ar = lag1[1], ma = 0), n = ntimesteps))
    d2  <- cbind(1:ntimesteps,
                 arima.sim(model = list(ar = lag1[2], ma = 0), n = ntimesteps))
  },
  function() {
    d1  <- cbind(1:ntimesteps,
                 arima.sim(model = list(ar = lag1[1], ma = 0), n = ntimesteps))
    d2  <- d1
    d2[,2] <- cbind(1:ntimesteps,
                    arima.sim(model = list(ar = lag1[2], ma = 0), n = ntimesteps))
  }
      
  times = 100000
)

get_minroots <- function(ar) {
  min(Mod(polyroot(c(1, -ar))))
}

ar1_ma0_sim <- function(minroots, n) {
  
  if (minroots <= 1) {
    stop("'ar' part of model is not stationary")
  }

  rand.gen <- rnorm
  innov <- rand.gen(n)
  n.start <- NA
  start.innov <- rand.gen(n.start)
  
  p <- 1 #length(model$ar)
  q <- 1 #length(model$ma)

  if (is.na(n.start)) {
    n.start <- p + q + ceiling(6 / log(minroots))
  }
  
  if (n.start < p + q)
    stop("burn-in 'n.start' must be as long as 'ar + ma'")
  
  if (!missing(start.innov) && length(start.innov) < n.start) {
    stop(sprintf(ngettext(n.start, "'start.innov' is too short: need %d point", 
                          "'start.innov' is too short: need %d points"), n.start), 
         domain = NA)
  }

  x <- ts(data = c(start.innov[seq_len(n.start)], innov[1L:n]),
          start = 1 - n.start)
  
  x <- filter(x, model$ar, method = "recursive")
  
  if (n.start > 0) {
    x <- x[-(seq_len(n.start))]
  }
  
  as.ts(x)
}



