check.data <- function (y, x1 = NULL, x2 = NULL) {
  y=check.datum(y)
  if (!is.null(x1)) {
    x1=check.datum(x1)
    if (y$dt != x1$dt)
      stop("The time series must have the same step size")
    if (y$n.obs != x1$n.obs)
      stop("The time series must have the same length (see merge command)")
  }
  if (!is.null(x2)) {
    x2=check.datum (x2)
    if (y$dt != x2$dt)
      stop("The time series must have the same step size")
    if (y$n.obs != x2$n.obs)
      stop("The time series must have the same length (see merge command)")
  }  
  return(list(y=y, x1=x1, x2=x2))
}

check.datum <- function (x) {
  if (NCOL(x) > 1) {
    if (class(x[, 1])[1] == "Date" | class(x[,1])[1] == "POSIXct") {
      diffs <- diff(x[, 1])
      if (all(diffs == diffs[1])) {
        t = x[, 1]
        dt = as.numeric(diffs[1])
      }
      else {
        t = 1:NROW(x)
        dt = diff(t)[1]
      }
    }
    else {
      dt = diff(x[, 1])[1]
      t = x[, 1]
    }
    epsilon <- 0.1*dt    
    if (any(abs(diff(t)-dt) > epsilon*dt))
        stop("The step size must be constant (see approx function to interpolate)")
  }
  else {
    stop("Error: the data must be in the form of an n x 2 matrix containing the 
      time steps in column 1 and the values in column 2")
  }  
  return (list(t=t, dt=dt, n.obs=NROW(x)))
}
