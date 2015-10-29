#' Check the format of time series
#' 
#' @param y time series y in matrix format (\code{n} rows x 2 columns). The
#'   first column should contain the time steps and the second column should
#'   contain the values.
#' @param x1 time series x1 in matrix format (\code{n} rows x 2 columns). The
#'   first column should contain the time steps and the second column should
#'   contain the values.
#' @param x2 time series x2 in matrix format (\code{n} rows x 2 columns). The
#'   first column should contain the time steps and the second column should
#'   contain the values.
#'   
#' @return Returns a named list containing:
#' \item{t}{time steps}
#' \item{dt}{size of a time step}
#' \item{n.obs}{number of observations}
#' 
#' @references
#' Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis. 
#' \emph{Bulletin of the American Meteorological Society} 79:61-78.
#' 
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#' 
#' @examples
#' # Not run:
#' t1=cbind(1:100, rnorm(100))
#' check.data(y=t1)
#' 
#' @export
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

#' Helper function
#' @param x TODO
#' @return TODO
#' @note This is not exported
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
