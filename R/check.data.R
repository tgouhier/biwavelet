#' Check the format of time series
#'
#' @param y Time series \code{y} in matrix format (\code{n} rows x 2 columns).
#'   The first column should contain the time steps and the second column should
#'   contain the values.
#' @param x1 Time series \code{x1} in matrix format (\code{n} rows x 2 columns).
#'   The first column should contain the time steps and the second column should
#'   contain the values.
#' @param x2 Time series \code{x2} in matrix format (\code{n} rows x 2 columns).
#'   The first column should contain the time steps and the second column should
#'   contain the values.
#'
#' @return Returns a named list containing:
#' \item{t}{Time steps}
#' \item{dt}{Size of a time step}
#' \item{n.obs}{Number of observations}
#'
#' @references
#' Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis.
#' \emph{Bulletin of the American Meteorological Society} 79:61-78.
#'
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#'
#' @examples
#' t1 <- cbind(1:100, rnorm(100))
#' check.data(y = t1)
#'
#' @export
check.data <- function(y, x1 = NULL, x2 = NULL) {

  y.check <- check.datum(y)
  x1.check <- NULL
  x2.check <- NULL

  if (!is.null(x1)) {
    x1.check <- check.datum(x1)
    if (any(diff(y[, 1]) != diff(x1[, 1]))) {
      stop("The time series must have the same step size")
    }
    if (y.check$n.obs != x1.check$n.obs) {
      stop("The time series must have the same length (see merge command)")
    }
  }

  if (!is.null(x2)) {
    x2.check <- check.datum(x2)
    if (any(diff(y[, 1]) != diff(x2[, 1]))) {
      stop("The time series must have the same step size")
    }
    if (y.check$n.obs != x2.check$n.obs) {
      stop("The time series must have the same length (see merge command)")
    }
  }
  return(list(y = y.check, x1 = x1.check, x2 = x2.check))
}

#' Helper function
#' @param x matrix
#' @return list(t, dt, n.obs)
#' @note This function is not exported
check.datum <- function(x) {
  if (NCOL(x) > 1) {
    t <- x[, 1]
    diffs <- diff(t)
    dt <- as.numeric(diffs[1])
    epsilon <- 0.1 * dt
    if (any(abs(diff(t) - dt) > epsilon)) {
      stop("The step size must be constant ",
           "(see approx function to interpolate)")
    } else {
      if (class(t) == "Date") {
        t <- seq_len(NROW(t))
        dt <- diff(t)[1]
      }
    }
  } else {
    stop("Error: the data must be in the form of an n x 2 matrix ",
         "containing the time steps in column 1 and the values in column 2")
  }
  return(list(t = t, dt = dt, n.obs = NROW(x)))
}
