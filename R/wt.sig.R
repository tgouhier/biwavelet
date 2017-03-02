#' Determine significance of wavelet transform
#'
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#'
#' Code based on wavelet MATLAB program written by Christopher Torrence
#' and Gibert P. Compo.
#'
#' @param d Time series in matrix format (\code{n} rows x 2 columns). The first
#'   column should contain the time steps and the second column should contain
#'   the values.
#' @param dt Length of a time step.
#' @param scale The wavelet scale.
#' @param sig.test Type of significance test. If set to 0, use a regular
#'   \eqn{\chi^2} test. If set to 1, then perform a time-average test. If set to
#'   2, then do a scale-average test.
#' @param sig.level Significance level.
#' @param dof Degrees of freedom for each point in wavelet power.
#' @param lag1 AR(1) coefficient of time series used to test for significant
#'   patterns.
#' @param mother Type of mother wavelet function to use. Can be set to
#'   \code{morlet}, \code{dog}, or \code{paul}.
#' @param param Nondimensional parameter specific to the wavelet function.
#' @param sigma2 Variance of time series
#' @param arima.method Fitting method. This parameter is passed as the
#'   \code{method} Parameter to the \code{\link{arima}} function.
#'
#' @return Returns a list containing:
#' \item{signif}{vector containing significance level for each scale}
#' \item{signif}{vector of red-noise spectrum for each period}
#'
#' @references
#' Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis.
#' \emph{Bulletin of the American Meteorological Society} 79:61-78.

#' @examples
#' # Not run: wt.sig(d, dt, scale, sig.test, sig.level, lag1,
#' #                 dof = -1, mother = "morlet", sigma2 = 1)
#'
#' @export
wt.sig <- function(d, dt, scale, sig.test = 0, sig.level = 0.95,
                   dof = 2, lag1 = NULL, mother = "morlet",
                   param = -1, sigma2 = NULL, arima.method = "CSS") {

  mother <- match.arg(tolower(mother), MOTHERS)

  x <- d[, 2]

  # Find the AR1 coefficient
  if (is.null(lag1)) {
    lag1 <- arima(x, order = c(1, 0, 0), method = arima.method)$coef[1]
  }

  J1 <- length(scale) - 1
  dj <- log(scale[2] / scale[1]) / log(2)

  if (is.null(sigma2)) {
    sigma2 <- 1
  }

  # using switch instead of if..elseif..else should be faster
  # see also http://stackoverflow.com/a/7826352/855435
  switch(mother,
    morlet = {
      if (param == -1) {
        param <- 6
      }
      k0 <- param
      fourier.factor <- 4 * pi / (k0 + sqrt(2 + k0 ^ 2))
      empir <- c(2, -1, -1, -1)
      if (k0 == 6) {
        empir[2:4] <- c(0.776, 2.32, 0.60)
      }
    },

    paul = {
      if (param == -1) {
        param <- 4
      }
      m <- param
      fourier.factor <- 4 * pi / (2 * m + 1)
      empir <- c(2, -1, -1, -1)
      if (m == 4) {
        empir[2:4] <- c(1.132, 1.17, 1.5)
      }
    },

    dog = {
      if (param == -1) {
        param <- 2
      }
      m <- param
      fourier.factor <- 2 * pi * sqrt(2 / (2 * m + 1));
      empir <- c(1, -1, -1, -1)
      if (m == 2) {
        empir[2:4] <- c(3.541, 1.43, 1.4)
      }
      if (m == 6) {
        empir[2:4] <- c(1.966, 1.37, 0.97)
      }
    },

    stop("Programming error! We should never reach this code.")
  )

  # Get the appropriate parameters
  period <- scale * fourier.factor
  dofmin <- empir[1]     # Degrees of freedom with no smoothing
  Cdelta <- empir[2]     # reconstruction factor
  gamma.fac <- empir[3]  # time-decorrelation factor
  dj0 <- empir[4]        # scale-decorrelation factor
  freq <- dt / period    # normalized frequency
  fft.theor <- (1 - lag1 ^ 2) / (1 - 2 * lag1 * cos(freq * 2 * pi) + lag1 ^ 2)
  fft.theor <- sigma2 * fft.theor  # include time-series variance
  signif <- fft.theor

  if (dof[1] == -1) {
    dof <- dofmin
  }

  # using switch instead of if..elseif..else should be faster
  # see also http://stackoverflow.com/a/7826352/855435
  switch( as.character(sig.test),

    # no smoothing, DOF = dofmin
    "0" = {
      dof <- dofmin
      chisquare <- qchisq(sig.level, dof) / dof
      signif <- fft.theor * chisquare
    },

    # time-averaged significance
    "1" = {
      if (length(dof) == 1) {
        dof <- rep(dof, J1 + 1)
      }
      truncate <- which(dof < 1)
      dof[truncate] <- rep(1, length(truncate))
      dof <- dofmin * sqrt(1 + (dof * dt / gamma.fac / scale) ^ 2)
      truncate <- which(dof < dofmin)
      dof[truncate] <- dofmin * rep(1, length(truncate)) # minimum DOF is dofmin
      for (a1 in seq_len(J1 + 1)) {
        chisquare <- qchisq(sig.level, dof[a1]) / dof[a1]
        signif[a1] <- fft.theor[a1] * chisquare
      }
    },

    # scale-averaged significance
    "2" = {
      if (length(dof) != 2) {
        stop("DOF must be set to [S1,S2], the range of scale-averages")
      }
      if (Cdelta == -1) {
        stop(paste("Cdelta & dj0 not defined for", mother,
                   "with param=", param))
      }
      s1 <- dof[1]
      s2 <- dof[2]
      avg <- which((scale >= s1) & (scale <= s2)) # scales between S1 & S2
      navg <- length(avg)
      if (navg == 0) {
        stop(paste("No valid scales between", s1, "and", s2))
      }
      Savg <- 1 / sum(1 / scale[avg])
      Smid <- exp((log(s1) + log(s2)) / 2) # power-of-two midpoint
      dof <- (dofmin * navg * Savg / Smid) * sqrt(1 + (navg * dj / dj0) ^ 2)
      fft.theor <- Savg * sum(fft.theor[avg] / scale[avg])
      chisquare <- qchisq(sig.level, dof) / dof
      signif <- (dj * dt / Cdelta / Savg) * fft.theor * chisquare
    },
    # otherwise
    stop("sig.test must be 0, 1, or 2")
  )

  list(signif = signif, fft.theor = fft.theor)
}
