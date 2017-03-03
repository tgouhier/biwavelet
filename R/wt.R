#' Compute wavelet transform
#'
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#'
#' Code based on wavelet MATLAB program written by Christopher Torrence
#' and Gibert P. Compo.
#'
#' @param d Time series in matrix format (\code{n} rows x 2 columns). The first
#'   column should contain the time steps and the second column should contain
#'   the values.
#' @param pad Pad the values will with zeros to increase the speed of the
#'   transform.
#' @param dt Length of a time step.
#' @param dj Spacing between successive scales.
#' @param s0 Smallest scale of the wavelet.
#' @param J1 Number of scales - 1.
#' @param max.scale Maximum scale. Computed automatically if left unspecified.
#' @param mother Type of mother wavelet function to use. Can be set to
#'   \code{morlet}, \code{dog}, or \code{paul}.
#' @param param Nondimensional parameter specific to the wavelet function.
#' @param lag1 AR(1) coefficient of time series used to test for significant
#'   patterns.
#' @param sig.level Significance level.
#' @param sig.test Type of significance test. If set to 0, use a regular
#'   \eqn{\chi^2} test. If set to 1, then perform a time-average test.
#'   If set to 2, then do a scale-average test.
#' @param do.sig Perform significance testing if \code{TRUE}.
#'
#' @param arima.method Fitting method. This parameter is passed as the
#' \code{method} Parameter to the \code{\link{arima}} function.
#'
#' @return Returns a \code{biwavelet} object containing:
#' \item{coi}{matrix containg cone of influence}
#' \item{wave}{matrix containing the wavelet transform}
#' \item{power}{matrix of power}
#' \item{power.corr}{matrix of bias-corrected power using the method described
#'   by \code{Liu et al. (2007)}}
#' \item{phase}{matrix of phases}
#' \item{period}{vector of periods}
#' \item{scale}{vector of scales}
#' \item{dt}{length of a time step}
#' \item{t}{vector of times}
#' \item{xaxis}{vector of values used to plot xaxis}
#' \item{s0}{smallest scale of the wavelet }
#' \item{dj}{spacing between successive scales}
#' \item{sigma2}{variance of time series}
#' \item{mother}{mother wavelet used}
#' \item{type}{type of \code{biwavelet} object created (\code{\link{wt}})}
#' \item{signif}{matrix containg significance levels}
#'
#' @references
#' Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis.
#' \emph{Bulletin of the American Meteorological Society} 79:61-78.
#'
#' Liu, Y., X. San Liang, and R. H. Weisberg. 2007. Rectification of the Bias in
#' the Wavelet Power Spectrum. \emph{Journal of Atmospheric and Oceanic
#' Technology} 24:2093-2102.
#'
#' @examples
#' t1 <- cbind(1:100, rnorm(100))
#'
#' ## Continuous wavelet transform
#' wt.t1 <- wt(t1)
#'
#' ## Plot power
#' ## Make room to the right for the color bar
#' par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
#' plot(wt.t1, plot.cb = TRUE, plot.phase = FALSE)
#'
#' @export
wt <- function(d, pad = TRUE, dt = NULL, dj = 1 / 12, s0 = 2 * dt,
               J1 = NULL, max.scale = NULL, mother = "morlet",
               param = -1, lag1 = NULL, sig.level = 0.95,
               sig.test = 0, do.sig = TRUE, arima.method = "CSS") {

  # Check data format
  checked <- check.datum(d)
  n.obs <- checked$n.obs
  dt <- checked$dt
  t <- checked$t
  xaxis <- d[, 1]
  x <- d[, 2] - mean(d[, 2])
  sigma2 <- var(d[,2])

  if (is.null(J1)) {
    if (is.null(max.scale)) {
      max.scale <- n.obs * 0.34 * dt # automaxscale
    }
    J1 <- round(log2(max.scale / s0) / dj)
  }

  # This could be made more efficient by removing the +1
  # but this can lead to insufficient padding in some instances.
  # Currently the padding is the same as that of Torrence & Compo (1998)
  if (pad) {
    x <- c(x, rep(0, 2 ^ ceiling(log2(n.obs) + 1) - n.obs))
  }

  n <- NROW(x)
  k <- seq_len(floor(.5 * n))
  k <- k * 2 * pi / (n * dt)
  k <- c(0, k, -k[ floor( .5 * (n - 1) ):1 ])
  f <- fft(x)
  invflen <- 1 / length(f)
  scale <- s0 * 2 ^ ((0:J1) * dj)
  wave <- matrix(0, nrow = J1 + 1, ncol = n)

  for (a1 in seq_len(J1 + 1)) {
    wb <- wt.bases(mother, k, scale[a1], param)
    wave[a1, ] <- fft(f * wb$daughter, inverse = TRUE) * invflen
  }

  period <- wb$fourier.factor * scale
  coi <- wb$coi * dt * c(1e-5,
                         seq_len(.5 * (n.obs + 1) - 1),
                         floor(.5 * n.obs - 1):1,
                         1e-5)

  wave <- wave[, seq_len(n.obs)] ## Get rid of padding before returning
  power.corr <- (abs(wave) ^ 2 * max.scale) /
                matrix(rep(period, length(t)), nrow = NROW(period))

  power <- abs(wave) ^ 2

  phase <- atan2(Im(wave), Re(wave))
  if (do.sig) {
    signif <- wt.sig(d = d, dt = dt, scale = scale, sig.test = sig.test,
                     sig.level = sig.level, lag1 = lag1, dof = -1,
                     mother = mother, sigma2 = 1,
                     arima.method = arima.method)$signif

    signif <- matrix(signif, nrow = length(signif), ncol = 1) %*% rep(1, n.obs)
    signif <- power / (sigma2 * signif)
  } else {
    signif <- NA
  }

  results <- list(coi = coi,
                  wave = wave,
                  power = power,
                  power.corr = power.corr,
                  phase = phase,
                  period = period,
                  scale = scale,
                  dt = dt,
                  t = t,
                  xaxis = xaxis,
                  s0 = s0,
                  dj = dj,
                  sigma2 = sigma2,
                  mother = mother,
                  type = "wt",
                  signif = signif)

  class(results) <- "biwavelet"
  return(results)
}
