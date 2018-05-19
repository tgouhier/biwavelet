#' Compute cross-wavelet
#'
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#' Code based on WTC MATLAB package written by Aslak Grinsted.
#'
#' @param d1 Time series 1 in matrix format (\code{n} rows x 2 columns). The
#'   first column should contain the time steps and the second column should
#'   contain the values.
#' @param d2 Time series 2 in matrix format (\code{n} rows x 2 columns). The
#'   first column should contain the time steps and the second column should
#'   contain the values.
#' @param pad Pad the values will with zeros to increase the speed of the
#'   transform.
#' @param dj Spacing between successive scales.
#' @param s0 Smallest scale of the wavelet.
#' @param J1 Number of scales - 1.
#' @param max.scale Maximum scale. Computed automatically if left unspecified.
#' @param mother Type of mother wavelet function to use. Can be set to
#'   \code{morlet}, \code{dog}, or \code{paul}. Significance testing is only
#'   available for \code{morlet} wavelet.
#' @param param Nondimensional parameter specific to the wavelet function.
#' @param lag1 Vector containing the AR(1) coefficient of each time series.
#' @param sig.level Significance level.
#' @param sig.test Type of significance test. If set to 0, use a regular
#'   \eqn{\chi^2} test. If set to 1, then perform a time-average test. If set to
#'   2, then do a scale-average test.
#' @param arima.method Fitting method. This parameter is passed as the
#' \code{method} parameter to the \code{\link{arima}} function.
#'
#' @return Returns a \code{biwavelet} object containing:
#' \item{coi}{matrix containg cone of influence}
#' \item{wave}{matrix containing the cross-wavelet transform}
#' \item{wave.corr}{matrix containing the bias-corrected cross-wavelet transform
#'   using the method described by \code{Veleda et al. (2012)}}
#' \item{power}{matrix of power}
#' \item{power.corr}{matrix of bias-corrected cross-wavelet power using the
#'   method described by \code{Veleda et al. (2012)}}
#' \item{phase}{matrix of phases}
#' \item{period}{vector of periods}
#' \item{scale}{vector of scales}
#' \item{dt}{length of a time step}
#' \item{t}{vector of times}
#' \item{xaxis}{vector of values used to plot xaxis}
#' \item{s0}{smallest scale of the wavelet }
#' \item{dj}{spacing between successive scales}
#' \item{d1.sigma}{standard deviation of time series 1}
#' \item{d2.sigma}{standard deviation of time series 2}
#' \item{mother}{mother wavelet used}
#' \item{type}{type of \code{biwavelet} object created (\code{\link{xwt}})}
#' \item{signif}{matrix containg significance levels}
#'
#' @references
#' Cazelles, B., M. Chavez, D. Berteaux, F. Menard, J. O. Vik, S. Jenouvrier,
#' and N. C. Stenseth. 2008. Wavelet analysis of ecological time series.
#' \emph{Oecologia} 156:287-304.
#'
#' Grinsted, A., J. C. Moore, and S. Jevrejeva. 2004. Application of the cross
#' wavelet transform and wavelet coherence to geophysical time series.
#' \emph{Nonlinear Processes in Geophysics} 11:561-566.
#'
#' Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis.
#' \emph{Bulletin of the American Meteorological Society} 79:61-78.
#'
#' Torrence, C., and P. J. Webster. 1998. The annual cycle of persistence in the
#' El Nino/Southern Oscillation. \emph{Quarterly Journal of the Royal
#' Meteorological Society} 124:1985-2004.
#'
#' Veleda, D., R. Montagne, and M. Araujo. 2012. Cross-Wavelet Bias Corrected by
#' Normalizing Scales. \emph{Journal of Atmospheric and Oceanic Technology}
#' 29:1401-1408.
#'
#' @example inst/doc/example-xwt.R
#' @export
xwt <- function(d1, d2, pad = TRUE, dj = 1 / 12, s0 = 2 * dt,
                J1 = NULL, max.scale = NULL, mother = "morlet",
                param = -1, lag1 = NULL, sig.level = 0.95,
                sig.test = 0, arima.method = "CSS") {

  mother <- match.arg(tolower(mother), MOTHERS)

  # Check data format
  checked <- check.data(y = d1, x1 = d2)
  xaxis <- d1[, 1]
  dt <- checked$y$dt

  t <- checked$y$t
  n <- checked$y$n.obs

  if (is.null(J1)) {
    if (is.null(max.scale)) {
      max.scale <- (n * 0.17) * 2 * dt ## automaxscale
    }
    J1 <- round(log2(max.scale / s0) / dj)
  }

  # Get AR(1) coefficients for each time series
  d1.ar1 <- arima(d1[,2], order = c(1, 0, 0), method = arima.method)$coef[1]
  d2.ar1 <- arima(d2[,2], order = c(1, 0, 0), method = arima.method)$coef[1]

  # Get CWT of each time series
  wt1 <- wt(d = d1, pad = pad, dj = dj, s0 = s0, J1 = J1,
            max.scale = max.scale, mother = mother, param = param,
            sig.level = sig.level, sig.test = sig.test, lag1 = lag1)

  wt2 <- wt(d = d2, pad = pad, dj = dj, s0 = s0, J1 = J1,
            max.scale = max.scale, mother = mother, param = param,
            sig.level = sig.level, sig.test = sig.test, lag1 = lag1)

  d1.sigma <- sd(d1[,2], na.rm = T)
  d2.sigma <- sd(d2[,2], na.rm = T)
  coi <- pmin(wt1$coi, wt2$coi, na.rm = T)

  # Cross-wavelet computation
  W.d1d2 <- wt1$wave * Conj(wt2$wave)

  # Power
  power <- abs(W.d1d2)

  # Bias-corrected cross-wavelet
  W.d1d2_corr <- (wt1$wave * Conj(wt2$wave) * max(wt1$period)) /
                 matrix(rep(wt1$period, length(t)), nrow = NROW(wt1$period))

  # Bias-corrected power
  power.corr <- abs(W.d1d2_corr)

  # Phase difference
  phase <- atan2(Im(W.d1d2), Re(W.d1d2))

  # Generate two null time series with the same AR(1) coefficient
  # as observed data
  P1 <- ar1.spectrum(d1.ar1, wt1$period / dt)
  P2 <- ar1.spectrum(d2.ar1, wt2$period / dt)

  # Significance
  signif <- switch(mother,
    morlet = {
      V <- 2
      Zv <- 3.9999
      signif <- d1.sigma * d2.sigma * sqrt(P1 * P2) * Zv / V
      signif <- matrix(signif, nrow = length(signif), ncol = 1) %*% rep(1, n)
      abs(W.d1d2) / signif
    },
    NA # for other mother wavelets
  )

  results <- list(coi = coi,
                  wave = W.d1d2,
                  wave.corr = W.d1d2_corr,
                  power = power,
                  power.corr = power.corr,
                  phase = phase,
                  period = wt1$period,
                  scale = wt1$scale,
                  dt = dt,
                  t = t,
                  xaxis = xaxis,
                  s0 = s0,
                  dj = dj,
                  d1.sigma = d1.sigma,
                  d2.sigma = d2.sigma,
                  mother = mother,
                  type = "xwt",
                  signif = signif)

  class(results) <- "biwavelet"
  return(results)
}
