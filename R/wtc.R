#' Compute wavelet coherence
#'
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#'
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
#'   \code{morlet}, \code{dog}, or \code{paul}.
#' @param param Nondimensional parameter specific to the wavelet function.
#' @param lag1 Vector containing the AR(1) coefficient of each time series.
#' @param sig.level Significance level.
#' @param sig.test Type of significance test. If set to 0, use a regular
#'   \eqn{\chi^2} test. If set to 1, then perform a time-average test. If set to
#'   2, then do a scale-average test.
#' @param nrands Number of Monte Carlo randomizations.
#' @param quiet Do not display progress bar.
#'
#' @return Return a \code{biwavelet} object containing:
#' \item{coi}{matrix containg cone of influence}
#' \item{wave}{matrix containing the cross-wavelet transform}
#' \item{wave.corr}{matrix containing the bias-corrected cross-wavelet transform
#'                  using the method described by \code{Veleda et al. (2012)}}
#' \item{power}{matrix of power}
#' \item{power.corr}{matrix of bias-corrected cross-wavelet power using the method described
#'                   by \code{Veleda et al. (2012)}}
#' \item{rsq}{matrix of wavelet coherence}
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
#' \item{type}{type of \code{biwavelet} object created (\code{\link{wtc}})}
#' \item{signif}{matrix containing \code{sig.level} percentiles of wavelet
#' coherence based on the Monte Carlo AR(1) time series}
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
#' @note The Monte Carlo randomizations can be extremely slow for large
#'   datasets. For instance, 1000 randomizations of a dataset consisting of 1000
#'   samples will take ~30 minutes on a 2.66 GHz dual-core Xeon processor.
#'
#' @examples
#' t1 <- cbind(1:100, rnorm(100))
#' t2 <- cbind(1:100, rnorm(100))
#'
#' ## Wavelet coherence
#' wtc.t1t2 <- wtc(t1, t2, nrands = 10)
#'
#' ## Plot wavelet coherence and phase difference (arrows)
#' ## Make room to the right for the color bar
#' par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
#' plot(wtc.t1t2, plot.cb = TRUE, plot.phase = TRUE)
#'
#' @export
wtc <- function(d1, d2, pad = TRUE, dj = 1 / 12,
                s0 = 2 * dt, # dt will be evaluated later (s0 is a promise)
                J1 = NULL, max.scale = NULL, mother = "morlet",
                param = -1, lag1 = NULL, sig.level = 0.95,
                sig.test = 0, nrands = 300, quiet = FALSE) {

  mother <- match.arg(tolower(mother), MOTHERS)

  # Check data format
  checked <- check.data(y = d1, x1 = d2)
  xaxis <- d1[, 1]
  dt <- checked$y$dt

  t <- checked$y$t
  n <- checked$y$n.obs

  if (is.null(J1)) {
    if (is.null(max.scale)) {
      max.scale <- (n * 0.17) * 2 * dt # automatic maxscale
    }
    J1 <- round(log2(max.scale / s0) / dj)
  }

  if (is.null(lag1)) {
    # Get AR(1) coefficients for each time series
    d1.ar1 <- arima(d1[,2], order = c(1, 0, 0))$coef[1]
    d2.ar1 <- arima(d2[,2], order = c(1, 0, 0))$coef[1]
    lag1 <- c(d1.ar1, d2.ar1)
  }

  # Get CWT of each time series
  wt1 <- wt(d = d1, pad = pad, dj = dj, s0 = s0, J1 = J1,
            max.scale = max.scale, mother = mother, param = param,
            sig.level = sig.level, sig.test = sig.test, lag1 = lag1[1])

  wt2 <- wt(d = d2, pad = pad, dj = dj, s0 = s0, J1 = J1,
            max.scale = max.scale, mother = mother, param = param,
            sig.level = sig.level, sig.test = sig.test, lag1 = lag1[2])

  # Standard deviation for each time series
  d1.sigma <- sd(d1[,2], na.rm = T)
  d2.sigma <- sd(d2[,2], na.rm = T)

  s.inv <- 1 / t(wt1$scale)
  s.inv <- matrix(rep(s.inv, n), nrow = NROW(wt1$wave))

  smooth.wt1 <- smooth.wavelet(s.inv * (abs(wt1$wave) ^ 2), dt, dj, wt1$scale)
  smooth.wt2 <- smooth.wavelet(s.inv * (abs(wt2$wave) ^ 2), dt, dj, wt2$scale)
  coi <- pmin(wt1$coi, wt2$coi, na.rm = T)

  # Cross-wavelet computation
  CW <- wt1$wave * Conj(wt2$wave)

  # Bias-corrected cross-wavelet
  CW.corr <- (wt1$wave * Conj(wt2$wave) * max(wt1$period)) /
              matrix(rep(wt1$period, length(t)), nrow = NROW(wt1$period))

  # Power
  power <- abs(CW) ^ 2

  # Bias-corrected power
  power.corr <- (abs(CW) ^ 2 * max.scale) /
                matrix(rep(wt1$period, length(t)), nrow = NROW(wt1$period))

  # Wavelet coherence
  smooth.CW <- smooth.wavelet(s.inv * (CW), dt, dj, wt1$scale)
  rsq <- abs(smooth.CW) ^ 2 / (smooth.wt1 * smooth.wt2)

  # Phase difference
  phase <- atan2(Im(CW), Re(CW))
  if (nrands > 0) {
    signif <- wtc.sig(nrands = nrands, lag1 = lag1,
                      dt = dt, ntimesteps = n, pad = pad, dj = dj, J1 = J1,
                      s0 = s0, max.scale = max.scale, mother = mother,
                      sig.level = sig.level, quiet = quiet)
  } else {
    signif <- NA
  }

  results <- list(coi = coi,
                  wave = CW,
                  wave.corr = CW.corr,
                  power = power,
                  power.corr = power.corr,
                  rsq = rsq,
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
                  type = "wtc",
                  signif = signif)

  class(results) <- "biwavelet"
  return(results)
}
