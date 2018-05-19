#' Compute partial wavelet coherence
#'
#' Compute partial wavelet coherence between \code{y} and \code{x1} by
#' partialling out the effect of \code{x2}
#'
#' @param y Time series \code{y} in matrix format (\code{n} rows x 2 columns).
#'   The first column should contain the time steps and the second column should
#'   contain the values.
#' @param x1 Time series \code{x1} in matrix format (\code{n} rows x 2 columns).
#'   The first column should contain the time steps and the second column should
#'   contain the values.
#' @param x2 Time series \code{x2} whose effects should be partialled out in
#'   matrix format (\code{n} rows x 2 columns). The first column should contain
#'   the time steps and the second column should contain the values.
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
#' @param nrands Number of Monte Carlo randomizations.
#' @param quiet Do not display progress bar.
#'
#' @return Return a \code{biwavelet} object containing:
#' \item{coi}{matrix containg cone of influence}
#' \item{wave}{matrix containing the cross-wavelet transform of \code{y} and
#'   \code{x1}}
#' \item{rsq}{matrix of partial wavelet coherence between \code{y} and \code{x1}
#'   (with \code{x2} partialled out)}
#' \item{phase}{matrix of phases between \code{y} and \code{x1}}
#' \item{period}{vector of periods}
#' \item{scale}{vector of scales}
#' \item{dt}{length of a time step}
#' \item{t}{vector of times}
#' \item{xaxis}{vector of values used to plot xaxis}
#' \item{s0}{smallest scale of the wavelet }
#' \item{dj}{spacing between successive scales}
#' \item{y.sigma}{standard deviation of \code{y}}
#' \item{x1.sigma}{standard deviation of \code{x1}}
#' \item{mother}{mother wavelet used}
#' \item{type}{type of \code{biwavelet} object created (\code{\link{pwtc}})}
#' \item{signif}{matrix containg \code{sig.level} percentiles of wavelet
#'   coherence based on the Monte Carlo AR(1) time series}
#'
#' @references
#' Aguiar-Conraria, L., and M. J. Soares. 2013. The Continuous Wavelet
#' Transform: moving beyond uni- and bivariate analysis. \emph{Journal of
#' Economic Surveys} In press.
#'
#' Cazelles, B., M. Chavez, D. Berteaux, F. Menard, J. O. Vik, S. Jenouvrier,
#' and N. C. Stenseth. 2008. Wavelet analysis of ecological time series.
#' \emph{Oecologia} 156:287-304.
#'
#' Grinsted, A., J. C. Moore, and S. Jevrejeva. 2004. Application of the cross
#' wavelet transform and wavelet coherence to geophysical time series.
#' \emph{Nonlinear Processes in Geophysics} 11:561-566.
#'
#' Ng, E. K. W., and J. C. L. Chan. 2012. Geophysical applications of partial
#' wavelet coherence and multiple wavelet coherence. \emph{Journal of
#' Atmospheric and Oceanic Technology} 29:1845-1853.
#'
#' Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis.
#' \emph{Bulletin of the American Meteorological Society} 79:61-78.
#'
#' Torrence, C., and P. J. Webster. 1998. The annual cycle of persistence in the
#' El Nino/Southern Oscillation. \emph{Quarterly Journal of the Royal
#' Meteorological Society} 124:1985-2004.
#'
#' @note The Monte Carlo randomizations can be extremely slow for large
#'   datasets. For instance, 1000 randomizations of a dataset consisting of 1000
#'   samples will take ~30 minutes on a 2.66 GHz dual-core Xeon processor.
#'
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#' Code based on WTC MATLAB package written by Aslak Grinsted.
#'
#' @example inst/doc/example-pwtc.R
#' @export
pwtc <- function(y, x1, x2, pad = TRUE, dj = 1 / 12, s0 = 2 * dt,
                 J1 = NULL, max.scale = NULL, mother = "morlet",
                 param = -1, lag1 = NULL, sig.level = 0.95,
                 sig.test = 0, nrands = 300, quiet = FALSE) {

  mother <- match.arg(tolower(mother), MOTHERS)

  checked <- check.data(y = y, x1 = x1, x2 = x2)
  xaxis <- y[, 1]
  dt <- checked$y$dt
  t <- checked$y$t
  n <- checked$y$n.obs

  if (is.null(J1)) {
    if (is.null(max.scale)) {
      max.scale <- (n * 0.17) * 2 * dt ## automatic maxscale
    }
    J1 <- round(log2(max.scale / s0) / dj)
  }

  # Get AR(1) coefficients for each time series
  y.ar1 <- arima(y[,2], order = c(1, 0, 0))$coef[1]
  x1.ar1 <- arima(x1[,2], order = c(1, 0, 0))$coef[1]

  # Get CWT of each time series
  wt.y <- wt(d = y, pad = pad, dj = dj, s0 = s0, J1 = J1,
             max.scale = max.scale, mother = mother, param = param,
             sig.level = sig.level, sig.test = sig.test, lag1 = lag1)

  wt.x1 <- wt(d = x1, pad = pad, dj = dj, s0 = s0, J1 = J1,
              max.scale = max.scale, mother = mother, param = param,
              sig.level = sig.level, sig.test = sig.test, lag1 = lag1)

  wt.x2 <- wt(d = x2, pad = pad, dj = dj, s0 = s0, J1 = J1,
              max.scale = max.scale, mother = mother, param = param,
              sig.level = sig.level, sig.test = sig.test, lag1 = lag1)

  # Standard deviation for each time series
  y.sigma <- sd(y[,2], na.rm = TRUE)
  x1.sigma <- sd(x1[,2], na.rm = TRUE)

  s.inv <- 1 / t(wt.y$scale)
  s.inv <- matrix(rep(s.inv, n), nrow = NROW(wt.y$wave))

  smooth.wt_y <- smooth.wavelet(
    s.inv * (abs(wt.y$wave) ^ 2), dt, dj, wt.y$scale)

  smooth.wt_x1 <- smooth.wavelet(
    s.inv * (abs(wt.x1$wave) ^ 2), dt, dj, wt.x1$scale)

  smooth.wt_x2 <- smooth.wavelet(
    s.inv * (abs(wt.x2$wave) ^ 2), dt, dj, wt.x2$scale)

  coi <- pmin(wt.y$coi, wt.x1$coi, wt.x2$coi, na.rm = T)

  # Cross-wavelet computation
  cw.yx1 <- wt.y$wave * Conj(wt.x1$wave)
  cw.yx2 <- wt.y$wave * Conj(wt.x2$wave)
  cw.x1x2 <- wt.x1$wave * Conj(wt.x2$wave)

  # Wavelet coherence
  smooth.cw_yx1 <- smooth.wavelet(s.inv * (cw.yx1), dt, dj, wt.y$scale)
  smooth.cw_yx2 <- smooth.wavelet(s.inv * (cw.yx2), dt, dj, wt.y$scale)
  smooth.cw_x1x2 <- smooth.wavelet(s.inv * (cw.x1x2), dt, dj, wt.y$scale)

  # Computing R^2
  rsq.yx1 <- abs(smooth.cw_yx1) ^ 2 / (smooth.wt_y * smooth.wt_x1)
  rsq.yx2 <- abs(smooth.cw_yx2) ^ 2 / (smooth.wt_y * smooth.wt_x2)
  rsq.x1x2 <- abs(smooth.cw_x1x2) ^ 2 / (smooth.wt_x1 * smooth.wt_x2)
  norm <- (1 - rsq.yx2) * (1 - rsq.x1x2)
  rsq <- abs(sqrt(rsq.yx1) - sqrt(rsq.yx2) * Conj(sqrt(rsq.x1x2))) ^ 2 / norm

  # Phase difference between y and x1
  phase <- atan2(Im(cw.yx1), Re(cw.yx1))
  if (nrands > 0) {
     signif <- wtc.sig(nrands = nrands, lag1 = c(y.ar1, x1.ar1),
                       dt = dt, n, pad = pad, dj = dj, J1 = J1,
                       s0 = s0, max.scale = max.scale, mother = mother,
                       sig.level = sig.level, quiet = quiet)
  } else {
    signif <- NA
  }

  results <- list(coi = coi,
                  wave = cw.yx1,
                  rsq = rsq,
                  phase = phase,
                  period = wt.y$period,
                  scale = wt.y$scale,
                  dt = dt,
                  t = t,
                  xaxis = xaxis,
                  s0 = s0,
                  dj = dj,
                  y.sigma = y.sigma,
                  x1.sigma = x1.sigma,
                  mother = mother,
                  type = "pwtc",
                  signif = signif)

  class(results) <- "biwavelet"
  return(results)
}
