#' Determine significance of wavelet coherence
#'
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#'
#' Code based on WTC MATLAB package written by Aslak Grinsted.
#'
#' @param nrands Number of Monte Carlo randomizations.
#' @param lag1 Vector containing the AR(1) coefficient of each time series.
#' @param dt Length of a time step.
#' @param ntimesteps Number of time steps in time series.
#' @param pad Pad the values will with zeros to increase the speed of the
#'   transform.
#' @param dj Spacing between successive scales.
#' @param s0 Smallest scale of the wavelet.
#' @param J1 Number of scales - 1.
#' @param max.scale Maximum scale.
#' @param mother Type of mother wavelet function to use. Can be set to
#'   \code{morlet}, \code{dog}, or \code{paul}.
#'   Significance testing is only available for \code{morlet} wavelet.
#' @param sig.level Significance level to compute.
#' @param quiet Do not display progress bar.
#'
#' @return Returns significance matrix containing the \code{sig.level}
#'   percentile of wavelet coherence at each time step and scale.
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
#' @note The Monte Carlo randomizations can be extremely slow for large
#'   datasets. For instance, 1000 randomizations of a dataset consisting of 1000
#'   samples will take ~30 minutes on a 2.66 GHz dual-core Xeon processor.
#'
#' @examples
#' # Not run: wtcsig <- wtc.sig(nrands, lag1 = c(d1.ar1, d2.ar1), dt,
#' #                            pad, dj, J1, s0, mother = "morlet")
#'
#' @export
wtc.sig <- function(nrands = 300, lag1, dt, ntimesteps, pad = TRUE,
                    dj = 1 / 12, s0, J1, max.scale = NULL,
                    mother = "morlet", sig.level = 0.95, quiet = FALSE) {

  if (nrands < 1) {
    return(NA)
  }

  mr1 <- get_minroots(lag1[1])
  mr2 <- get_minroots(lag1[2])
  ntseq <- seq_len(ntimesteps)

  d1 <- cbind(ntseq, ar1_ma0_sim(mr1, lag1[1], ntimesteps))

  wt1 <- wt(d = d1, pad = pad, dj = dj, dt = dt, s0 = s0, J1 = J1,
           max.scale = max.scale, mother = mother, do.sig = FALSE)

  s.inv <- 1 / t(wt1$scale)
  s.inv <- matrix(rep(s.inv, ntimesteps), nrow = NROW(wt1$wave))

  rand.rsq <- array(dim = c(NROW(wt1$wave), NCOL(wt1$wave), nrands), NA)
  if (!quiet) {
    prog.bar <- txtProgressBar(min = 0, max = nrands, style = 3)
  }

  for (r in seq_len(nrands)) {

    # Generate time series
    d1 <- cbind(ntseq, ar1_ma0_sim(mr1, lag1[1], ntimesteps))
    d2 <- cbind(ntseq, ar1_ma0_sim(mr2, lag1[2], ntimesteps))

    # Wavelet transforms
    wt1 <- wt(d = d1, pad = pad, dj = dj, dt = dt, s0 = s0, J1 = J1,
              max.scale = max.scale, mother = mother, do.sig = FALSE)
    wt2 <- wt(d = d2, pad = pad, dj = dj, dt = dt, s0 = s0, J1 = J1,
              max.scale = max.scale, mother = mother, do.sig = FALSE)

    # Smoothed cross wavelet transform
    smooth.CW <- smooth.wavelet(s.inv * wt1$wave * Conj(wt2$wave),
                                dt, dj, wt1$scale)

    sw1 <- smooth.wavelet(s.inv * (abs(wt1$wave) ^ 2), dt, dj, wt1$scale)
    sw2 <- smooth.wavelet(s.inv * (abs(wt2$wave) ^ 2), dt, dj, wt2$scale)

    rand.rsq[, , r] <- abs(smooth.CW) ^ 2 / (sw1 * sw2)

    if (!quiet) {
      setTxtProgressBar(prog.bar, r)
    }
  }

  if (!quiet) {
    close(prog.bar)
  }

  # The original slow implementation was using "apply" and "quantile" functions
  # apply(rand.rsq, MARGIN = c(1,2), quantile, sig.level, na.rm = TRUE)

  # This has been replaced with a C++ implementation taken from WGCNA package
  result <- matrix(nrow = nrow(rand.rsq), ncol = ncol(rand.rsq))
  for (i in seq_len(ncol(rand.rsq))) {
    # TODO: can be facter if we remove as.matrix()
    result[,i] <- rcpp_row_quantile(as.matrix(rand.rsq[,i,]), sig.level)
  }
  return(result)
}

#' Helper function (not exported)
#' @param ar The 'ar' part of AR(1)
#' @return double
get_minroots <- function(ar) {
  min(Mod(polyroot(c(1, -ar))))
}

#' Slightly faster \code{\link{arima.sim}} implementation which assumes AR(1)
#' and \code{ma=0}.
#'
#' @param minroots Output from \code{\link{get_minroots}} function.
#' @param ar The 'ar' part of AR(1)
#' @param n Length of output series, before un-differencing. A strictly positive
#'   integer.
#' @seealso \code{\link{arima.sim}}
ar1_ma0_sim <- function(minroots, ar, n) {

  if (minroots <= 1) {
    stop("'ar' part of model is not stationary")
  }

  nstart <- 2 + ceiling(6 / log(minroots))

  x <- ts(data = rnorm(n + nstart), start = 1 - nstart)
  x <- filter(x, ar, method = "recursive")
  x[-seq_len(nstart)]
  # maybe also this: as.ts(x)
}
