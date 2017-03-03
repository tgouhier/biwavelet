#' Power spectrum of a random red noise process
#'
#' Generate the power spectrum of a random time series with a specific AR(1)
#' coefficient.
#'
#' @param ar1 First order coefficient desired.
#' @param periods Periods of the time series at which the spectrum should be
#'   computed.
#' @return Returns the power spectrum as a vector of real numbers.
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
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#' Code based on WTC MATLAB package written by Aslak Grinsted.
#'
#' @examples
#' p <- ar1.spectrum(0.5, 1:25)
#'
#' @export
ar1.spectrum <- function(ar1, periods) {
  (1 - ar1 ^ 2) / abs(1 - ar1 * exp(-2 * 1i * pi * (1 / periods))) ^ 2
}
