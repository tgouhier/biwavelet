% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wt.sig.R
\name{wt.sig}
\alias{wt.sig}
\title{Determine significance of wavelet transform}
\usage{
wt.sig(
  d,
  dt,
  scale,
  sig.test = 0,
  sig.level = 0.95,
  dof = 2,
  lag1 = NULL,
  mother = "morlet",
  param = -1,
  sigma2 = NULL,
  arima.method = "CSS"
)
}
\arguments{
\item{d}{Time series in matrix format (\code{n} rows x 2 columns). The first
column should contain the time steps and the second column should contain
the values.}

\item{dt}{Length of a time step.}

\item{scale}{The wavelet scale.}

\item{sig.test}{Type of significance test. If set to 0, use a regular
\eqn{\chi^2} test. If set to 1, then perform a time-average test. If set to
2, then do a scale-average test.}

\item{sig.level}{Significance level.}

\item{dof}{Degrees of freedom for each point in wavelet power.}

\item{lag1}{AR(1) coefficient of time series used to test for significant
patterns.}

\item{mother}{Type of mother wavelet function to use. Can be set to
\code{morlet}, \code{dog}, or \code{paul}.}

\item{param}{Nondimensional parameter specific to the wavelet function.}

\item{sigma2}{Variance of time series}

\item{arima.method}{Fitting method. This parameter is passed as the
\code{method} Parameter to the \code{\link{arima}} function.}
}
\value{
Returns a list containing:
\item{signif}{vector containing significance level for each scale}
\item{signif}{vector of red-noise spectrum for each period}
}
\description{
Determine significance of wavelet transform
}
\examples{
# Not run: wt.sig(d, dt, scale, sig.test, sig.level, lag1,
#                 dof = -1, mother = "morlet", sigma2 = 1)

}
\references{
Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis.
\emph{Bulletin of the American Meteorological Society} 79:61-78.
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)

Code based on wavelet MATLAB program written by Christopher Torrence
and Gibert P. Compo.
}
