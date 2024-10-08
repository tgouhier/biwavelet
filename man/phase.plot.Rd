% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phase.plot.R
\name{phase.plot}
\alias{phase.plot}
\title{Plot phases with arrows}
\usage{
phase.plot(
  x,
  y,
  phases,
  arrow.len = min(par()$pin[2]/30, par()$pin[1]/40),
  arrow.col = "black",
  arrow.lwd = arrow.len * 0.3
)
}
\arguments{
\item{x}{X-coordinates}

\item{y}{Y-coordinates}

\item{phases}{Phases}

\item{arrow.len}{Size of the arrows. Default is based on plotting region.}

\item{arrow.col}{Arrow line color.}

\item{arrow.lwd}{Width/thickness of arrows.}
}
\description{
Plot phases with arrows
}
\note{
Arrows pointing to the right mean that \code{x} and \code{y} are in phase.

Arrows pointing to the left mean that \code{x} and \code{y} are in anti-phase.

Arrows pointing up mean that \code{x} leads \code{y} by \eqn{\pi/2}.

Arrows pointing down mean that \code{y} leads \code{x} by \eqn{\pi/2}.
}
\examples{
# Code to help interpret arrow direction
a <- 0.5 * pi # phase difference
f <- 10
t <- 1:200
# x leads y by a = 0.5 * pi
x <- sin(t / max(t) * f * 2 * pi)
y <- sin(t / max(t) * f * 2 * pi - a)
par(mfrow = c(2, 1))
plot(t, x, t = "l")
lines(t, y, col = "red")
my_xwt <- xwt(cbind(t, x), cbind(t, y))
plot(my_xwt, plot.phase = TRUE)
# arrows pointing up indicating x leads y
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)

Huidong Tian provided a much better implementation of the phase.plot function
that allows for more accurate phase arrows.

Original code based on WTC MATLAB package written by Aslak Grinsted.
}
