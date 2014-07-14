\name{phase.plot}
\alias{phase.plot}
\title{
Plot phases
}
\description{
Plot phases with arrows
}
\usage{
phase.plot (x, y, phases, arrow.size = 0.05, arrow.size.head = 0.05, 
            arrow.lwd = 2, arrow.col="black")
}
\arguments{
  \item{x}{
x-coordinates
}
  \item{y}{
y-coordinates
}
  \item{phases}{
phases
}
  \item{arrow.size}{
arrow size. Default is 0.05
}
  \item{arrow.size.head}{
size of the arrow head. Default is 0.05
}
  \item{arrow.lwd}{
arrow line width. Default is 2.
}
  \item{arrow.col}{
arrow line color. Default is \code{black}.
}
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)

Code based on WTC MATLAB package written by Aslak Grinsted.
}
\note{
Arrows pointing to the right mean that \code{x} and \code{y} are in phase.

Arrows pointing to the left mean that \code{x} and \code{y} are in anti-phase.

Arrows pointing up mean that \code{y} leads \code{x} by \eqn{\pi/2}.

Arrows pointing down mean that \code{x} leads \code{y} by \eqn{\pi/2}.
}
\examples{
## Not run: phase.plot(x, y, phases)
}
