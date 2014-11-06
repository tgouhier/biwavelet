\name{phase.plot}
\alias{phase.plot}
\title{
Plot phases
}
\description{
Plot phases with arrows
}
\usage{
phase.plot (x, y, phases, arrow.len=min(par()$pin[2]/30,par()$pin[1]/40), 
            arrow.col="black", arrow.lwd=arrow.len*0.3)
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
  \item{arrow.len}{
size of the arrows. Default is based on plotting region (min(par()$pin[2]/30,par()$pin[1]/40).
}
  \item{arrow.lwd}{
width/thickness of arrows. Default is arrow.len*0.3.
}
  \item{arrow.col}{
arrow line color. Default is \code{black}.
}
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)

Huidong Tian provided a much better implementation of the phase.plot function
that allows for better looking phase arrows.

Original code based on WTC MATLAB package written by Aslak Grinsted.
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
