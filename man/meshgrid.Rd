\name{meshgrid}
\alias{meshgrid}

\title{
Rectangular grid in two dimensions
}
\description{
Replicates the grid vectors \code{xv} and \code{yv} to generate a full grid
}
\usage{
meshgrid (xv, yv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xv}{
vector of numeric values
}
  \item{yv}{
vector of numeric values
}
}
\value{
Returns a list containing the full grid components \code{xv} and \code{yv}:

\item{x }{replicated values of xv}
\item{y }{replicated values of yv}
}
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)
}

\examples{
xv=runif(10)
yv=runif(5)
g=meshgrid(xv, yv)
}
