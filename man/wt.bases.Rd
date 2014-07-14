\name{wt.bases}
\alias{wt.bases}
\title{
Compute wavelet
}
\description{
Computes the wavelet as a function of Fourier frequency.
}
\usage{
wt.bases (mother = c("morlet", "paul", "dog"), k, scale, param = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mother}{
type of mother wavelet function to use. Can be set to \code{morlet}, \code{dog}, 
or \code{paul}. Default is \code{morlet}.
}
  \item{k}{
vector of frequencies at which to calculate the wavelet.
}
  \item{scale}{
the wavelet scale.
}
  \item{param}{
nondimensional parameter specific to the wavelet function.
}
}
\value{
Returns a list containing:

\item{daughter }{wavelet function}
\item{fourier.factor }{ratio of fourier period to scale}
\item{coi }{cone of influence}
\item{dof }{degrees of freedom for each point in wavelet power}
}
\references{
Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis. 
\emph{Bulletin of the American Meteorological Society} 79:61-78.
}

\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)

Code based on wavelet MATLAB program written by Christopher Torrence and
Gibert P. Compo.
}

\examples{
## Not run: wb=wt.bases(mother, k, scale[a1], param)
}
