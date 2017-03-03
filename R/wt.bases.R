#' Compute wavelet
#'
#' Computes the wavelet as a function of Fourier frequency.
#'
#' @param mother Type of mother wavelet function to use. Can be set to
#'   \code{morlet}, \code{dog}, or \code{paul}.
#' @param ... See parameters \code{k}, \code{scale} and \code{param} in
#'   functions: \code{\link{wt.bases.morlet}}, \code{\link{wt.bases.paul}} and
#'             \code{\link{wt.bases.dog}}
#'
#' @return Returns a list containing:
#' \item{daughter}{wavelet function}
#' \item{fourier.factor}{ratio of fourier period to scale}
#' \item{coi}{cone of influence}
#' \item{dof}{degrees of freedom for each point in wavelet power}
#'
#' @references
#' Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis.
#' \emph{Bulletin of the American Meteorological Society} 79:61-78.
#'
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#'
#' Code based on wavelet MATLAB program written by Christopher Torrence and
#' Gibert P. Compo.
#'
#' @examples
#' # Not run: wb <- wt.bases(mother, k, scale[a1], param)
#'
#' @export
wt.bases <- function(mother = "morlet", ...) {
  switch(mother,

         #morlet = wt.bases.morlet(...), # R
         morlet = rcpp_wt_bases_morlet(...), # Rcpp

         #paul = wt.bases.paul(...), # R
         paul = rcpp_wt_bases_paul(...), # Rcpp

         #dog = wt.bases.dog(...), # R
         dog = rcpp_wt_bases_dog(...), # Rcpp

         stop(paste("mother wavelet parameter must be one of:",
              paste(MOTHERS, collapse = ", ")))
  )
}

# Helper functions ########################

#' Helper method (not exported)
#'
#' @param k Vector of frequencies at which to calculate the wavelet.
#' @param scale The wavelet scale.
#' @param param Nondimensional parameter specific to the wavelet function.
#'
#' @return Returns a list containing:
#' \item{daughter}{wavelet function}
#' \item{fourier.factor}{ratio of fourier period to scale}
#' \item{coi}{cone of influence}
#' \item{dof}{degrees of freedom for each point in wavelet power}
wt.bases.morlet <- function(k, scale, param = -1) {
  k0 <- ifelse(param == -1, 6, param)
  expnt <- -(scale * k - k0) ^ 2 / 2 * (k > 0)
  norm <- sqrt(scale * k[2]) * (pi ^ (-0.25)) * sqrt(length(k))
  fourier.factor <- 4 * pi / (k0 + sqrt(2 + k0 ^ 2))
  list(daughter = norm * exp(expnt) * (k > 0),
       fourier.factor = fourier.factor,
       coi = fourier.factor / sqrt(2),
       dof = 2)
}

#' Helper method (not exported)
#'
#' @param k Vector of frequencies at which to calculate the wavelet.
#' @param scale The wavelet scale.
#' @param param Nondimensional parameter specific to the wavelet function.
#'
#' @return Returns a list containing:
#' \item{daughter}{wavelet function}
#' \item{fourier.factor}{ratio of fourier period to scale}
#' \item{coi}{cone of influence}
#' \item{dof}{degrees of freedom for each point in wavelet power}
wt.bases.paul <- function(k, scale, param = -1) {
  m <- ifelse(param == -1, 4, param)
  expnt <- -(scale * k) * (k > 0)
  norm <- sqrt(scale * k[2]) *
          (2 ^ m / sqrt(m * prod(2:(2 * m - 1)))) * sqrt(length(k))
  fourier.factor <- 4 * pi / (2 * m + 1)
  list(daughter = norm * ((scale * k) ^ m) * exp(expnt) * (k > 0),
       fourier.factor = fourier.factor,
       coi = fourier.factor / sqrt(2),
       dof = 2)
}

#' Helper method (not exported)
#'
#' @param k Vector of frequencies at which to calculate the wavelet.
#' @param scale The wavelet scale.
#' @param param Nondimensional parameter specific to the wavelet function.
#'
#' @return Returns a list containing:
#' \item{daughter}{wavelet function}
#' \item{fourier.factor}{ratio of fourier period to scale}
#' \item{coi}{cone of influence}
#' \item{dof}{degrees of freedom for each point in wavelet power}
wt.bases.dog <- function(k, scale, param = -1) {
  m <- ifelse(param == -1, 2, param)
  expnt <- -(scale * k) ^ 2 / 2
  norm <- sqrt(scale * k[2] / gamma(m + 0.5)) * sqrt(length(k))
  fourier.factor <- 2 * pi * sqrt(2 / (2 * m + 1))
  list(daughter = -norm * (1i ^ m) * ((scale * k) ^ m) * exp(expnt),
       fourier.factor = fourier.factor,
       coi = fourier.factor / sqrt(2),
       dof = 1)
}
