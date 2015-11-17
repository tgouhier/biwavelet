#' Compute wavelet
#' 
#' Computes the wavelet as a function of Fourier frequency.
#' 
#' @param mother type of mother wavelet function to use. Can be set to
#'   \code{morlet}, \code{dog}, or \code{paul}. Default is \code{morlet}.
#' @param k vector of frequencies at which to calculate the wavelet.
#' @param scale the wavelet scale.
#' @param param nondimensional parameter specific to the wavelet function.
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
wt.bases <- function(mother = "morlet", k, scale, param = -1) {

  mother <- match.arg(tolower(mother), MOTHERS)
  n <- length(k)

  # using switch instead of if..elseif..else should be faster
  # see also http://stackoverflow.com/a/7826352/855435
  switch(mother,
    morlet = {
      if (param == -1) {
        param <- 6
      }
      k0 <- param

      # TODO refactor this line to a nicer expression
      expnt <- -(scale * k - k0) ^ 2 / 2 * (k > 0)

      norm <- sqrt(scale * k[2]) * (pi ^ (-0.25)) * sqrt(n)
      daughter <- norm * exp(expnt)
      daughter <- daughter * (k > 0)
      fourier.factor <- 4 * pi / (k0 + sqrt(2 + k0 ^ 2))
      coi <- fourier.factor / sqrt(2)
      dof <- 2
    },

    paul = {
      if (param == -1) {
        param <- 4
      }
      m <- param
      expnt <- -(scale * k) * (k > 0)

      norm <- sqrt(scale * k[2]) *
        (2 ^ m / sqrt(m * prod(2:(2 * m - 1)))) * sqrt(n)

      daughter <- norm * ((scale * k) ^ m) * exp(expnt)
      daughter <- daughter * (k > 0)
      fourier.factor <- 4 * pi / (2 * m + 1)
      coi <- fourier.factor * sqrt(2)
      dof <- 2
    },

    dog = {
      if (param == -1) {
        param <- 2
      }
      m <- param
      expnt <- -(scale * k) ^ 2 / 2
      norm <- sqrt(scale * k[2] / gamma(m + 0.5)) * sqrt(n)
      daughter <- -norm * (1i ^ m) * ((scale * k) ^ m) * exp(expnt)
      fourier.factor <- 2 * pi * sqrt(2 / (2 * m + 1))
      coi <- fourier.factor / sqrt(2)
      dof <- 1
    },

    stop(paste("mother wavelet parameter must be one of:",
               paste(MOTHERS, collapse = ", ")))
  )

  list(daughter = daughter,
       fourier.factor = fourier.factor,
       coi = coi,
       dof = dof)
}
