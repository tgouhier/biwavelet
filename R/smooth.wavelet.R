#' Smooth wavelet in both the time and scale domains
#' 
#' The time smoothing uses a filter given by the absolute value of the wavelet 
#' function at each scale, normalized to have a total weight of unity, which is
#' a Gaussian function for the Morlet wavelet. The scale smoothing is done with
#' a boxcar function of width 0.6, which corresponds to the decorrelation scale
#' of the Morlet wavelet.
#' 
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#' 
#' Code based on WTC MATLAB package written by Aslak Grinsted.
#' 
#' @references
#' Torrence, C., and P. J. Webster. 1998. The annual cycle of persistence in the
#' El Nino/Southern Oscillation.
#' \emph{Quarterly Journal of the Royal Meteorological Society} 124:1985-2004.
#' 
#' @note This function is used internally for computing wavelet coherence.
#'   It is only appropriate for the morlet wavelet.
#'   
#' @examples
#' ## Not run: smooth.wt1=smooth.wavelet(wave, dt, dj, scale)
#' 
#' @param wave wavelet coefficients
#' @param dt size of time steps
#' @param dj number of octaves per scale
#' @param scale wavelet scales
#' 
#' @return Returns the smoothed wavelet.
#' @export
#' @importFrom stats fft
smooth.wavelet <- function (wave, dt, dj, scale) {

  m <- NCOL(wave)
  n <- NROW(wave)
  
  twave <- matrix(nrow = n, ncol = m, 0)
  
  ## zero-pad to power of 2... Speeds up fft calcs if n is large
  npad <- 2^ceiling(log2(m)) # new size after padding
  
  # original slower code: k = seq(1, as.integer(npad/2), 1)
  k <- 1:as.integer(npad/2) # faster
  
  k <- k*(2*pi/npad)
  
  # original slower code: k = c(0, k, -k[seq(as.integer((npad-1)/2), 1, -1)])
  k <- c(0, k, -k[as.integer((npad-1)/2):1]) # faster
  
  k2 <- k^2
  snorm <- scale / dt
  smooth <- numeric(length = length(k2))
  for (ii in 1:n) {
    F <- exp(- 0.5 * (snorm[ii]^2) * k2)
    wave.pad <- rep(0i, times = length(F))
    wave.pad[1:m] <- wave[ii,]
    smooth <- fft(F*fft(wave.pad), inverse = TRUE) / npad
    twave[ii, ] <- smooth[1:m]
  }

  if (is.double(wave)) {
    twave <- Re(twave)
  }

  ## scale smoothing (boxcar with width of 0.6)
  dj0 <- 0.6
  dj0steps <- dj0 / (dj*2)
  dj0steps.mod <- dj0steps %% 1
  dj0steps.len <- 2*round(dj0steps)
  ker <- c(dj0steps.mod, rep(1, length = dj0steps.len - 1), dj0steps.mod)
  ker <- ker / (dj0steps.len - 1 + 2*dj0steps.mod)
  keep.start <- floor(length(ker)/2) + 1
  swave <- convolve2D(twave, rev(ker), type = "o")
  swave <- swave[keep.start:(keep.start+n-1),]
  
  return (swave)
}

