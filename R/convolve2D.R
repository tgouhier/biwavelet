#' Fast column-wise convolution of a matrix
#' 
#' Use the Fast Fourier Transform to perform convolutions between a sequence and
#' each column of a matrix.
#' 
#' @author Brandon Whitcher
#' @note This function was copied from \code{waveslim} to limit package dependencies.
#' @details This is a corrupted version of convolve made by replacing \code{fft} with
#'   \code{mvfft} in a few places. It would be nice to submit this to the R Developers
#'   for inclusion.
#'   
#' @param x M \code{x} n matrix.
#' @param y numeric sequence of length N.
#' @param conj logical; if \code{TRUE}, take the complex conjugate before
#'   back-transforming. Default is \code{TRUE} and used for usual convolution.
#'   
#' @param type character; one of \code{circular}, \code{open} (beginning of word
#'   is ok).
#'   
#'   For \code{circular}, the two sequences are treated as circular, i.e.,
#'   periodic.
#'   
#'   For \code{open} and \code{filter}, the sequences are padded with zeros
#'   (from left and right) first; \code{filter} returns the middle sub-vector of
#'   open, namely, the result of running a weighted mean of \code{x} with
#'   weights \code{y}.
#'   
#' @return M \code{x} n matrix
#' @importFrom stats fft mvfft
convolve2D <- function(x, y, conj = TRUE, type = c("circular", "open")) {
  type <- match.arg(type)
  n <- nrow(x)
  ny <- length(y)
  Real <- is.numeric(x) && is.numeric(y)
  if (type == "circular") {
    if (ny != n) {
      stop("length mismatch in convolution")
    }
  }
  else {
    n1 <- ny - 1
    x <- rbind(matrix(0, n1, ncol(x)), x)
    y <- c(y, rep.int(0, n - 1))
    n <- length(y)
  }
  
  x <- mvfft(mvfft(x) * (if (conj) Conj(fft(y)) else fft(y)),
             inverse = TRUE)
  
  (if (Real) Re(x) else x) / n
}
