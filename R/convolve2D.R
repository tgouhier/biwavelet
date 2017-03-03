#' Fast column-wise convolution of a matrix
#'
#' Use the Fast Fourier Transform to perform convolutions between a sequence and
#' each column of a matrix.
#'
#' @author Brandon Whitcher
#' @note This function was copied from \code{waveslim} to limit package
#'   dependencies.
#' @details This is a corrupted version of convolve made by replacing
#'   \code{\link{fft}} with \code{\link{mvfft}} in a few places. It would be
#'   nice to submit this to the R Developers for inclusion.
#'
#' @param x M \code{x} n matrix.
#' @param y Numeric sequence of length N.
#' @param conj Logical; if \code{TRUE}, take the complex conjugate before
#'   back-transforming. \code{TRUE} is used for usual convolution.
#'
#' @param type Character; one of \code{circular}, \code{open} (beginning of word
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
#'
#' @export
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
  } else {
    n1 <- ny - 1
    x <- rbind(matrix(0, n1, ncol(x)), x)
    y <- c(y, rep.int(0, n - 1))
    n <- length(y)
  }

  x <- mvfft(mvfft(x) * (if (conj) Conj(fft(y)) else fft(y)),
             inverse = TRUE)

  (if (Real) Re(x) else x) / n
}

#' Speed-optimized version of convolve2D
#'
#' Equivalent to \code{convolve2D(x, y, type = "open")}. The motivation for this
#' function was that convolution is called many times in a loop from
#' \code{\link{smooth.wavelet}}, always with the \code{type = "open"} parameter.
#'
#' @author Viliam Simko
#' @inheritParams convolve2D
#' @seealso \code{\link{convolve2D}}
convolve2D_typeopen <- function(x, y) {
  n <- nrow(x)
  Real <- is.numeric(x) && is.numeric(y)

  # assuming type = "open"
  x <- rbind(matrix(0, length(y) - 1, ncol(x)), x)
  y <- c(y, rep.int(0, n - 1))
  n <- length(y)

  x <- mvfft(mvfft(x) * Conj(fft(y)), inverse = TRUE)

  (if (Real) Re(x) else x) / n
}
