convolve2D <- function (x, y, conj = TRUE, type = c("circular", "open")) 
{
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
    n <- length(y <- c(y, rep.int(0, n - 1)))
  }
  x <- mvfft(mvfft(x) * 
    (if (conj) Conj(fft(y))
     else fft(y)), inverse = TRUE)
  (if (Real) 
    Re(x)
   else x)/n
}
