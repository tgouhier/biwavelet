#' This function calls the C++ implementation of row quantile.
#'
#' @param data Numeric matrix whose row quantiles are wanted.
#' @param p Probability with value in [0,1]
#' @return A vector of length \code{nrows(data)}, where each element represents
#'   row quantile.
row_quantile <- function(data, p) {

  # fail fast
  if (length(p) > 1) {
    stop("This function only calculates one quantile at a time, for now.");
  }

  data <- as.matrix(data)
  nrow <- nrow(data);

  .C("rowQuantileC", data = as.double(data),
    nrow = nrow, ncol = ncol(data), p = as.double(p),
    quantiles = rep(0, nrow), NAOK = TRUE
  )$quantiles
}
