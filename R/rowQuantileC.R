#' This function calls the C++ implementation of row quantile.
#' @param data Numeric matrix whose row quantiles are wanted.
#' @param p Probability with value in [0,1]
#' @return A vector of length nrows(data), where each element represents row
#'   quantile.
rowQuantileC <- function(data, p) {

  data <- as.matrix(data)
  ncol <- ncol(data);
  nrow <- nrow(data);
  quantiles <- rep(0, nrow);

  p <- as.numeric(as.character(p));
  if (length(p) > 1) {
    stop("This function only calculates one quantile at a time, for now.");
  }

  if ( (p < 0) || (p > 1) ) {
    stop(paste("Probability", p,
               "is out of the allowed range between 0 and 1."));
  }

  res <- .C("rowQuantileC", data = as.double(data),
            nrow = as.integer(nrow), ncol = as.integer(ncol),
            p = as.double(p), quantiles = as.double(quantiles), NAOK = TRUE)

  res$quantiles;
}
