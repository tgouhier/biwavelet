#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
  #include "quantile.h"
}

//' Row-wise quantile of a matrix
//'
//' This is a C++ speed-optimized version. It is equivalent to R version
//' \code{quantile(data, q, na.rm = TRUE)}
//'
//' @author Viliam Simko
//'
//' @param data Numeric matrix whose row quantiles are wanted.
//' @param q Probability with value in [0,1]
//' @return A vector of length \code{nrows(data)}, where each element represents
//'   row quantile.
// [[Rcpp::export]]
NumericVector rcpp_row_quantile(NumericMatrix data, const double q) {

  const size_t rowLen = data.ncol();
  const size_t nrow = data.nrow();

  // fail fast
  if ((q < 0) || (q > 1)) {
    stop("value 'q' is out of range 0 to 1");
    return NumericVector(nrow, NA_REAL);
  }

  // a vector of NAs is returned for matrices without columns
  if (rowLen == 0) {
    return NumericVector(nrow, NA_REAL);
  }

  // here we allocate space for the result
  NumericVector result(nrow);

  // buffer for a row copy (needed by the quantile function)

  // VLA (variable size arrays) are forbidden in ISO C++,
  // the following line won't work with -Werror-Wall -pedantic flags
  // double rowData[rowLen];

  // therefore use heap allocation
  // (just don't forget to call delete before return)
  double* rowData = new double[rowLen];

  for (size_t row = 0; row < nrow; row++) {
    for (size_t col = 0; col < rowLen; col++) {
      rowData[col] = data(row, col);
    }
    result[row] = quantile(rowData, rowLen, q);
  }

  delete[] rowData;
  return result;
}

/*** R
data <- matrix(rnorm(25), 5, 5)
rcpp_row_quantile(data, .75)
sapply(1:5, function(x) quantile(data[x,], .75))
*/
