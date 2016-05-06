#include "array.h"
#include <Rcpp.h>
using namespace Rcpp;

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

  // fail fast
  if ((q < 0) || (q > 1)) {
    stop("value 'q' is out of range 0 to 1");
  }

  const int nr = data.nrow();
  const int nc = data.ncol();

  // here we allocate space for the result
  NumericVector result(nr);

  dArray dataArrayWrapper;
  dataArrayWrapper.wrap(data.begin(), nr, nc);

  dArray resultArrayWrapper;
  resultArrayWrapper.wrap(result.begin(), nr);

  // compute the quantiles from data to result
  dataArrayWrapper.rowQuantile(q, resultArrayWrapper);

  return result;
}

/*** R
data <- matrix(rnorm(25), 5, 5)
rcpp_row_quantile(data, .75)
*/
