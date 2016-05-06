#include "array.h"
#include <Rcpp.h>
using namespace Rcpp;

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
