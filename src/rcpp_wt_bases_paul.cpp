#include <Rcpp.h>

using namespace Rcpp;

const double SQRT_ONE_HALF = sqrt(.5); // NOTE: sqrt(1/2) = 1/sqrt(2)

// Used the following R code to generate:
// options(digits = 22)
// cat(sapply(0:10, function(m) 4*pi / (2*m + 1)), sep=",\n")
const double map_m_to_ffact[] = {
  12.56637061435917246399,
  4.188790204786390525271,
  2.513274122871834492798,
  1.795195802051310351999,
  1.396263401595463582439,
  1.142397328578106607821,
  0.9666438934122439929908,
  0.8377580409572781272587,
  0.739198271432892517474,
  0.6613879270715353753118,
  0.5983986006837701543404
};

// Used the following R code to generate:
// options(scipen = 999)
// cat(sapply(0:10, function(m) prod(2:(2 * m - 1)) ), sep = ",\n")
const uint64_t map_m_to_prod[] = {
  0,                    // m = 0
  2,                    // m = 1
  6,                    // m = 2
  120,                  // m = 3
  5040,                 // m = 4
  362880,               // m = 5
  39916800,             // m = 6
  6227020800,           // m = 7
  1307674368000,        // m = 8
  355687428096000,      // m = 9
  121645100408832000    // m = 10
};

//' Optimized "wt.bases.paul" function.
//'
//' This si a C++ version optimized for speed.
//' Computes the wavelet as a function of Fourier frequency
//' for "paul" mother wavelet.
//'
//' @author Viliam Simko
//'
//' @param k vector of frequencies at which to calculate the wavelet.
//' @param scale the wavelet scale.
//' @param param nondimensional parameter specific to the wavelet function.
//' @return Returns a list containing:
//' \item{daughter}{wavelet function}
//' \item{fourier.factor}{ratio of fourier period to scale}
//' \item{coi}{cone of influence}
//' \item{dof}{degrees of freedom for each point in wavelet power}
//'
//' @note This c++ implementation is approx. 59% faster than the original R code
// [[Rcpp::export]]
List rcpp_wt_bases_paul(const NumericVector k,
                        const double scale,
                        const int param = -1) {

  const int m = (param == -1 ? 4 : param);
  if(m < 0 || m > 10) {
    stop("Parameter 'm' must be within 0..10");
    return List::create();
  }

  const int klen = k.length();

  // R: expnt <- -(scale * k) * (k > 0)
  // R: exp_expnt_kgtzero <- exp(expnt) * (k > 0)
  NumericVector exp_expnt_kgtzero(klen);
  for(int i = 0; i < klen; ++i) {
    exp_expnt_kgtzero[i] = (k[i] > 0) ? exp(-scale * k[i]) : 0;
  }

  // R: prod(2:(2 * m - 1))
  // We are actually computing factorial: (2m-1)!
  // In C++, this would be:
  // long prod = 2;
  // for(int i=3; i <= (2*m-1); ++i) {
  //   prod *= i;
  // }
  // ... but we can precompute the values in a table map_m_to_prod
  const uint64_t prod = map_m_to_prod[m]; // only works for m = 0..10
  // R: fourier.factor <- 4 * pi / (2 * m + 1)
  // .. we can precompute the values in a table
  const double ffact = map_m_to_ffact[m]; // only works for m = 0..10

  NumericVector daughter;
  if(klen < 2) {
    daughter = NA_REAL; // becomes NA_real_ in R
  } else {
    // R: sqrt(scale * k[2]) * sqrt(length(k)) * (2^m) / sqrt(m * prod(2:(2 * m - 1)))
    // Note: k[2] in R is k[1] in c++ because vectors are indexed from 0
    const double norm =
      sqrt(scale * k[1]) * sqrt((double) klen) *
      pow(2.0, m) / sqrt((double) m * prod);

    // R: daughter = norm * ((scale * k) ^ m) * exp(expnt) * (k > 0)
    daughter = norm * pow(scale * k, m) * exp_expnt_kgtzero;
  }

  return List::create(
    _["daughter"] = daughter,
    _["fourier.factor"] = ffact,
    _["coi"] = ffact * SQRT_ONE_HALF,
    _["dof"] = 2
  );
}

/*** R
library(microbenchmark)
library(ggplot2)
library(dplyr)

k <- 1:10
s <- 2
p <- 1

out <- microbenchmark(
  wt.bases.paul(k, s, p),
  rcpp_wt_bases_paul(k, s, p),
  times = 100000L
)
options(microbenchmark.unit = "t")
print(out)
options(microbenchmark.unit = "relative")
print(out)
autoplot(out)
*/
