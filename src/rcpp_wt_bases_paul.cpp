#include <Rcpp.h>
using namespace Rcpp;

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
//' @note This c++ implementation is approx. 73% faster than the original R code
// [[Rcpp::export]]
List rcpp_wt_bases_paul(const NumericVector k,
                       const int scale,
                       const int param = -1) {
  
  const int m = param == -1 ? 4 : param;
  const int klen = k.length();

  // R: expnt <- -(scale * k) * (k > 0)
  // R: exp_expnt_kgtzero <- exp(expnt) * (k > 0)
  NumericVector exp_expnt_kgtzero(klen);
  for(int i=0; i<klen; ++i) {
    exp_expnt_kgtzero[i] = (k[i] > 0) ? exp(-scale * k[i]) : 0;
  }
  
  // R: prod(2:(2 * m - 1))
  long prod = 2;
  for(int i=3; i <= (2*m-1); ++i) {
    prod *= i;
  }
  
  // R: sqrt(scale * k[2]) * sqrt(length(k)) * (2^m) / sqrt(m * prod(2:(2 * m - 1)))
  // Note: k[2] in R is k[1] in c++ because vectors are indexed from 0
  const double norm =
    sqrt(scale * k[1]) * sqrt(klen) * pow(2, m) / sqrt(m * prod);

  // R: fourier.factor <- 4 * pi / (2 * m + 1)
  const double ffact = 4 * PI / (2 * m + 1);
  
  return List::create(
    // R: daughter = norm * ((scale * k) ^ m) * exp(expnt) * (k > 0)
    _["daughter"] = norm * pow(scale * k, m) * exp_expnt_kgtzero,
    _["fourier.factor"] = ffact,
    _["coi"] = ffact * sqrt(.5), // because sqrt(.5) = 1 / sqrt(2)
    _["dof"] = 2
  );
}
;

/*** R
library(biwavelet)
library(microbenchmark)

k <- 1:10
scale <- 2
param <- -1

microbenchmark(
  biwavelet:::wt.bases.paul(k, scale, param),
  rcpp_wt_bases_paul(k, scale, param),
  times = 100000
)
*/
