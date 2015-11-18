#include <Rcpp.h>
using namespace Rcpp;

const Rcomplex POW_1i[4] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
const double PI2 = 2 * PI;
const double SQRT_ONE_HALF = sqrt(.5); // NOTE: sqrt(1/2) = 1/sqrt(2)

//' Optimized "wt.bases.dog" function.
//' 
//' This si a C++ version optimized for speed.
//' Computes the wavelet as a function of Fourier frequency
//' for "dog" mother wavelet.
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
//' @note This c++ implementation is approx. 36% faster than the original R code
// [[Rcpp::export]]
List rcpp_wt_bases_dog(const NumericVector k,
                       const int scale,
                       const int param = -1) {
  
  const int m = param == -1 ? 2 : param;
  const NumericVector expnt = -0.5 * pow(scale * k, 2);
  const double ffact = PI2 / sqrt(m + 0.5);

  // R: -norm * (1i ^ m) * ((scale * k) ^ m) * exp(expnt)
  const NumericVector daughter_real =
    -sqrt(k.length() * scale * k[1] / R::gammafn(m + 0.5))
    * pow(scale * k, m) * exp(expnt);
  
  return List::create(
    _["daughter"] = POW_1i[m % 4] * as<ComplexVector>(daughter_real),
    _["fourier.factor"] = ffact,
    _["coi"] = ffact * SQRT_ONE_HALF,
    _["dof"] = 1
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
  wt.bases.dog(k, scale, param),
  rcpp_wt_bases_dog(k, scale, param),
  times = 100000
)

*/
