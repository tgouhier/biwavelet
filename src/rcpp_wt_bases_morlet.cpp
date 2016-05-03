#include <Rcpp.h>
using namespace Rcpp;

const double PI4 = 4 * PI;
const double POWPI14 = pow(PI, -0.25);
const double SQRT_ONE_HALF = sqrt(.5); // NOTE: sqrt(1/2) = 1/sqrt(2)

//' Optimized "wt.bases.morlet" function.
//' 
//' This si a C++ version optimized for speed.
//' Computes the wavelet as a function of Fourier frequency
//' for "morlet" mother wavelet.
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
//' @note This c++ implementation is approx. 60% faster than the original R code
// [[Rcpp::export]]
List rcpp_wt_bases_morlet(const NumericVector k,
                       const double scale,
                       const int param = -1) {
  
  const int k0 = param == -1 ? 6 : param;
  const int klen = k.length();

  // R: expnt <- -(scale * k - k0) ^ 2 / 2 * (k > 0)
  // R: exp_expnt_kgtzero <- exp(expnt) * (k > 0)
  NumericVector exp_expnt_kgtzero(klen);
  for(int i=0; i<klen; ++i) {
    exp_expnt_kgtzero[i] = (k[i] > 0) ? exp(-pow(scale*k[i]-k0, 2) * .5) : 0;
  }
  
  NumericVector daughter;
  if(klen < 2) {
    daughter = NA_REAL; // becomes NA_real_ in R
  } else {
    // R: norm <- sqrt(scale * k[2]) * sqrt(length(k)) * (pi ^ (-1/4))
    // Note: k[2] in R is k[1] in c++ because vectors are indexed from 0
    const double norm = sqrt(scale * k[1]) * sqrt(klen) * POWPI14;
    
    // R: daughter = norm * exp(expnt) * (k > 0)
    daughter = norm * exp_expnt_kgtzero;
  }
  
  // R: fourier.factor <- 4 * pi / (k0 + sqrt(2 + k0 ^ 2))
  const double ffact = PI4 / (k0 + sqrt(2 + k0*k0));

  return List::create(
    _["daughter"] = daughter,
    _["fourier.factor"] = ffact,
    _["coi"] = ffact * SQRT_ONE_HALF,
    _["dof"] = 2
  );
}
;

/*** R
library(microbenchmark)
library(ggplot2)
library(dplyr)

k <- 1:10
s <- 2
p <- 4

out <- microbenchmark(
  wt.bases.morlet(k, s, p),
  rcpp_wt_bases_morlet(k, s, p),
  times = 2000L
)
options(microbenchmark.unit = "t")
print(out)
options(microbenchmark.unit = "relative")
print(out)
autoplot(out)
*/
