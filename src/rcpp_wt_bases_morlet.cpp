#include <Rcpp.h>
using namespace Rcpp;

const double POWPI14 = pow(PI, -0.25);
const double SQRT_ONE_HALF = sqrt(.5); // NOTE: sqrt(1/2) = 1/sqrt(2)

// Used the following R code to generate:
// options(digits = 22)
// cat(sapply(0:10, function(m) 4 * pi / (m + sqrt(2 + m^2))), sep = ",\n")
const double map_m_to_ffact[] = {
  8.885765876316732203577,
  4.599610878225720789203,
  2.824227347583196046088,
  1.989412230649865165333,
  1.524556400231851682747,
  1.232462020317988793394,
  1.033043647749253723944,
  0.8886210242379330992435,
  0.779356281534662409527,
  0.6938746418987127295708,
  0.6252079667084186054282
};

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

  const int m = param == -1 ? 6 : param;
  if(m < 0 || m > 10) {
    stop("Parameter 'm' must be within 0..10");
    return List::create();
  }

  const int klen = k.length();

  // R: expnt <- -(scale * k - m) ^ 2 / 2 * (k > 0)
  // R: exp_expnt_kgtzero <- exp(expnt) * (k > 0)
  NumericVector exp_expnt_kgtzero(klen);
  for(int i=0; i<klen; ++i) {
    exp_expnt_kgtzero[i] = (k[i] > 0) ? exp(-pow(scale*k[i]-m, 2) * .5) : 0;
  }

  NumericVector daughter;
  if(klen < 2) {
    daughter = NA_REAL; // becomes NA_real_ in R
  } else {
    // R: norm <- sqrt(scale * k[2]) * sqrt(length(k)) * (pi ^ (-1/4))
    // Note: k[2] in R is k[1] in c++ because vectors are indexed from 0
    const double norm = sqrt(scale * k[1]) * sqrt((double) klen) * POWPI14;

    // R: daughter = norm * exp(expnt) * (k > 0)
    daughter = norm * exp_expnt_kgtzero;
  }

  // R: fourier.factor <- 4 * pi / (m + sqrt(2 + m ^ 2))
  // using precomputed table
  const double ffact = map_m_to_ffact[m];

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
p <- 4

out <- microbenchmark(
  wt.bases.morlet(k, s, p),
  rcpp_wt_bases_morlet(k, s, p),
  times = 100000L
)
options(digits = 6)
options(microbenchmark.unit = "t")
print(out)
options(microbenchmark.unit = "relative")
print(out)
autoplot(out)
*/
