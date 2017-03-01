#include <Rcpp.h>
using namespace Rcpp;

const Rcomplex POW_1i[4] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
const double SQRT_ONE_HALF = sqrt(.5); // NOTE: sqrt(1/2) = 1/sqrt(2)

// Used the following R code to generate:
// options(digits = 22)
// cat(sapply(0:10, function(m) 2*pi / sqrt(m + 0.5)), sep = ",\n")
const double map_m_to_ffact[] = {
  8.885765876316732203577,
  5.130199320647456318056,
  3.973835306318440174778,
  3.358503816725427970624,
  2.961921958772244511948,
  2.679159216971467749602,
  2.464468037602168148936,
  2.29429488381819046694,
  2.155114780738997648513,
  2.038534499517198561591,
  1.939033082660811313502
};

//' Optimized "wt.bases.dog" function.
//'
//' This is a C++ version optimized for speed.
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
//' @note This c++ implementation is approx. 50% faster than the original R code
// [[Rcpp::export]]
List rcpp_wt_bases_dog(const NumericVector k,
                       const double scale,
                       const int param = -1) {

  const int m = param == -1 ? 2 : param;
  if(m < 0 || m > 10) {
    stop("Parameter 'm' must be within 0..10");
    return List::create();
  }

  const NumericVector expnt = -0.5 * pow(scale * k, 2);
  const int klen = k.length();

  // using precomputed table
  const double ffact = map_m_to_ffact[m];

  ComplexVector daughter;
  if(klen < 2) {
    daughter = NA_REAL; // becomes NA_complex_ in R
  } else {
    // R: -norm * (1i ^ m) * ((scale * k) ^ m) * exp(expnt)
    const NumericVector daughter_real =
      -sqrt(klen * scale * k[1] / R::gammafn(m + 0.5))
      * pow(scale * k, m) * exp(expnt);

    daughter = POW_1i[m % 4] * as<ComplexVector>(daughter_real);
  }

  return List::create(
    _["daughter"] = daughter,
    _["fourier.factor"] = ffact,
    _["coi"] = ffact * SQRT_ONE_HALF,
    _["dof"] = 1
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
  wt.bases.dog(k, s, p),
  rcpp_wt_bases_dog(k, s, p),
  times = 100000L
)
options(microbenchmark.unit = "t")
print(out)
options(microbenchmark.unit = "relative")
print(out)
autoplot(out)
*/
