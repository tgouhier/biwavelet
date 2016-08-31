#' @docType package
#' @name biwavelet-package
#' @aliases biwavelet
#' @useDynLib biwavelet
#' @exportPattern ^[[:alpha:]]+
#' @importFrom Rcpp evalCpp
#' @importFrom stats arima as.dist filter qchisq quantile rnorm sd ts var weighted.mean
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @title
#' Conduct Univariate and Bivariate Wavelet Analyses
#'
#' @description
#' This is a port of the WTC MATLAB package written by Aslak Grinsted and the
#' wavelet program written by Christopher Torrence and Gibert P. Compo. This
#' package can be used to perform univariate and bivariate (cross-wavelet,
#' wavelet coherence, wavelet clustering) wavelet analyses.
#'
#' @details
#' As of biwavelet version 0.14, the bias-corrected wavelet and cross-wavelet spectra
#' are automatically computed and plotted by default using the methods
#' described by Liu et al. (2007) and Veleda et al. (2012). This correction
#' is needed because the traditional approach for computing the power spectrum
#' (e.g., Torrence and Compo 1998) leads to an artificial and systematic reduction
#' in power at lower periods.
#'
#' @author Tarik C. Gouhier
#'
#' Maintainer: Tarik C. Gouhier <tarik.gouhier@@gmail.com>
#'
#' Code based on WTC MATLAB package written by Aslak Grinsted and the wavelet
#' MATLAB program written by Christopher Torrence and Gibert P. Compo.
#'
#' @references
#' Cazelles, B., M. Chavez, D. Berteaux, F. Menard, J. O. Vik, S. Jenouvrier,
#' and N. C. Stenseth. 2008. Wavelet analysis of ecological time series.
#' \emph{Oecologia} 156:287-304.
#'
#' Grinsted, A., J. C. Moore, and S. Jevrejeva. 2004. Application of the cross
#' wavelet transform and wavelet coherence to geophysical time series.
#' \emph{Nonlinear Processes in Geophysics} 11:561-566.
#'
#' Liu, Y., X. San Liang, and R. H. Weisberg. 2007. Rectification of the Bias in
#' the Wavelet Power Spectrum. \emph{Journal of Atmospheric and Oceanic Technology}
#' 24:2093-2102.
#'
#' Rouyer, T., J. M. Fromentin, F. Menard, B. Cazelles, K. Briand, R. Pianet,
#' B. Planque, and N. C. Stenseth. 2008. Complex interplays among population
#' dynamics, environmental forcing, and exploitation in fisheries.
#' \emph{Proceedings of the National Academy of Sciences} 105:5420-5425.
#'
#' Rouyer, T., J. M. Fromentin, N. C. Stenseth, and B. Cazelles. 2008.
#' Analysing multiple time series and extending significance testing in
#' wavelet analysis. \emph{Marine Ecology Progress Series} 359:11-23.
#'
#' Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis.
#' \emph{Bulletin of the American Meteorological Society} 79:61-78.
#'
#' Torrence, C., and P. J. Webster. 1998. The annual cycle of persistence in the
#' El Nino/Southern Oscillation.
#' \emph{Quarterly Journal of the Royal Meteorological Society} 124:1985-2004.
#'
#' Veleda, D., R. Montagne, and M. Araujo. 2012. Cross-Wavelet Bias Corrected by Normalizing Scales.
#' \emph{Journal of Atmospheric and Oceanic Technology} 29:1401-1408.
#'
#' @keywords wavelet
#' @keywords coherence
#' @keywords cross-wavelet
#'
#' @examples
#' # As of biwavelet version 0.14, the bias-corrected wavelet and cross-wavelet spectra
#' # are automatically computed and plotted by default using the methods
#' # described by Liu et al. (2007) and Veleda et al. (2012). This correction
#' # is needed because the traditional approach for computing the power spectrum
#' # (e.g., Torrence and Compo 1998) leads to an artificial and systematic reduction
#' # in power at low periods.
#'
#' # EXAMPLE OF BIAS CORRECTION:
#' require(biwavelet)
#' # Generate a synthetic time series 's' with the same power at three distinct periods
#' t1=sin(seq(from=0, to=2*5*pi, length=1000))
#' t2=sin(seq(from=0, to=2*15*pi, length=1000))
#' t3=sin(seq(from=0, to=2*40*pi, length=1000))
#' s=t1+t2+t3
#'
#' # Compare non-corrected vs. corrected wavelet spectrum
#' wt1=wt(cbind(1:1000, s))
#' par(mfrow=c(1,2))
#' plot(wt1, type="power.corr.norm", main="Bias-corrected")
#' plot(wt1, type="power.norm", main="Not-corrected")
#'
#' # ADDITIONAL EXAMPLES
#' t1 <- cbind(1:100, rnorm(100))
#' t2 <- cbind(1:100, rnorm(100))
#'
#' # Continuous wavelet transform
#' wt.t1 <- wt(t1)
#'
#' # Plot power
#' # Make room to the right for the color bar
#' par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
#' plot(wt.t1, plot.cb=TRUE, plot.phase=FALSE)
#'
#' # Compute cross-wavelet
#' xwt.t1t2 <- xwt(t1, t2)
#'
#' # Plot cross wavelet power and phase difference (arrows)
#' plot(xwt.t1t2, plot.cb=TRUE)
#'
#' # Wavelet coherence; nrands should be large (>= 1000)
#' wtc.t1t2=wtc(t1, t2, nrands=10)
#' # Plot wavelet coherence and phase difference (arrows)
#' # Make room to the right for the color bar
#' par(oma=c(0, 0, 0, 1), mar=c(5, 4, 4, 5) + 0.1)
#' plot(wtc.t1t2, plot.cb=TRUE)
#'
#' # Perform wavelet clustering of three time series
#' t1=cbind(1:100, sin(seq(from=0, to=10*2*pi, length.out=100)))
#' t2=cbind(1:100, sin(seq(from=0, to=10*2*pi, length.out=100)+0.1*pi))
#' t3=cbind(1:100, rnorm(100))
#'
#' # Compute wavelet spectra
#' wt.t1=wt(t1)
#' wt.t2=wt(t2)
#' wt.t3=wt(t3)
#'
#' # Store all wavelet spectra into array
#' w.arr=array(NA, dim=c(3, NROW(wt.t1$wave), NCOL(wt.t1$wave)))
#' w.arr[1, , ]=wt.t1$wave
#' w.arr[2, , ]=wt.t2$wave
#' w.arr[3, , ]=wt.t3$wave
#'
#' # Compute dissimilarity and distance matrices
#' w.arr.dis <- wclust(w.arr)
#' plot(hclust(w.arr.dis$dist.mat, method = "ward.D"), sub = "", main = "",
#'      ylab = "Dissimilarity", hang = -1)
NULL

.onAttach <- function(libname, pkgname) {

  # lazily evaluated promise of supported mother wavelets
  delayedAssign("MOTHERS", c("morlet", "paul", "dog"),
                assign.env = as.environment("package:biwavelet"))

  # just to show a startup message
  message <- paste("biwavelet", utils::packageVersion("biwavelet"), "loaded.")
  packageStartupMessage(message, appendLF = TRUE)
}

# Datasets #############################

#' Supported mother wavelets
#'
#' The list of supported mother wavelets is used in multiple places
#' therefore, we provide it as a lazily evaluated promise.
MOTHERS <- c("morlet", "paul", "dog")

#' @docType data
#' @name enviro.data
#' @title Multivariate ENSO (MEI), NPGO, and PDO indices
#' @description Monthly indices of ENSO, NPGO, and PDO from 1950 to 2009
#' @usage data(enviro.data)
#' @format A data frame with 720 observations on the following 6 variables.
#' \describe{
#'  \item{\code{month}}{a numeric vector containing the month}
#'  \item{\code{year}}{a numeric vector containing the year}
#'  \item{\code{date}}{a numeric vecor containing the date}
#'  \item{\code{mei}}{a numeric vector containing the MEI index}
#'  \item{\code{npgo}}{a numeric vector containing the NPGO index}
#'  \item{\code{pdo}}{a numeric vector containing the PDO index}
#' }
#'
#' @source
#' MEI: \url{http://www.esrl.noaa.gov/psd/enso/mei}
#' NPGO: \url{http://www.o3d.org/npgo}
#' PDO: \url{http://jisao.washington.edu/pdo}
#'
#' @references
#' Di Lorenzo, E., N. Schneider, K. M. Cobb, P. J. S. Franks, K. Chhak, A. J. Miller,
#' J. C. McWilliams, S. J. Bograd, H. Arango, E. Curchitser, T. M. Powell, and
#' P. Riviere. 2008. North Pacific Gyre Oscillation links ocean climate and
#' ecosystem change. \emph{Geophys. Res. Lett.} 35:L08607.
#'
#' Mantua, N. J., and S. R. Hare. 2002. The Pacific decadal oscillation.
#' \emph{Journal of Oceanography} 58:35-44.
#'
#' Zhang, Y., J. M. Wallace, and D. S. Battisti. 1997. ENSO-like interdecadal
#' variability: 1900-93. \emph{Journal of Climate} 10:1004-1020.
#'
#' @examples
#' data(enviro.data)
#' head(enviro.data)
#'
#' @keywords dataset
NULL
