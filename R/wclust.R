#' Compute dissimilarity between multiple wavelet spectra
#' 
#' @param w.arr \code{N x p x t} array of wavelet spectra where \code{N} is the
#'   number of wavelet spectra to be compared, \code{p} is the number of periods
#'   in each wavelet spectrum and \code{t} is the number of time steps in each
#'   wavelet spectrum.
#' @return Returns a list containing:
#'   \item{diss.mat}{square dissimilarity matrix}
#'   \item{dist.mat}{(lower triangular) distance matrix}
#' 
#' @references
#' Rouyer, T., J. M. Fromentin, F. Menard, B. Cazelles, K. Briand, R. Pianet, 
#' B. Planque, and N. C. Stenseth. 2008. Complex interplays among population 
#' dynamics, environmental forcing, and exploitation in fisheries. 
#' \emph{Proceedings of the National Academy of Sciences} 105:5420-5425.
#' 
#' Rouyer, T., J. M. Fromentin, N. C. Stenseth, and B. Cazelles. 2008. 
#' Analysing multiple time series and extending significance testing in 
#' wavelet analysis. \emph{Marine Ecology Progress Series} 359:11-23. 
#' 
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#' @examples
#' t1 <- cbind(1:100, sin(seq(from = 0, to = 10*2*pi, length.out = 100)))
#' t2 <- cbind(1:100, sin(seq(from = 0, to = 10*2*pi, length.out = 100)+0.1*pi))
#' t3 <- cbind(1:100, rnorm(100))
#' 
#' ## Compute wavelet spectra
#' wt.t1 <- wt(t1)
#' wt.t2 <- wt(t2)
#' wt.t3 <- wt(t3)
#' 
#' ## Store all wavelet spectra into array
#' w.arr <- array(NA, dim = c(3, NROW(wt.t1$wave), NCOL(wt.t1$wave)))
#' w.arr[1, , ] <- wt.t1$wave
#' w.arr[2, , ] <- wt.t2$wave
#' w.arr[3, , ] <- wt.t3$wave
#' 
#' ## Compute dissimilarity and distance matrices
#' w.arr.dis <- wclust(w.arr)
#' plot(hclust(w.arr.dis$dist.mat, method = "ward"), sub = "", main = "", 
#'      ylab = "Dissimilarity", hang = -1)
#'
#' @export
wclust <- function (w.arr) {
  s=dim(w.arr)
  nW=s[1]
  dist.matrix=matrix(NA, nrow=nW, ncol=nW)
  k=1
  nWaves <- 1:nW
  prog.bar=txtProgressBar(min = 0, length(nWaves)^2, style = 3)
  for (n in nWaves) {
    for (j in nWaves) {
      dist.matrix[n,j]=wdist(w.arr[n, ,], w.arr[j, ,])
      k=k+1
      setTxtProgressBar(prog.bar, k)
    }
  }
  close(prog.bar)
  return (list(diss.mat=dist.matrix, dist.mat=as.dist(dist.matrix)))
}
