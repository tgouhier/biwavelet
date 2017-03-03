#' Parallel \code{\link{wtc.sig}}
#'
#' Parallelized Monte Carlo simulation equivalent to \code{\link{wtc.sig}}.
#'
#' @examples
#' # Not run: library(foreach)
#' # library(doParallel)
#' # cl <- makeCluster(4, outfile="") # number of cores. Notice 'outfile'
#' # registerDoParallel(cl)
#' # wtc_sig_parallel(your parameters go here)
#' # stopCluster(cl)
#'
#' @seealso \code{\link[foreach]{foreach}}
#' @seealso \code{\link{wtc.sig}}
#'
#' @inheritParams wtc.sig
#' @export
#' @import foreach
wtc_sig_parallel <- function(nrands = 300, lag1, dt, ntimesteps, pad = TRUE,
                    dj = 1 / 12, s0, J1, max.scale = NULL,
                    mother = "morlet", sig.level = 0.95, quiet = TRUE) {

  if (nrands < 1) {
    return(NA)
  }

  mr1 <- get_minroots(lag1[1])
  mr2 <- get_minroots(lag1[2])
  ntseq <- seq_len(ntimesteps)

  d1 <- cbind(ntseq, ar1_ma0_sim(mr1, lag1[1], ntimesteps))

  wt1 <- wt(d = d1, pad = pad, dj = dj, dt = dt, s0 = s0, J1 = J1,
            max.scale = max.scale, mother = mother, do.sig = FALSE)

  s.inv <- 1 / t(wt1$scale)
  s.inv <- matrix(rep(s.inv, ntimesteps), nrow = NROW(wt1$wave))

  NUMC <- getDoParWorkers()
  CLEN <- ceiling(nrands / NUMC)

  if (!quiet) {
    prog.bar <- txtProgressBar(min = 0, max = NUMC, style = 3)
  }

  CID <- NULL # this is only necessary for R CMD check --as-cran
  rand.rsq <- foreach(CID = seq_len(NUMC),
                      .final = function(x) {
                        array(unlist(x),
                              dim = c(NROW(wt1$wave), NCOL(wt1$wave), nrands))
                      }) %dopar% {

    if (CID == NUMC) {
      CLEN <- nrands - CLEN * (NUMC - 1)
    }

    out <- vector("list", CLEN)
    for (r in seq_len(CLEN)) {

      # Generate time series
      d1 <- cbind(ntseq, ar1_ma0_sim(mr1, lag1[1], ntimesteps))
      d2 <- cbind(ntseq, ar1_ma0_sim(mr2, lag1[2], ntimesteps))

      # Wavelet transforms
      wt1 <- wt(d = d1, pad = pad, dj = dj, dt = dt, s0 = s0, J1 = J1,
                max.scale = max.scale, mother = mother, do.sig = FALSE)
      wt2 <- wt(d = d2, pad = pad, dj = dj, dt = dt, s0 = s0, J1 = J1,
                max.scale = max.scale, mother = mother, do.sig = FALSE)

      # Smoothed cross wavelet transform
      smooth.CW <- smooth.wavelet(s.inv * wt1$wave * Conj(wt2$wave),
                                  dt, dj, wt1$scale)

      sw1 <- smooth.wavelet(s.inv * (abs(wt1$wave) ^ 2), dt, dj, wt1$scale)
      sw2 <- smooth.wavelet(s.inv * (abs(wt2$wave) ^ 2), dt, dj, wt2$scale)

      #rand.rsq[, , r] <- abs(smooth.CW) ^ 2 / (sw1 * sw2)
      out[[r]] <- abs(smooth.CW) ^ 2 / (sw1 * sw2)
    }

    # returned to the combine function
    return(out)
  }

  # TODO: the progressbar should be implemented properly
  if (!quiet) {
    setTxtProgressBar(prog.bar, NUMC)
  }

  if (!quiet) {
    close(prog.bar)
  }

  # The original slow implementation was using "apply" and "quantile" functions
  # apply(rand.rsq, MARGIN = c(1,2), quantile, sig.level, na.rm = TRUE)
  # This has been replaced with a C++ implementation taken from WGCNA package
  j <- NULL # this is only necessary for R CMD check --as-cran
  foreach(j = seq_len(ncol(rand.rsq)), .combine = cbind) %dopar% {
    # TODO: can be facter if we remove as.matrix()
    rcpp_row_quantile(as.matrix(rand.rsq[,j,]), sig.level)
  }
}
