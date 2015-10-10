wtc.sig <-
function (nrands=300, lag1, dt, ntimesteps, pad=TRUE, dj=1/12, s0, J1, max.scale=NULL,
          mother=c("morlet", "paul", "dog"), sig.level=0.95, quiet=FALSE) {
  
  mothers=c("morlet", "paul", "dog")
  mother=match.arg(tolower(mother), mothers)
  ms = s0 * (2^(J1 * dj)) / dt ## maxscale in units of samples
  n = ceiling(ms * 6)
  d1 = cbind(1:ntimesteps, arima.sim(model = list(ar = lag1[1], ma = 0), n = ntimesteps))
  # d1 = d1 - mean(d1)
  wt1 = wt(d = d1, pad = pad, dj = dj, dt = dt, s0 = s0, J1 = J1, 
           max.scale = max.scale, mother = mother, do.sig=FALSE)
  s.inv=1/t(wt1$scale)
  s.inv=matrix(rep(s.inv, ntimesteps), nrow=NROW(wt1$wave))

  if (nrands < 1) {
    wtcsig=NA
  }
  else {
    rand.rsq=array(dim=c(NROW(wt1$wave), NCOL(wt1$wave), nrands), NA)
    if (!quiet) {
      prog.bar=txtProgressBar(min = 0, max = nrands, style = 3)
    }
    for (r in 1:nrands) {
      ## Generate time series
      d1 =  cbind(1:ntimesteps, arima.sim(model = list(ar = lag1[1], ma = 0), n = ntimesteps))
      # d1 = d1 - mean(d1)
      d2 =  cbind(1:ntimesteps, arima.sim(model = list(ar = lag1[2], ma = 0), n = ntimesteps))
      # d2 = d2 - mean(d2)
      # Wavelet transforms
      wt1 = wt(d = d1, pad = pad, dj = dj, dt = dt, s0 = s0, J1 = J1, 
               max.scale = max.scale, mother = mother, do.sig=FALSE)
      wt2 = wt(d = d2, pad = pad, dj = dj, dt = dt, s0 = s0, J1 = J1, 
               max.scale = max.scale, mother = mother, do.sig=FALSE)
      # Smoothed cross wavelet transform
      smooth.CW=smooth.wavelet(s.inv * (wt1$wave * Conj(wt2$wave)), dt, dj, wt1$scale)
      rand.rsq[, , r]=abs(smooth.CW)^2/
        (smooth.wavelet(s.inv * (abs(wt1$wave)^2), dt , dj, wt1$scale) *
           smooth.wavelet(s.inv * (abs(wt2$wave)^2), dt, dj, wt2$scale))
      if (!quiet)
        setTxtProgressBar(prog.bar, r)
    }
    close(prog.bar)
    wtcsig=apply(rand.rsq, MARGIN=c(1,2), quantile, sig.level, na.rm = TRUE)
  }
  return (wtcsig)
}
