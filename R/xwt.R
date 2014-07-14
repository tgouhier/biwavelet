xwt <-
function (d1, d2, pad=TRUE, dj=1/12, s0=2*dt, J1=NULL, max.scale=NULL, 
                mother=c("morlet", "paul", "dog"), param=-1, lag1=NULL, sig.level=0.95, sig.test=0) {
  
  mothers=c("morlet", "paul", "dog")
  mother=match.arg(tolower(mother), mothers)
  
  # Check data format
  checked=check.data(y=d1, x1=d2)
  xaxis = d1[, 1]
  dt=checked$y$dt
  dt.t2=checked$x1$dt
  t=checked$y$t
  n = checked$y$n.obs

  if (is.null(J1)) {
    if (is.null(max.scale)) {
      max.scale=(n * 0.17) * 2 * dt ## automaxscale
    }
    J1=round(log2(max.scale / s0) / dj)
  }
  # Get AR(1) coefficients for each time series
  d1.ar1=arima(d1[,2], order=c(1, 0, 0))$coef[1]
  d2.ar1=arima(d2[,2], order=c(1, 0, 0))$coef[1]  
  # Get CWT of each time series
  wt1=wt(d=d1, pad=pad, dj=dj, s0=s0, J1=J1, max.scale=max.scale, mother=mother,
         param=param, sig.level=sig.level, sig.test=sig.test, lag1=lag1)
  wt2=wt(d=d2, pad=pad, dj=dj, s0=s0, J1=J1, max.scale=max.scale, mother=mother,
         param=param, sig.level=sig.level, sig.test=sig.test, lag1=lag1)
  d1.sigma=sd(d1[,2], na.rm=T)
  d2.sigma=sd(d2[,2], na.rm=T)  
  coi=pmin(wt1$coi, wt2$coi, na.rm=T)
  # Cross-wavelet
  W.d1d2=wt1$wave*Conj(wt2$wave)
  ## Power
  power = abs(W.d1d2)
  # Bias-corrected cross-wavelet
  W.d1d2.corr=(wt1$wave*Conj(wt2$wave)*max(wt1$period))/matrix(rep(wt1$period, length(t)), nrow=NROW(wt1$period))
  # Bias-corrected power
  power.corr = abs(W.d1d2.corr)
  # Phase difference
  phase=atan2(Im(W.d1d2), Re(W.d1d2))
  # Generate two null time series with the same AR(1) coefficient as observed data
  P1=ar1.spectrum(d1.ar1, wt1$period/dt)
  P2=ar1.spectrum(d2.ar1, wt2$period/dt)
  # Significance
  if (mother=='morlet') {
    V=2
    Zv=3.9999
    signif=d1.sigma * d2.sigma * sqrt(P1 * P2) * Zv / V
    signif=matrix(signif, nrow=length(signif), ncol=1) %*% rep(1, n)
    signif=abs(W.d1d2) / signif
  }
  else
    signif=NA
  
  results=list(coi=coi, 
               wave=W.d1d2,
               wave.corr=W.d1d2.corr,
               power=power,
               power.corr=power.corr,
               phase=phase,
               period=wt1$period, 
               scale=wt1$scale, 
               dt=dt,
               t=t,
               xaxis=xaxis,
               s0=s0, 
               dj=dj, 
               d1.sigma=d1.sigma,
               d2.sigma=d2.sigma,
               mother=mother,
               type = "xwt",
               signif=signif)
  class(results)="biwavelet"
  return (results)
}

