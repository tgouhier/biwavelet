wtc <-
  function (d1, d2, pad=TRUE, dj=1/12, s0=2*dt, J1=NULL, max.scale=NULL, 
            mother=c("morlet", "paul", "dog"), param=-1, lag1=NULL, 
            sig.level=0.95, sig.test = 0, nrands=300, quiet = FALSE) {
    
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
        max.scale=(n * 0.17) * 2 * dt ## automatic maxscale
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
    # Standard deviation for each time series
    d1.sigma=sd(d1[,2], na.rm=T)
    d2.sigma=sd(d2[,2], na.rm=T)  
    
    s.inv=1/t(wt1$scale)
    s.inv=matrix(rep(s.inv, n), nrow=NROW(wt1$wave))
    
    smooth.wt1=smooth.wavelet(s.inv * (abs(wt1$wave)^2), dt, dj, wt1$scale)
    smooth.wt2=smooth.wavelet(s.inv * (abs(wt2$wave)^2), dt, dj, wt2$scale)
    coi=pmin(wt1$coi, wt2$coi, na.rm=T)
    # Cross-wavelet
    CW=wt1$wave*Conj(wt2$wave)
    # Bias-corrected cross-wavelet
    CW.corr=(wt1$wave*Conj(wt2$wave)*max(wt1$period))/matrix(rep(wt1$period, length(t)), nrow=NROW(wt1$period))
    # Power
    power=abs(CW)^2
    # Bias-corrected power
    power.corr = (abs(CW)^2*max.scale)/matrix(rep(wt1$period, length(t)), nrow=NROW(wt1$period))    
    # Wavelet coherence
    smooth.CW=smooth.wavelet(s.inv * (CW), dt, dj, wt1$scale)
    rsq=abs(smooth.CW)^2 / (smooth.wt1 * smooth.wt2)
    # Phase difference
    phase=atan2(Im(CW), Re(CW))
    if (nrands > 0) {
       signif=wtc.sig(nrands = nrands, lag1 = c(d1.ar1, d2.ar1), dt = dt, n, pad = pad, 
                     dj = dj, J1 = J1, s0 = s0, max.scale = max.scale, mother = mother, 
                      sig.level = sig.level, quiet = quiet)
    }
    else
      signif=NA
    
    results=list(coi = coi, 
                 wave = CW,
                 wave.corr = CW.corr,
                 power = power,
                 power.corr = power.corr,
                 rsq = rsq, 
                 phase = phase,
                 period = wt1$period, 
                 scale = wt1$scale, 
                 dt = dt,
                 t = t,
                 xaxis = xaxis,
                 s0 = s0, 
                 dj = dj, 
                 d1.sigma = d1.sigma,
                 d2.sigma = d2.sigma,
                 mother = mother,
                 type = "wtc",
                 signif = signif)
    class(results)="biwavelet"
    return (results)
  }
