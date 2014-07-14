pwtc <-
  function (y, x1, x2, pad=TRUE, dj=1/12, s0=2*dt, J1=NULL, max.scale=NULL, 
            mother=c("morlet", "paul", "dog"), param=-1, lag1=NULL, 
            sig.level=0.95, sig.test = 0, nrands=300, quiet = FALSE) {
    
    mothers=c("morlet", "paul", "dog")
    mother=match.arg(tolower(mother), mothers)
    
    checked=check.data(y=y, x1=x1, x2=x2)
    xaxis = y[, 1]
    dt=checked$y$dt
    t=checked$y$t
    n = checked$y$n.obs
    
    if (is.null(J1)) {
      if (is.null(max.scale)) {
        max.scale=(n * 0.17) * 2 * dt ## automatic maxscale
      }
      J1=round(log2(max.scale / s0) / dj)
    }
    # Get AR(1) coefficients for each time series
    y.ar1=arima(y[,2], order=c(1, 0, 0))$coef[1]
    x1.ar1=arima(x1[,2], order=c(1, 0, 0))$coef[1]  
    x2.ar1=arima(x2[,2], order=c(1, 0, 0))$coef[1]  
    # Get CWT of each time series
    wt.y=wt(d=y, pad=pad, dj=dj, s0=s0, J1=J1, max.scale=max.scale, mother=mother,
           param=param, sig.level=sig.level, sig.test=sig.test, lag1=lag1)
    wt.x1=wt(d=x1, pad=pad, dj=dj, s0=s0, J1=J1, max.scale=max.scale, mother=mother,
           param=param, sig.level=sig.level, sig.test=sig.test, lag1=lag1)
    wt.x2=wt(d=x2, pad=pad, dj=dj, s0=s0, J1=J1, max.scale=max.scale, mother=mother,
             param=param, sig.level=sig.level, sig.test=sig.test, lag1=lag1)
    # Standard deviation for each time series
    y.sigma=sd(y[,2], na.rm=TRUE)
    x1.sigma=sd(x1[,2], na.rm=TRUE)  
    x2.sigma=sd(x2[,2], na.rm=TRUE)  
    
    s.inv=1/t(wt.y$scale)
    s.inv=matrix(rep(s.inv, n), nrow=NROW(wt.y$wave))
    
    smooth.wt.y=smooth.wavelet(s.inv * (abs(wt.y$wave)^2), dt, dj, wt.y$scale)
    smooth.wt.x1=smooth.wavelet(s.inv * (abs(wt.x1$wave)^2), dt, dj, wt.x1$scale)
    smooth.wt.x2=smooth.wavelet(s.inv * (abs(wt.x2$wave)^2), dt, dj, wt.x2$scale)
    coi=pmin(wt.y$coi, wt.x1$coi, wt.x2$coi, na.rm=T)
    
    # Cross-wavelet
    cw.yx1=wt.y$wave*Conj(wt.x1$wave)
    cw.yx2=wt.y$wave*Conj(wt.x2$wave)
    cw.x1x2=wt.x1$wave*Conj(wt.x2$wave)
    
    # Wavelet coherence
    smooth.cw.yx1=smooth.wavelet(s.inv * (cw.yx1), dt, dj, wt.y$scale)
    smooth.cw.yx2=smooth.wavelet(s.inv * (cw.yx2), dt, dj, wt.y$scale)
    smooth.cw.x1x2=smooth.wavelet(s.inv * (cw.x1x2), dt, dj, wt.y$scale)

    # R^2
    rsq.yx1=abs(smooth.cw.yx1)^2 / (smooth.wt.y * smooth.wt.x1)
    rsq.yx2=abs(smooth.cw.yx2)^2 / (smooth.wt.y * smooth.wt.x2)
    rsq.x1x2=abs(smooth.cw.x1x2)^2 / (smooth.wt.x1 * smooth.wt.x2)
    rsq=abs(sqrt(rsq.yx1) - sqrt(rsq.yx2)*Conj(sqrt(rsq.x1x2)))^2/((1-rsq.yx2) * (1-rsq.x1x2))

    # Phase difference between y and x1
    phase=atan2(Im(cw.yx1), Re(cw.yx1))
    if (nrands > 0) {
       signif=wtc.sig(nrands = nrands, lag1 = c(y.ar1, x1.ar1), dt = dt, n, pad = pad, 
                     dj = dj, J1 = J1, s0 = s0, max.scale = max.scale, mother = mother, 
                      sig.level = sig.level, quiet = quiet)
    }
    else
      signif=NA
    
    results=list(coi = coi, 
                 wave = cw.yx1,
                 rsq = rsq, 
                 phase = phase,
                 period = wt.y$period, 
                 scale = wt.y$scale, 
                 dt = dt,
                 t = t,
                 xaxis = xaxis,
                 s0 = s0, 
                 dj = dj, 
                 y.sigma = y.sigma,
                 x1.sigma = x1.sigma,
                 mother = mother,
                 type = "pwtc",
                 signif = signif)
    class(results)="biwavelet"
    return (results)
  }
