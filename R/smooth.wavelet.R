smooth.wavelet <-
function (wave, dt, dj, scale) {
  m=NCOL(wave)
  n=NROW(wave)
  twave=matrix(nrow=n, ncol=m, 0)
  ## zero-pad to power of 2... Speeds up fft calcs if n is large
  npad=2^ceiling(log2(m))
  
  k = seq(1, as.integer(npad/2), 1)
  k = k*(2*pi/npad)
  k = c(0, k, -k[seq(as.integer((npad-1)/2), 1, -1)])
  k2 = k^2
  snorm = scale / dt
  smooth=numeric(length = length(k2))
  for (ii in seq(1, n, 1)) {
    F=exp(- 0.5 * (snorm[ii]^2) * k2)
    wave.pad=rep(0i, times = length(F))
    wave.pad[1:m] = wave[ii,]
    smooth=fft(F*fft(wave.pad), inverse=TRUE)/npad
    twave[ii, ]=smooth[1:m]
  }    
  if (is.double(wave))
    twave=Re(twave)
  ## scale smoothing (boxcar with width of 0.6)
  dj0=0.6
  dj0steps=dj0/(dj*2)
  dj0steps.mod=dj0steps %% 1
  dj0steps.len=2*round(dj0steps)
  ker=c(dj0steps.mod, rep(length=dj0steps.len-1, 1), dj0steps.mod)
  ker=ker/(dj0steps.len-1+2*dj0steps.mod)
  keep.start=floor(length(ker)/2)+1
  swave=convolve2D(twave, rev(ker), type="o")
  swave=swave[keep.start:(keep.start+n-1),]
  
  return (swave)
}

