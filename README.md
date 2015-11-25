
# biwavelet R package

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/biwavelet)](http://cran.r-project.org/package=biwavelet)
[![Build Status](https://travis-ci.org/vsimko/biwavelet.svg)](https://travis-ci.org/vsimko/biwavelet)
[![codecov.io](https://codecov.io/github/vsimko/biwavelet/coverage.svg?branch=master)](https://codecov.io/github/vsimko/biwavelet?branch=master)
![CRAN Downloads](http://cranlogs-dev.r-pkg.org/badges/biwavelet)

The biwavelet R package is a port of the WTC MATLAB program written by Aslak Grinsted and the wavelet program written by Christopher Torrence and Gibert P. Compo. This package can be used to perform univariate and bivariate wavelet analyses. Wavelet analyses are resolved in the time and frequency domains, and thus ideal for identifying changes over time in the contribution of each frequency (or period) of a time series.

Since version 0.14, biwavelet also plots the **bias-corrected wavelet** and **cross-wavelet power spectrum** using the methods described by **Liu et al. (2007)** and **Veleda et al. (2012)**. This correction is needed because the traditional approach for computing the power spectrum (e.g., Torrence and Compo 1998) leads to an artificial and systematic reduction in power at lower periods. To demonstrate this bias, we can construct a time series by summing three sinusoidal waves each characterized by the same power at a different period:

```{r}
t1 <- sin(seq(from = 0, to = 2 * 5 * pi, length = 1000)) 
t2 <- sin(seq(from = 0, to = 2 * 15 * pi, length = 1000)) 
t3 <- sin(seq(from = 0, to = 2 * 40 * pi, length = 1000)) 
timeseries <- t1 + t2 + t3
```

The wavelet spectrum of the time series should show peaks of identical power at each of the three dominant periods.
However, the traditional approach leads to a consistent reduction in power at low periodicities:

```{r}
wt1 <- wt(cbind(1:1000, s))
par(mfrow = c(1,2)) 
plot(wt1, type = "power.corr.norm", main = "Bias-corrected wavelet power") 
plot(wt1, type = "power.norm", main = "Biased wavelet power")
```
