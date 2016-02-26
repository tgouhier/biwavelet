library(biwavelet)

t1 <- cbind(1:100, sin(seq(0, 10 * 2 * pi, length.out = 100)))
t2 <- cbind(1:100, sin(seq(0, 10 * 2 * pi, length.out = 100) + 0.1 * pi))

# Compute wavelet spectra
wt.t1 <- wt(t1)
wt.t2 <- wt(t2)

# Compute dissimilarity
wdist(wt.t1$wave, wt.t2$wave)
