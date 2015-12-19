library(biwavelet)

t1 <- cbind(1:100, rnorm(100))
t2 <- cbind(1:100, rnorm(100))

# Continuous wavelet transform
wt.t1 <- wt(t1)

# Plot power
# Make room to the right for the color bar
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
plot(wt.t1, plot.cb = TRUE, plot.phase = FALSE)

# Cross-wavelet transform
xwt.t1t2 <- xwt(t1, t2)

# Plot cross-wavelet
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
plot(xwt.t1t2, plot.cb = TRUE)

# Example of bias-correction
t1 <- sin(seq(0, 2 * 5 * pi, length.out = 1000))
t2 <- sin(seq(0, 2 * 15 * pi, length.out = 1000))
t3 <- sin(seq(0, 2 * 40 * pi, length.out = 1000))

# This aggregate time series should have the same power
# at three distinct periods
s <- t1 + t2 + t3

# Compare plots to see bias-effect on CWT:
# biased power spectrum artificially
# reduces the power of higher-frequency fluctuations.
wt1 <- wt(cbind(1:1000, s))
par(mfrow = c(1,2))
plot(wt1, type = "power.corr.norm", main = "Bias-corrected")
plot(wt1, type = "power.norm", main = "Biased")

# Compare plots to see bias-effect on XWT:
# biased power spectrum artificially
# reduces the power of higher-frequency fluctuations.
x1 <- xwt(cbind(1:1000, s), cbind(1:1000, s))
par(mfrow = c(1,2))

plot(x1, type = "power.corr.norm", main = "Bias-corrected")
plot(x1, type = "power.norm", main = "Biased")
