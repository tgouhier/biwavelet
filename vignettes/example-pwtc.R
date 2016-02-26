library(biwavelet)

y <- cbind(1:100, rnorm(100))
x1 <- cbind(1:100, rnorm(100))
x2 <- cbind(1:100, rnorm(100))

# Partial wavelet coherence of y and x1
pwtc.yx1 <- pwtc(y, x1, x2, nrands = 0)

# Partial wavelet coherence of y and x2
pwtc.yx2 <- pwtc(y, x2, x1, nrands = 0)

# Plot partial wavelet coherence and phase difference (arrows)
# Make room to the right for the color bar
par(mfrow = c(2,1), oma = c(4, 0, 0, 1),
    mar = c(1, 4, 4, 5), mgp = c(1.5, 0.5, 0))

plot(pwtc.yx1, xlab = "", plot.cb = TRUE,
     main = "Partial wavelet coherence of y and x1 | x2")

plot(pwtc.yx2, plot.cb = TRUE,
     main = "Partial wavelet coherence of y and x2 | x1")
