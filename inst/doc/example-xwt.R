library(biwavelet)

t1 <- cbind(1:100, rnorm(100))
t2 <- cbind(1:100, rnorm(100))

# Compute Cross-wavelet
xwt.t1t2 <- xwt(t1, t2)
plot(xwt.t1t2, plot.cb = TRUE, plot.phase = TRUE,
     main = "Plot cross-wavelet and phase difference (arrows)")

# Real data
data(enviro.data)

# Cross-wavelet of MEI and NPGO
xwt.mei.npgo <- xwt(subset(enviro.data, select = c("date", "mei")),
                    subset(enviro.data, select = c("date", "npgo")))
                    
# Make room to the right for the color bar
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
plot(xwt.mei.npgo, plot.cb = TRUE, plot.phase = TRUE)
