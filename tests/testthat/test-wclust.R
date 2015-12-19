context("Wavelet clustering (wclust)")

# Setup variables =====

# Sample wavelets
t1 <- cbind(1:100, sin(seq(0, 10 * 2 * pi, length.out = 100)))
t2 <- cbind(1:100, sin(seq(0, 10 * 2 * pi, length.out = 100) + 0.1 * pi))
t3 <- cbind(1:100, rnorm(100)) # white noise

# Compute wavelet spectra
wt.t1 <- wt(t1)
wt.t2 <- wt(t2)
wt.t3 <- wt(t3)

# Store all wavelet spectra into array
w.arr <- array(dim = c(3, NROW(wt.t1$wave), NCOL(wt.t1$wave)))
w.arr[1, , ] <- wt.t1$wave
w.arr[2, , ] <- wt.t2$wave
w.arr[3, , ] <- wt.t3$wave

# Tests ==========

test_that("Basic test of wclust without progress bar", {
  # Compute dissimilarity and distance matrices
  c <- wclust(w.arr, quiet = TRUE)

  expect_true(is.matrix(c$diss.mat))
  expect_equal(dim(c$diss.mat), c(3,3))
  expect_equal(class(c$dist.mat), "dist")
})

test_that("Progressbar should not cause errors", {
  expect_output(
    # Compute dissimilarity and distance matrices
    out <- wclust(w.arr),
    regexp = "\\|=+=\\| 100%"
  )
})
