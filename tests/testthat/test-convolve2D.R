context("Fast column-wise convolution of a matrix (convolve2D)")

# Tests ==========

test_that("Original convolve2D should be equal to convolve2D_typeopen", {
  x <- matrix(1, nrow = 10, ncol = 3)
  y <- c(1,1,2)
  expect_equal(
    convolve2D_typeopen(x,y),
    convolve2D(x,y, type = "open"))
})

test_that("Circular convolve2D should identify length mismatch", {
  x <- matrix(1, nrow = 10, ncol = 3)
  y <- c(1,1,2)
  expect_error(
    convolve2D(x,y, type = "circular"),
    regexp = "length mismatch")
})
