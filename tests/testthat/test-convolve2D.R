context("Fast column-wise convolution of a matrix (convolve2D)")

# Tests ==========

test_that("Original convolve2D should be equal to convolve2D_typeopen", {
  expect_equal(
    convolve2D_typeopen(x,y),
    convolve2D(x,y, type = "open"))
})


