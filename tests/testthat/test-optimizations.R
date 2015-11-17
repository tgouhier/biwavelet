context("Performance optimizations")

test_that("Optimized version of wt.bases.dog is equal to original", {
  for (param in 0:4) {
    expect_equal(
      biwavelet:::wt.bases.dog(1:10, 2, param),
      biwavelet:::rcpp_wt_bases_dog(1:10, 2, param)
    )
  }
})

test_that("replacing seq(1,N,1) with 1:N", {

  m <- J1 <- npad <- n.obs <- num_waves <- 17

  # in wt.bases
  expect_equal(
    2:(2 * m - 1),
    seq(2, 2 * m - 1, 1) )

  # in wt.sig and wt
  expect_equal(
    1:(J1 + 1),
    seq(1, J1 + 1, 1) )

  # in smooth.wavelet
  expect_equal(
    1:as.integer(npad / 2),
    seq(1, as.integer(npad / 2), 1) )

  # in smooth.wavelet
  expect_equal(
    as.integer((npad - 1) / 2):1,
    seq(as.integer((npad - 1) / 2), 1, -1) )

  # in wt
  expect_equal(
    1:((n.obs + 1) / 2 - 1),
    seq(1, (n.obs + 1) / 2 - 1) )

  # in wt
  expect_equal(
    floor(n.obs / 2 - 1):1,
    seq(floor(n.obs / 2 - 1), 1, -1) )

  # in wclust
  expect_equal(
    1:num_waves,
    seq(from = 1, to = num_waves, by = 1) )

})
