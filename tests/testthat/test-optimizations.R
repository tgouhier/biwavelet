context("Performance optimizations")

test_that("Optimized version of wt.bases.dog is equal to original", {
  for (param in 1:10) {

    # integer scale
    expect_equal(
      biwavelet:::wt.bases.dog(1:10, 2, param),
      biwavelet:::rcpp_wt_bases_dog(1:10, 2, param)
    )

    # floating point scale
    expect_equal(
      biwavelet:::wt.bases.dog(1:10, 1 / 3, param),
      biwavelet:::rcpp_wt_bases_dog(1:10, 1 / 3, param)
    )
  }
})

test_that("Optimized version of wt.bases.paul is equal to original", {
  for (param in 1:10) {

    # integer scale
    expect_equal(
      biwavelet:::wt.bases.paul(1:10, 2, param),
      biwavelet:::rcpp_wt_bases_paul(1:10, 2, param)
    )

    # floating point scale
    expect_equal(
      biwavelet:::wt.bases.paul(1:10, 1 / 3, param),
      biwavelet:::rcpp_wt_bases_paul(1:10, 1 / 3, param)
    )
  }
})

test_that("Optimized version of wt.bases.morlet is equal to original", {
  for (param in 1:10) {

    # integer scale
    expect_equal(
      biwavelet:::wt.bases.morlet(1:10, 2, param),
      biwavelet:::rcpp_wt_bases_morlet(1:10, 2, param)
    )

    # floating point scale
    expect_equal(
      biwavelet:::wt.bases.morlet(1:10, 1 / 3, param),
      biwavelet:::rcpp_wt_bases_morlet(1:10, 1 / 3, param)
    )
  }
})

test_that("Parameter m outside supported interval should fail", {

  FMIN <- -2 # lower bound when the function should fail
  FMAX <- 11 # upper bound when the function should fail
  ERRMSG <- "must be within"

  # dog
  expect_error(biwavelet:::rcpp_wt_bases_dog(1:10, 2, FMIN), regexp = ERRMSG)
  expect_error(biwavelet:::rcpp_wt_bases_dog(1:10, 2, FMAX), regexp = ERRMSG)

  # morlet
  expect_error(biwavelet:::rcpp_wt_bases_morlet(1:10, 2, FMIN), regexp = ERRMSG)
  expect_error(biwavelet:::rcpp_wt_bases_morlet(1:10, 2, FMAX), regexp = ERRMSG)

  # paul
  expect_error(biwavelet:::rcpp_wt_bases_paul(1:10, 2, FMIN), regexp = ERRMSG)
  expect_error(biwavelet:::rcpp_wt_bases_paul(1:10, 2, FMAX), regexp = ERRMSG)
})

test_that("Default 'param' values for rcpp_wt_bases", {
  expect_equal(
    biwavelet:::rcpp_wt_bases_dog(1:10, 2, -1),
    biwavelet:::rcpp_wt_bases_dog(1:10, 2, 2))
  expect_equal(
    biwavelet:::rcpp_wt_bases_morlet(1:10, 2, -1),
    biwavelet:::rcpp_wt_bases_morlet(1:10, 2, 6))
  expect_equal(
    biwavelet:::rcpp_wt_bases_paul(1:10, 2, -1),
    biwavelet:::rcpp_wt_bases_paul(1:10, 2, 4))
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
