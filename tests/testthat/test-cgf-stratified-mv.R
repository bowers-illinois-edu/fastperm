## The joint CGF is the enabler for the multivariate (quadratic-form, max)
## saddlepoints. Its contract: the gradient and Hessian at s = 0 must be the
## exact permutation mean vector and covariance of the representation vector --
## the multivariate analogue of the 1-D cumulant check. This pins that the
## projected elementary-symmetric DP really computes the joint cumulants.

mk_treat <- function(strata, m) {
  z <- numeric(length(strata))
  for (ix in split(seq_along(strata), strata)) z[ix[seq_len(m)]] <- 1
  z
}

test_that("joint CGF gradient and Hessian at 0 are the mean vector and covariance", {
  set.seed(7)
  strata <- factor(rep(1:8, each = 6))
  n <- length(strata); d <- 4L
  scores <- matrix(rnorm(n * d), n, d)
  ## centre within stratum so the statistics have permutation mean ~ 0
  for (k in seq_len(d)) scores[, k] <- scores[, k] - ave(scores[, k], strata)
  treatment <- mk_treat(strata, 3)

  cg <- fastperm_linear_cgf_mv(scores, treatment, strata)
  K <- cg$cgf
  h <- 1e-4

  ## numerical gradient at 0 -> mean vector
  grad <- vapply(seq_len(d), function(k) {
    e <- numeric(d); e[k] <- h
    (K(e) - K(-e)) / (2 * h)
  }, numeric(1))
  expect_equal(unname(grad), unname(cg$mean), tolerance = 1e-5)

  ## numerical Hessian at 0 -> covariance matrix
  hess <- matrix(0, d, d)
  for (k in seq_len(d)) for (l in seq_len(d)) {
    ek <- numeric(d); ek[k] <- h; el <- numeric(d); el[l] <- h
    hess[k, l] <- (K(ek + el) - K(ek - el) - K(-ek + el) + K(-ek - el)) / (4 * h^2)
  }
  expect_equal(hess, unname(cg$cov), tolerance = 1e-4)
})

test_that("the projected CGF agrees with the univariate CGF on the projection", {
  ## K_T(s) must equal the 1-D CGF of the projected scores at theta = 1, since
  ## projecting onto s and tilting at theta = 1 is the same operation
  set.seed(8)
  strata <- factor(rep(1:5, each = 5))
  n <- length(strata); d <- 3L
  scores <- matrix(rnorm(n * d), n, d)
  treatment <- mk_treat(strata, 2)
  s <- c(0.7, -0.4, 0.2)

  cg_mv <- fastperm_linear_cgf_mv(scores, treatment, strata)
  cg_1d <- fastperm_linear_cgf(as.numeric(scores %*% s), treatment, strata)
  expect_equal(cg_mv$cgf(s), cg_1d$cgf(1), tolerance = 1e-10)
})
