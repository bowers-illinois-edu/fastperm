## test-saddlepoint-quadratic.R --- M1 milestone of the Route B quadratic-form
## saddlepoint: the Gaussian-d weighted-chi-square tail of
## Q = (T - mu)' M^{-1} (T - mu), inverted by Lugannani-Rice.
##
## These tests encode the statistical content M1 must satisfy, not just that it
## runs:
##   * the rotation reduction Q = ||d_y||^2 is an exact identity (the design
##     decision that lets the metric be preprocessing, not an inversion argument);
##   * M1 reproduces CompQuadForm::imhof, the exact tail of sum_k lambda_k chi^2_1
##     under Gaussian d (M1 is a saddlepoint *to that* reference);
##   * with metric = V_d the rotated covariance is the identity, so Q collapses to
##     chi^2_rank -- the Hansen-Bowers omnibus balance statistic;
##   * the closed-form enabler moments equal the exact orbit moments.
##
## M1 is validated against the Gaussian-d reference (imhof / chi-square), NOT
## against the exact permutation law. Closing the gap to enumeration is M2, whose
## tests will use the same enumerate_Q oracle.
##
## API under test (to be implemented after this checkpoint):
##   fastperm_spa_quadratic(scores, treatment, strata,
##                          metric = "cov",                 # V_d, or a p x p matrix
##                          method = c("gaussian", "saddlepoint"))
##   -> list(p.value, statistic = Q_obs, rank, eigenvalues, mean = E[Q],
##           metric, method)

## --- fixtures ---------------------------------------------------------------
## one stratum, n = 10, m = 5: orbit choose(10, 5) = 252, large enough that the
## Gaussian-d approximation is reasonable and imhof / chi-square are meaningful
## references; treated = the five high units, so Q_obs sits in the upper tail.
sc_1s <- cbind(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
               c(2, 1, 4, 3, 6, 5, 8, 7, 9, 10))
tr_1s <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
st_1s <- rep(1, 10)

## two strata, n = 8, two treated per stratum
sc_2s <- cbind(c(1, 2, 3, 4, 2, 4, 6, 8),
               c(1, 1, 3, 5, 1, 2, 2, 7))
tr_2s <- c(0, 0, 1, 1, 0, 0, 1, 1)
st_2s <- c(1, 1, 1, 1, 2, 2, 2, 2)


## --- oracle-level checks (independent of the M1 implementation) --------------

test_that("the rotation reduction Q = ||d_y||^2 is an exact identity", {
  ## Load-bearing simplification: because R R' = M^{-1}, the quadratic form is a
  ## squared norm in rotated coordinates, so the saddlepoint core can be built on
  ## the sum of squares with the metric as preprocessing. If this drifts from
  ## machine precision the whole design is wrong.
  for (M in list(diag(2), stats::cov(sc_1s))) {
    en <- enumerate_Q(sc_1s, tr_1s, st_1s, M)
    expect_equal(en$Q, en$Q_rot, tolerance = 1e-9)
  }
})

test_that("the enabler's permutation covariance matches the orbit covariance", {
  ## fastperm_linear_cgf_mv returns the closed-form Strasser-Weber moments; they
  ## must equal the population moments over the fully enumerated orbit.
  en <- enumerate_Q(sc_2s, tr_2s, st_2s, diag(2))
  cg <- fastperm_linear_cgf_mv(sc_2s, tr_2s, st_2s)
  expect_equal(unname(cg$cov),  unname(en$Sigma), tolerance = 1e-9)
  expect_equal(unname(cg$mean), unname(en$mu),    tolerance = 1e-9)
})

test_that("E[Q] equals trace(M^{-1} Sigma), the sum of the chi-square weights", {
  ## Exact-moment check on the oracle: the mean of the enumerated Q is the trace
  ## of the rotated covariance, i.e. the sum of the weighted-chi-square weights.
  M  <- stats::cov(sc_1s)
  en <- enumerate_Q(sc_1s, tr_1s, st_1s, M)
  A  <- chol2inv(chol(M))
  expect_equal(mean(en$Q), sum(diag(A %*% en$Sigma)), tolerance = 1e-9)
})


## --- Lugannani-Rice near the mean: the Daniels removable-singularity limit ----

test_that(".quad_lr_upper returns the Daniels skewness limit at the mean, not 1/2", {
  ## At q = E[Q] the Lugannani-Rice 1/u - 1/w term is a removable 0/0 whose finite
  ## limit (Daniels 1987) keeps the third cumulant:
  ##   P(Q >= E[Q]) = 1/2 - K'''(0) / (6 sqrt(2 pi) K''(0)^{3/2}).
  ## For a right-skewed weighted chi-square (K'''(0) > 0) this is strictly below
  ## the 1/2 a normal approximation gives. The old fallback returned the plain
  ## normal tail and dropped this skewness term; this is the regression guard.
  lambda <- c(3, 1, 0.5)                        # distinct weights -> right-skewed Q
  m0 <- sum(lambda)
  k2 <- 2 * sum(lambda^2)                       # K''(0)
  k3 <- 8 * sum(lambda^3)                       # K'''(0)
  daniels <- 0.5 - k3 / (6 * sqrt(2 * pi) * k2^(3 / 2))
  expect_equal(.quad_lr_upper(lambda, m0), daniels, tolerance = 1e-9)
  expect_lt(.quad_lr_upper(lambda, m0), 0.5)
})

test_that(".quad_lr_upper is continuous and monotone across the mean", {
  ## The fix must not introduce a jump at the |theta|-small band edge: the
  ## Lugannani-Rice value approaches the Daniels limit as q -> E[Q] from both
  ## sides. (The old code jumped from ~1/2 in the band to the LR curve outside.)
  lambda <- c(3, 1, 0.5)
  m0  <- sum(lambda); sd0 <- sqrt(2 * sum(lambda^2))
  qs  <- m0 + seq(-0.3, 0.3, length.out = 41) * sd0
  ps  <- vapply(qs, function(q) .quad_lr_upper(lambda, q), numeric(1))
  expect_true(all(diff(ps) < 0))                # strictly decreasing in q
  expect_lt(max(abs(diff(ps))), 0.05)           # no jump at the band edge
})


## --- M1 specification (fails until fastperm_spa_quadratic is implemented) -----

test_that("M1 with metric = 'cov' reduces to the chi-square tail", {
  ## M = V_d gives rotated covariance I, so under Gaussian d, Q ~ chi^2_rank.
  ## This is the Hansen-Bowers omnibus balance statistic; M1 must match pchisq.
  res <- fastperm_spa_quadratic(sc_1s, tr_1s, st_1s,
                                metric = "cov", method = "gaussian")
  expect_equal(res$rank, 2L)
  ## Lugannani-Rice approximates the exact chi-square tail; this fixture puts
  ## Q_obs near the 3% tail, where LR carries about 1% relative error. 2e-2 is
  ## the calibrated relative tolerance -- loose enough for the method's own
  ## accuracy, tight enough to catch a wrong formula (which would be off by tens
  ## of percent or more).
  expect_equal(res$p.value,
               stats::pchisq(res$statistic, df = res$rank, lower.tail = FALSE),
               tolerance = 2e-2)
})

test_that("M1 reproduces imhof for a non-identity metric", {
  skip_if_not_installed("CompQuadForm")
  ## M = Sigma_x (here the sample covariance of the scores): Q is a genuine
  ## weighted sum of chi-squares; M1's saddlepoint tail must match imhof's exact
  ## inversion of that same Gaussian-d distribution.
  res <- fastperm_spa_quadratic(sc_1s, tr_1s, st_1s,
                                metric = stats::cov(sc_1s), method = "gaussian")
  imhof <- CompQuadForm::imhof(res$statistic, lambda = res$eigenvalues)$Qq
  ## imhof is exact for the weighted chi-square; M1 is the Lugannani-Rice
  ## saddlepoint to that same distribution, so they agree up to LR's ~1% tail
  ## relative error.
  expect_equal(res$p.value, imhof, tolerance = 2e-2)
})

test_that("M1 reports a statistic and weights consistent with the oracle", {
  M   <- stats::cov(sc_2s)
  en  <- enumerate_Q(sc_2s, tr_2s, st_2s, M)
  res <- fastperm_spa_quadratic(sc_2s, tr_2s, st_2s,
                                metric = M, method = "gaussian")
  expect_equal(res$statistic, en$q_obs, tolerance = 1e-9)
  expect_true(res$p.value > 0 && res$p.value < 1)
  ## eigenvalues of the rotated covariance sum to E[Q]
  expect_equal(sum(res$eigenvalues), mean(en$Q), tolerance = 1e-9)
})

test_that("the saddlepoint (M2) method is not yet available", {
  ## Documents the staged plan: the non-normal permutation saddlepoint is M2.
  ## Until it lands, requesting it must fail loudly rather than silently fall
  ## back to the Gaussian-d answer.
  expect_error(
    fastperm_spa_quadratic(sc_1s, tr_1s, st_1s, method = "saddlepoint"),
    "not yet implemented")
})
