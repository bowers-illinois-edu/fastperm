## The saddlepoint is built on the EXACT permutation cumulant generating
## function, so its first two cumulants must equal the exact finite-population
## (Strasser-Weber) permutation mean and variance of the linear statistic. These
## tests pin that down: the elementary-symmetric-DP CGF must reproduce the
## closed-form moments and the exact support of the statistic.

## mark m treated units in each stratum (length-n 0/1 vector, original order)
mk_treat <- function(strata, m) {
  z <- numeric(length(strata))
  for (ix in split(seq_along(strata), strata)) z[ix[seq_len(m)]] <- 1
  z
}

## independent closed-form within-stratum permutation moments of T = sum_i z_i c_i
## (treated count fixed per stratum) -- the ground truth the CGF must reproduce
sw_moments <- function(score, treatment, strata) {
  mu <- 0; v <- 0
  for (ix in split(seq_along(score), strata)) {
    c <- score[ix]; n <- length(c); m <- sum(treatment[ix])
    mu <- mu + (m / n) * sum(c)
    if (n > 1 && m > 0 && m < n) {
      w <- m * (n - m) / (n * (n - 1))
      v <- v + w * sum((c - mean(c))^2)
    }
  }
  c(mean = mu, var = v)
}

test_that("CGF cumulants match the closed-form Strasser-Weber moments", {
  set.seed(1)
  strata <- factor(rep(1:6, each = 5))
  score <- rnorm(30)
  treatment <- mk_treat(strata, 2)

  cg <- fastperm_linear_cgf(score, treatment, strata)
  truth <- sw_moments(score, treatment, strata)

  ## the design's reported mean/variance are the exact moments ...
  expect_equal(cg$mean, unname(truth["mean"]), tolerance = 1e-10)
  expect_equal(cg$variance, unname(truth["var"]), tolerance = 1e-10)
  ## ... and the CGF derivatives at theta = 0 reproduce them (the DP is correct)
  expect_equal(cg$d1(0), unname(truth["mean"]), tolerance = 1e-6)
  expect_equal(cg$d2(0), unname(truth["var"]), tolerance = 1e-5)
})

test_that("centred scores give permutation mean zero", {
  set.seed(2)
  strata <- factor(rep(1:4, each = 6))
  raw <- rnorm(24)
  ## centre within stratum, as riposte's representations are
  score <- raw - ave(raw, strata)
  treatment <- mk_treat(strata, 3)
  cg <- fastperm_linear_cgf(score, treatment, strata)
  expect_equal(cg$mean, 0, tolerance = 1e-10)
})

test_that("support equals the extreme attainable statistic values", {
  ## one stratum, 4 units, 2 treated: max = sum of the two largest scores
  score <- c(-3, -1, 2, 5)
  treatment <- c(0, 0, 1, 1)
  strata <- factor(rep(1, 4))
  cg <- fastperm_linear_cgf(score, treatment, strata)
  expect_equal(cg$supmax, 2 + 5)     # two largest
  expect_equal(cg$supmin, -3 + -1)   # two smallest
})
