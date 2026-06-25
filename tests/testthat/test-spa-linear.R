## fastperm_spa_linear() is the user-facing convenience: it builds the CGF for a
## linear statistic, evaluates the observed value, and returns a saddlepoint
## p-value. These tests check the contract a permutation p-value must satisfy and
## that the one- and two-sided alternatives correspond to the right tails.

mk_treat <- function(strata, m) {
  z <- numeric(length(strata))
  for (ix in split(seq_along(strata), strata)) z[ix[seq_len(m)]] <- 1
  z
}

test_that("the observed statistic and moments are reported correctly", {
  score <- c(1, 2, 3, 4)          # one stratum, 2 treated
  treatment <- c(0, 1, 0, 1)      # treated = units 2 and 4
  strata <- factor(rep(1, 4))
  res <- fastperm_spa_linear(score, treatment, strata, alternative = "greater")
  expect_equal(res$statistic, 2 + 4)            # sum of treated scores
  expect_equal(res$mean, (2 / 4) * sum(score))  # (m/n) * sum
})

test_that("p-values are valid probabilities and alternatives point the right way", {
  set.seed(4)
  strata <- factor(rep(1:4, each = 6))
  raw <- rnorm(24)
  score <- raw - ave(raw, strata)               # centred -> mean 0
  ## a strongly (but not maximally) upper-tail observed assignment: treat the
  ## 1st, 2nd, and 4th highest-score units in each stratum (skip the 3rd so the
  ## observed statistic stays strictly inside the support).
  treatment <- numeric(24)
  for (ix in split(seq_along(strata), strata))
    treatment[ix[order(score[ix], decreasing = TRUE)[c(1, 2, 4)]]] <- 1

  g <- fastperm_spa_linear(score, treatment, strata, alternative = "greater")$p.value
  l <- fastperm_spa_linear(score, treatment, strata, alternative = "less")$p.value
  t2 <- fastperm_spa_linear(score, treatment, strata, alternative = "two.sided")$p.value

  expect_true(g > 0 && g <= 1 && l > 0 && l <= 1 && t2 > 0 && t2 <= 1)
  expect_lt(g, l)                       # upper-tail evidence: greater p < less p
  expect_equal(g + l, 1, tolerance = 0.02)        # complementary one-sided tails
  expect_equal(t2, 2 * min(g, l), tolerance = 0.02)  # two-sided ~ 2x smaller tail
})

test_that("greater/less match saddlepoint_tail on the built CGF", {
  set.seed(5)
  strata <- factor(rep(1:5, each = 5))
  score <- rnorm(25)
  treatment <- mk_treat(strata, 2)
  cg <- fastperm_linear_cgf(score, treatment, strata)
  t_obs <- sum(treatment * score)
  res_g <- fastperm_spa_linear(score, treatment, strata, alternative = "greater")
  res_l <- fastperm_spa_linear(score, treatment, strata, alternative = "less")
  expect_equal(res_g$p.value, saddlepoint_tail(cg, t_obs, lower.tail = FALSE))
  expect_equal(res_l$p.value, saddlepoint_tail(cg, t_obs, lower.tail = TRUE))
})
