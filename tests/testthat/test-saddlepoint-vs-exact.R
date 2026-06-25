## The whole point of the saddlepoint is to approximate the EXACT permutation
## tail probability without enumerating re-randomizations. Here the strata are
## small enough to enumerate the entire permutation orbit, so the tail is exact
## (zero Monte-Carlo noise) and any gap is pure saddlepoint error. We require
## agreement to a few digits across the body and tail, and that the log-space
## tail equals the direct tail.

## exact within-stratum permutation distribution of T = sum_i z_i c_i, by
## convolving each stratum's vector of all m-subset sums (orbit = product)
exact_T <- function(score, treatment, strata) {
  per <- lapply(split(seq_along(score), strata), function(ix) {
    m <- sum(treatment[ix]); combn(score[ix], m, FUN = sum)
  })
  Reduce(function(a, v) as.vector(outer(a, v, "+")), per)
}

test_that("saddlepoint tail agrees with the exact enumerated tail", {
  set.seed(3)
  ## 5 strata of 5, 2 treated each: orbit = choose(5,2)^5 = 10^5 = 100000
  strata <- factor(rep(1:5, each = 5))
  score <- rnorm(25)
  treatment <- numeric(25)
  for (ix in split(seq_along(strata), strata)) treatment[ix[1:2]] <- 1

  cg <- fastperm_linear_cgf(score, treatment, strata)
  Tall <- exact_T(score, treatment, strata)
  exact_upper <- function(t) mean(Tall >= t)

  ## compare across body-to-tail quantiles of the exact distribution
  thr <- as.numeric(quantile(Tall, probs = c(.5, .75, .9, .95, .99), type = 1))
  spa <- vapply(thr, function(t) saddlepoint_tail(cg, t), numeric(1))
  exa <- vapply(thr, exact_upper, numeric(1))
  expect_lt(max(abs(spa - exa)), 0.01)            # a few digits even at n=25

  ## log-space tail must equal the direct tail wherever the direct one is usable
  lg <- vapply(thr, function(t) saddlepoint_tail(cg, t, log.p = TRUE), numeric(1))
  expect_equal(exp(lg), spa, tolerance = 1e-8)
})

test_that("tails outside the support are 0 (upper) and 1 (lower)", {
  score <- c(-2, -1, 1, 2, 3)
  treatment <- c(0, 0, 0, 1, 1)
  strata <- factor(rep(1, 5))
  cg <- fastperm_linear_cgf(score, treatment, strata)
  expect_equal(saddlepoint_tail(cg, cg$supmax + 1), 0)               # impossible upper
  expect_equal(saddlepoint_tail(cg, cg$supmin - 1, lower.tail = TRUE), 0)
  expect_equal(saddlepoint_tail(cg, cg$supmin), 1)                   # everything >= min
})
