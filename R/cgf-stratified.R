## cgf-stratified.R --- exact cumulant generating function (CGF) of a linear
## statistic under the within-stratum permutation null. This is the Tier-1
## saddlepoint input.
##
## For a linear statistic T = sum_i z_i c_i referred to the permutation null that
## fixes each stratum's treated count, the contribution of stratum b is the sum
## of the scores over a uniformly chosen m_b-subset of the stratum. Its moment
## generating function is the elementary symmetric polynomial of the tilted
## weights exp(theta c_i), normalised by choose(n_b, m_b); strata are independent,
## so the log-MGFs add. We build the exact CGF K(theta) = sum_b log M_b(theta),
## together with the design's exact mean, variance, and support. The saddlepoint
## inversion in saddlepoint.R turns this into a tail p-value with no
## re-randomization. The mean and variance returned here are the exact
## finite-population (Strasser-Weber) permutation moments, and K'(0), K''(0)
## reproduce them; the saddlepoint's accuracy beyond the normal approximation
## comes from the higher cumulants the full CGF carries.

## stable log(exp(a) + exp(b)), vectorised
.logaddexp <- function(a, b) {
  hi <- pmax(a, b); lo <- pmin(a, b)
  res <- hi + log1p(exp(lo - hi))
  inf <- is.infinite(hi)                     # both -Inf -> -Inf; +Inf -> +Inf
  res[inf] <- hi[inf]
  res
}

## log of the degree-m elementary symmetric polynomial of the weights exp(x),
## via the O(n*m) dynamic-programming recursion in log space (the weights span
## many orders of magnitude at the saddlepoint, so we never form them directly).
## log e_k lives at index k+1; processing one unit updates indices high-to-low.
.log_esym <- function(x, m) {
  logp <- c(0, rep(-Inf, m))                 # log e_0 = 0
  for (xi in x) logp[2:(m + 1)] <- .logaddexp(logp[2:(m + 1)], xi + logp[1:m])
  logp[m + 1]
}

#' Exact CGF of a linear statistic under the within-stratum permutation null
#'
#' Builds the cumulant generating function of the linear statistic
#' `T = sum_i treatment_i * score_i` referred to the permutation null that holds
#' each stratum's treated count fixed (a simple random sample of treated units
#' within each stratum). The returned object carries the exact permutation mean,
#' variance, and support of `T`, plus closures for `K(theta)` and its first two
#' derivatives, which [saddlepoint_tail()] inverts into a tail probability.
#'
#' @param score numeric vector of per-unit scores `c_i` (for a centred
#'   representation, the within-stratum-centred values).
#' @param treatment 0/1 vector; its per-stratum sums give the treated counts that
#'   define the permutation null.
#' @param strata factor or vector of stratum labels, the same length as `score`.
#' @return an object of class `fastperm_cgf`: a list with `cgf`, `d1`, `d2`
#'   (closures in `theta`), and the scalars `mean`, `variance`, `supmin`,
#'   `supmax`.
#' @examples
#' score <- c(1, 2, 3, 4)
#' cg <- fastperm_linear_cgf(score, c(0, 1, 0, 1), rep(1, 4))
#' cg$mean
#' @export
fastperm_linear_cgf <- function(score, treatment, strata) {
  score <- as.numeric(score); treatment <- as.numeric(treatment)
  strata <- as.factor(strata)
  n <- length(score)
  if (length(treatment) != n || length(strata) != n)
    stop("`score`, `treatment`, and `strata` must describe the same units.",
         call. = FALSE)

  idx <- split(seq_len(n), strata)
  ## unname so the accumulated scalars (mean, variance, support, CGF) do not
  ## inherit stratum labels from the split
  mb  <- unname(vapply(idx, function(ix) sum(treatment[ix]), numeric(1)))
  nb  <- unname(lengths(idx))

  ## exact moments and support, accumulated over the independent strata
  mean0 <- 0; var0 <- 0; smax <- 0; smin <- 0
  for (b in seq_along(idx)) {
    ix <- idx[[b]]; c <- score[ix]; m <- mb[b]; nbb <- nb[b]
    mean0 <- mean0 + (m / nbb) * sum(c)
    if (nbb > 1 && m > 0 && m < nbb) {
      w <- m * (nbb - m) / (nbb * (nbb - 1))      # finite-population factor
      var0 <- var0 + w * sum((c - mean(c))^2)
    }
    if (m > 0L) {                                 # support: top/bottom m scores
      sc <- sort(c)                               # ascending
      smax <- smax + sum(sc[(nbb - m + 1L):nbb])
      smin <- smin + sum(sc[1L:m])
    }
  }

  ## K(theta) = sum_b log M_b(theta). A stratum with m in {0, n_b} contributes a
  ## deterministic amount; the rest use the elementary-symmetric MGF.
  cgf <- function(theta) {
    acc <- 0
    for (b in seq_along(idx)) {
      ix <- idx[[b]]; m <- mb[b]; nbb <- nb[b]
      if (m == 0L) next
      if (m == nbb) { acc <- acc + theta * sum(score[ix]); next }
      acc <- acc + .log_esym(theta * score[ix], m) - lchoose(nbb, m)
    }
    acc
  }
  ## central differences; K is smooth and O(1), accurate well past the tail
  ## digits the saddlepoint needs (checked against the closed-form moments)
  h <- 1e-3
  d1 <- function(theta) (cgf(theta + h) - cgf(theta - h)) / (2 * h)
  d2 <- function(theta) (cgf(theta + h) - 2 * cgf(theta) + cgf(theta - h)) / h^2

  structure(list(cgf = cgf, d1 = d1, d2 = d2, mean = mean0, variance = var0,
                 supmin = smin, supmax = smax), class = "fastperm_cgf")
}

#' @export
print.fastperm_cgf <- function(x, ...) {
  cat("<fastperm_cgf> linear statistic under within-stratum permutation\n")
  cat(sprintf("  mean %.4g, sd %.4g, support [%.4g, %.4g]\n",
              x$mean, sqrt(x$variance), x$supmin, x$supmax))
  invisible(x)
}
