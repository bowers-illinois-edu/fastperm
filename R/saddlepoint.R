## saddlepoint.R --- Lugannani-Rice inversion of a fastperm_cgf into a tail
## probability. This is Tier 1: closed-form, no re-randomization.
##
## Given the exact permutation CGF (cgf-stratified.R), solve the saddlepoint
## equation K'(theta_hat) = t and apply the Lugannani-Rice formula. We evaluate
## the small side of the distribution in log space (Mills-ratio form) so deep
## tail probabilities -- far smaller than any feasible Monte-Carlo resolution --
## stay accurate, and obtain the large side by complement.

## solve K'(theta_hat) = t; K is convex so K' is increasing and the root is
## unique for t strictly inside the support. The bracket is capped so a t at the
## support edge returns a large finite theta instead of diverging.
.solve_saddle <- function(cg, t) {
  mu0 <- cg$mean; d1 <- cg$d1
  if (abs(t - mu0) < 1e-9) return(0)
  if (t > mu0) {
    lo <- 0; hi <- 1
    while (d1(hi) < t && hi < 1e6) hi <- hi * 2
    if (d1(hi) < t) return(hi)
  } else {
    hi <- 0; lo <- -1
    while (d1(lo) > t && lo > -1e6) lo <- lo * 2
    if (d1(lo) > t) return(lo)
  }
  stats::uniroot(function(th) d1(th) - t, c(lo, hi), tol = 1e-10)$root
}

## log P(T >= t) for t >= mu (the small upper tail), Lugannani-Rice in Mills-ratio
## form so it never underflows. Near theta_hat = 0 fall back to the normal tail.
.lr_upper_log <- function(cg, t) {
  th <- .solve_saddle(cg, t)
  if (abs(th) < 1e-6)
    return(stats::pnorm((t - cg$mean) / sqrt(cg$d2(0)), lower.tail = FALSE, log.p = TRUE))
  w <- sign(th) * sqrt(2 * (th * t - cg$cgf(th))); u <- th * sqrt(cg$d2(th))
  R <- exp(stats::pnorm(w, lower.tail = FALSE, log.p = TRUE) - stats::dnorm(w, log = TRUE))
  stats::dnorm(w, log = TRUE) + log(R + 1 / u - 1 / w)   # tail = phi(w)[R + 1/u - 1/w]
}

## log P(T <= t) for t <= mu (the small lower tail), the reflected Mills form
.lr_lower_log <- function(cg, t) {
  th <- .solve_saddle(cg, t)
  if (abs(th) < 1e-6)
    return(stats::pnorm((t - cg$mean) / sqrt(cg$d2(0)), lower.tail = TRUE, log.p = TRUE))
  w <- sign(th) * sqrt(2 * (th * t - cg$cgf(th))); u <- th * sqrt(cg$d2(th))
  M <- exp(stats::pnorm(w, log.p = TRUE) - stats::dnorm(w, log = TRUE))  # Phi(w)/phi(w)
  stats::dnorm(w, log = TRUE) + log(M + 1 / w - 1 / u)   # tail = phi(w)[M + 1/w - 1/u]
}

#' Saddlepoint tail probability of a linear statistic
#'
#' Inverts a [fastperm_linear_cgf()] object into a one-sided permutation tail
#' probability `P(T >= t)` (default) or `P(T <= t)`, using the Lugannani-Rice
#' saddlepoint approximation. Probabilities outside the statistic's support are
#' returned exactly (0 or 1). The small tail is computed in log space, so
#' `log.p = TRUE` returns accurate values far past the reach of simulation.
#'
#' @param cgf a `fastperm_cgf` object from [fastperm_linear_cgf()].
#' @param t the value at which to evaluate the tail.
#' @param lower.tail if `FALSE` (default) return `P(T >= t)`; if `TRUE`, `P(T <= t)`.
#' @param log.p if `TRUE`, return the natural log of the probability.
#' @return a probability in `[0, 1]` (or its log if `log.p = TRUE`).
#' @examples
#' cg <- fastperm_linear_cgf(c(1, 2, 3, 4), c(0, 1, 0, 1), rep(1, 4))
#' saddlepoint_tail(cg, 6)              # P(T >= 6)
#' @export
saddlepoint_tail <- function(cgf, t, lower.tail = FALSE, log.p = FALSE) {
  cg <- cgf
  if (!inherits(cg, "fastperm_cgf"))
    stop("`cgf` must be a fastperm_cgf object from fastperm_linear_cgf().", call. = FALSE)
  mu0 <- cg$mean
  if (!lower.tail) {                               # P(T >= t)
    lp <- if (t > cg$supmax) -Inf
          else if (t <= cg$supmin) 0
          else if (t >= mu0) .lr_upper_log(cg, t)
          else log1p(-exp(.lr_lower_log(cg, t)))   # > 0.5, by complement
  } else {                                         # P(T <= t)
    lp <- if (t < cg$supmin) -Inf
          else if (t >= cg$supmax) 0
          else if (t <= mu0) .lr_lower_log(cg, t)
          else log1p(-exp(.lr_upper_log(cg, t)))
  }
  if (log.p) lp else exp(lp)
}
