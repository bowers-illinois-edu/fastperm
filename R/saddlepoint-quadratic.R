## saddlepoint-quadratic.R --- saddlepoint tail probability of a quadratic form
## in a vector of linear statistics under the within-stratum permutation null.
## This is the Route B entry point. The statistic is
##   Q = (T - mu)' M^{-1} (T - mu),   T_c = sum_i z_i score_{i,c},
## a Mahalanobis-type distance of the treated-group score vector from its
## permutation mean in a user-chosen metric M. With M the permutation covariance
## this is the omnibus covariate-balance statistic of Hansen and Bowers (2008).
##
## The metric enters only as a rotation. For any R with R R' = M^{-1},
##   Q = || R'(T - mu) ||^2,
## the squared norm of the rotated, centred linear-statistic vector. So the core
## works on a sum of squares and the metric is preprocessing (an exact identity,
## checked to machine precision against full enumeration in the tests).
##
## Two inversions, staged:
##  * method = "gaussian" (M1, implemented here). Treat the centred vector as
##    Gaussian with the exact permutation covariance; then Q is a weighted sum of
##    independent chi-square_1 variables, sum_k lambda_k chi^2_1, with lambda the
##    eigenvalues of the rotated covariance. Its CGF K(theta) = -1/2 sum log(1 -
##    2 theta lambda_k) is closed form, inverted by Lugannani-Rice. This
##    reproduces CompQuadForm::imhof, and the chi-square tail when M is the
##    permutation covariance, but inherits the Gaussian approximation -- which is
##    inaccurate for the small, skewed permutation distributions that balance
##    testing lives in.
##  * method = "saddlepoint" (M2, forthcoming). The non-normal permutation
##    saddlepoint from the full joint CGF (cgf-stratified-mv.R), which closes the
##    gap to the exact permutation distribution.

## M^{-1/2} as a p x r matrix R with R R' = M^{-1}, via the spectral
## decomposition; numerically-zero eigenvalues are dropped so a rank-deficient
## metric projects onto its range (r = rank of M).
.quad_inv_sqrt <- function(M) {
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  keep <- e$values > sqrt(.Machine$double.eps) * max(e$values, 1)
  if (!any(keep)) stop("the metric has no positive eigenvalues.", call. = FALSE)
  e$vectors[, keep, drop = FALSE] %*%
    diag(1 / sqrt(e$values[keep]), sum(keep), sum(keep))
}

## P(Q >= q) for Q = sum_k lambda_k chi^2_1 (lambda > 0) by Lugannani-Rice. The
## CGF K(theta) = -1/2 sum log(1 - 2 theta lambda_k) has a pole at
## theta = 1/(2 max lambda), so the saddlepoint root is bracketed below the pole
## (q above the mean) or on the negative axis (q below the mean).
##
## The correction term 1/u - 1/w is a removable 0/0 at q = E[Q] and loses
## precision to cancellation (w -> u) over a neighbourhood of the mean, not just
## at the point. There we replace it by its Taylor series in theta, c0 + c1 theta
## (Daniels 1987), with
##   c0 = -K'''(0) / (6 K''(0)^{3/2}),
##   c1 = 5 K'''(0)^2 / (24 K''(0)^{5/2}) - K''''(0) / (8 K''(0)^{3/2}).
## The leading term gives P(Q >= E[Q]) = 1/2 + phi(0) c0 = 1/2 - K'''(0) /
## (6 sqrt(2 pi) K''(0)^{3/2}) -- the skewness-corrected limit a plain normal tail
## (1/2 at the mean) would miss. The series is used while |w| < 0.1, where the raw
## difference is unreliable and the linear term is accurate to ~1e-3; the direct
## formula is used in the tails, where it is exact and the series would drift.
.quad_lr_upper <- function(lambda, q) {
  K   <- function(th) -0.5 * sum(log(1 - 2 * th * lambda))
  Kp  <- function(th) sum(lambda / (1 - 2 * th * lambda))
  Kpp <- function(th) sum(2 * lambda^2 / (1 - 2 * th * lambda)^2)
  m0  <- sum(lambda)
  k2  <- 2  * sum(lambda^2)                               # K''(0)
  k3  <- 8  * sum(lambda^3)                               # K'''(0)
  k4  <- 48 * sum(lambda^4)                               # K''''(0)
  c0  <- -k3 / (6 * k2^(3 / 2))
  c1  <- 5 * k3^2 / (24 * k2^(5 / 2)) - k4 / (8 * k2^(3 / 2))
  if (q <= 0) return(1)                                   # Q is non-negative
  if (abs(q - m0) < 1e-9 * (1 + abs(m0)))
    return(0.5 + stats::dnorm(0) * c0)                    # q = E[Q]: w = 0
  if (q > m0) {
    pole <- 1 / (2 * max(lambda))
    th <- stats::uniroot(function(t) Kp(t) - q,
                         c(1e-12, pole * (1 - 1e-12)))$root
  } else {
    lo <- -1                                              # expand down to bracket
    while (Kp(lo) > q && lo > -1e10) lo <- lo * 2
    th <- stats::uniroot(function(t) Kp(t) - q, c(lo, -1e-12))$root
  }
  w    <- sign(th) * sqrt(2 * (th * q - K(th)))
  term <- if (abs(w) < 0.1) c0 + c1 * th                  # near-mean series
          else 1 / (th * sqrt(Kpp(th))) - 1 / w           # direct, in the tails
  p <- stats::pnorm(w, lower.tail = FALSE) + stats::dnorm(w) * term
  min(1, max(0, p))
}

#' Saddlepoint tail probability of a quadratic form under the permutation null
#'
#' Computes the within-stratum permutation tail probability of the quadratic form
#' `Q = (T - mu)' M^{-1} (T - mu)`, where `T_c = sum_i treatment_i * score_{i,c}`
#' is a vector of linear statistics, `mu` is its exact permutation mean, and `M`
#' is a metric. This is the omnibus covariate-balance statistic of Hansen and
#' Bowers (2008) when `M` is the permutation covariance of `T`. The permutation
#' null fixes each stratum's treated count.
#'
#' The metric enters only through the rotation `R` with `R R' = M^{-1}`: exactly,
#' `Q = || R'(T - mu) ||^2`, a sum of squares whose Gaussian-approximation weights
#' are the eigenvalues of the rotated permutation covariance.
#'
#' @param scores numeric matrix of per-unit scores (units x representations).
#' @param treatment 0/1 vector; its per-stratum sums fix the treated counts.
#' @param strata factor or vector of stratum labels, the same length as the rows
#'   of `scores`.
#' @param metric either the string `"cov"` (use the exact permutation covariance
#'   of `T`, giving the Hansen-Bowers omnibus statistic) or a symmetric
#'   positive-definite `p x p` matrix (for example a covariate covariance).
#' @param method `"gaussian"` (the weighted-chi-square tail under a Gaussian
#'   approximation to the centred statistic, inverted by Lugannani-Rice) or
#'   `"saddlepoint"` (the non-normal permutation saddlepoint; not yet
#'   implemented).
#' @return a list with `p.value`, the observed `statistic` (`Q`), the `rank` of
#'   the metric, the chi-square `eigenvalues` (weights), the permutation `mean`
#'   of `Q`, the resolved `metric`, and the `method`.
#' @references Hansen, B. B. and Bowers, J. (2008). Covariate balance in simple,
#'   stratified and clustered comparative studies. Statistical Science 23,
#'   219--236.
#' @examples
#' scores <- cbind(c(1, 2, 3, 4, 5, 6), c(1, 1, 2, 2, 5, 9))
#' z <- c(0, 0, 0, 1, 1, 1)
#' # omnibus balance test (metric = the permutation covariance)
#' fastperm_spa_quadratic(scores, z, rep(1, 6))$p.value
#' @export
fastperm_spa_quadratic <- function(scores, treatment, strata,
                                   metric = "cov",
                                   method = c("gaussian", "saddlepoint")) {
  method <- match.arg(method)
  if (identical(method, "saddlepoint"))
    stop("the 'saddlepoint' method (M2, the non-normal permutation saddlepoint) ",
         "is not yet implemented; use method = \"gaussian\".", call. = FALSE)

  scores <- as.matrix(scores)
  p <- ncol(scores)

  ## the joint CGF carries the exact permutation mean vector and covariance V_d
  cg <- fastperm_linear_cgf_mv(scores, treatment, strata)
  mu <- cg$mean; Vd <- cg$cov

  ## resolve the metric M
  if (is.character(metric) && length(metric) == 1L && metric == "cov") {
    M <- Vd
  } else if (is.matrix(metric) && nrow(metric) == p && ncol(metric) == p) {
    M <- metric
  } else {
    stop("`metric` must be \"cov\" or a ", p, " x ", p, " matrix.", call. = FALSE)
  }

  ## rotate by M^{-1/2}: Q = ||rot' (T - mu)||^2 and the Gaussian-approximation
  ## weights are the eigenvalues of the rotated permutation covariance
  rot     <- .quad_inv_sqrt(M)
  Sigma_y <- crossprod(rot, Vd %*% rot)                  # rot' V_d rot (r x r)
  lambda  <- Re(eigen((Sigma_y + t(Sigma_y)) / 2, symmetric = TRUE,
                      only.values = TRUE)$values)
  lambda  <- lambda[lambda > sqrt(.Machine$double.eps) * max(lambda, 1)]

  t_obs <- as.numeric(crossprod(scores, as.numeric(treatment)))   # T at observed z
  d_y   <- as.numeric(crossprod(rot, t_obs - mu))
  q_obs <- sum(d_y^2)

  list(p.value     = .quad_lr_upper(lambda, q_obs),
       statistic   = q_obs,
       rank        = length(lambda),
       eigenvalues = lambda,
       mean        = sum(lambda),
       metric      = M,
       method      = method)
}
