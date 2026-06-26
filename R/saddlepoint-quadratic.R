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

## ---- M2: the non-Gaussian permutation saddlepoint --------------------------
## K_Q(theta) = log E_{W ~ N(0, I_r)}[exp(K_e(sqrt(2 theta) W))], the exact CGF of
## Q = ||e||^2 via the Hubbard-Stratonovich identity (theta > 0). K_e is the
## centred joint CGF of the rotated, centred statistic e = T_y - mu_y. The
## Gaussian expectation is evaluated by tensor Gauss-Hermite quadrature; the
## resulting exact K_Q is inverted by Lugannani-Rice (and Barndorff-Nielsen r*).

## physicists' Gauss-Hermite nodes/weights (Golub-Welsch): the standard rule for a
## Gaussian expectation, vendored to keep the runtime MIT-licensed.
.gh_nodes <- function(n) {
  i <- seq_len(n - 1); a <- sqrt(i / 2)
  J <- matrix(0, n, n); J[cbind(i, i + 1)] <- a; J[cbind(i + 1, i)] <- a
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = (e$vectors[1, o])^2 * sqrt(pi))
}

## tensor N(0, I_r) rule: E[f(W)] ~ sum_k exp(lw_k) f(Z[k, ]). Physicists' nodes x
## map to z = sqrt(2) x, weights w to log(w) - r/2 log(pi).
.tensor_normal <- function(n, r) {
  gh <- .gh_nodes(n); z1 <- sqrt(2) * gh$x; lw1 <- log(gh$w) - 0.5 * log(pi)
  grid <- as.matrix(expand.grid(rep(list(seq_len(n)), r)))
  list(Z = matrix(z1[grid], nrow(grid), r),
       lw = rowSums(matrix(lw1[grid], nrow(grid), r)))
}

## .log_esym, vectorised across columns: Wb is n_b x N (tilted weights, one column
## per quadrature node), returns the length-N vector of log e_m(exp(column)). Same
## O(n_b * m) log-space recursion as .log_esym, run on all nodes at once.
.log_esym_cols <- function(Wb, m) {
  N <- ncol(Wb)
  logp <- matrix(-Inf, m + 1L, N); logp[1L, ] <- 0           # log e_0 = 0
  for (i in seq_len(nrow(Wb))) {
    shifted <- logp[seq_len(m), , drop = FALSE] + rep(Wb[i, ], each = m)
    logp[2L:(m + 1L), ] <- .logaddexp(logp[2L:(m + 1L), , drop = FALSE], shifted)
  }
  logp[m + 1L, ]
}

## K_Q(theta) closure: the centred CGF of the rotated statistic evaluated at every
## node (vectorised), then the log of the Gauss-Hermite-weighted sum.
.make_KQ_quad <- function(scoresY, idx, mb, nb, muy, tns) {
  Z <- tns$Z; lw <- tns$lw
  function(theta) {
    if (theta <= 0) return(NA_real_)
    S    <- sqrt(2 * theta) * Z                             # nodes x r tilt vectors
    Wall <- scoresY %*% t(S)                                # n x N projected weights
    KT   <- numeric(nrow(Z))
    for (b in seq_along(idx)) {
      m <- mb[b]; if (m == 0L) next
      ix <- idx[[b]]; nbb <- nb[b]
      if (m == nbb) { KT <- KT + colSums(Wall[ix, , drop = FALSE]); next }
      KT <- KT + .log_esym_cols(Wall[ix, , drop = FALSE], m) - lchoose(nbb, m)
    }
    Ke <- KT - as.numeric(S %*% muy)                        # centre: K_e = K_T - s'mu_y
    mx <- max(lw + Ke); mx + log(sum(exp(lw + Ke - mx)))    # logsumexp
  }
}

## upper-tail P(Q >= q) from the numeric CGF K_Q (q > E[Q] => saddle theta > 0; Q
## has bounded support, so K_Q has no pole). Derivatives by finite difference (GH
## is near machine precision). Returns Lugannani-Rice and Barndorff-Nielsen r*.
## Below the mean, and in the narrow band where the saddle falls under the FD
## floor, falls back to the Gaussian-d tail (accurate there, not a rejection
## region). h1, h2 are the first- and second-derivative FD steps.
.quad_spa_upper <- function(KQ, q, m0, lambda, h1 = 1e-4, h2 = 2e-3,
                            theta_max = 200) {
  gauss <- .quad_lr_upper(lambda, q)
  if (q <= 0) return(c(lr = 1, rstar = 1))
  if (q <= m0) return(c(lr = gauss, rstar = gauss))
  KQp  <- function(th) (KQ(th + h1) - KQ(th - h1)) / (2 * h1)
  KQpp <- function(th) (KQ(th + h2) - 2 * KQ(th) + KQ(th - h2)) / h2^2
  lo <- 2 * h2                                             # keep theta - h2 > 0
  if (KQp(lo) >= q) return(c(lr = gauss, rstar = gauss))   # saddle under FD floor
  hi <- lo
  while (KQp(hi) < q && hi < theta_max) { lo <- hi; hi <- hi * 2 }
  if (KQp(hi) < q) {                                       # q at/beyond support max
    warning("observed Q is at or beyond the saddlepoint's reliable range; ",
            "falling back to the Gaussian-d tail.", call. = FALSE)
    return(c(lr = gauss, rstar = gauss))
  }
  th  <- stats::uniroot(function(t) KQp(t) - q, c(lo, hi))$root
  kpp <- KQpp(th)
  if (!is.finite(th) || kpp <= 0) {                        # FD noise (coarse/extreme)
    warning("the saddlepoint second derivative is unreliable here (coarse lattice ",
            "or extreme tail); falling back to the Gaussian-d tail.", call. = FALSE)
    return(c(lr = gauss, rstar = gauss))
  }
  w  <- sign(th) * sqrt(2 * (th * q - KQ(th)))
  u  <- th * sqrt(kpp)
  lr <- stats::pnorm(w, lower.tail = FALSE) + stats::dnorm(w) * (1 / u - 1 / w)
  rs <- stats::pnorm(w + log(u / w) / w, lower.tail = FALSE)
  if (!is.finite(lr) || !is.finite(rs)) return(c(lr = gauss, rstar = gauss))
  c(lr = min(1, max(0, lr)), rstar = min(1, max(0, rs)))
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
#' @param method `"gaussian"` (M1: the weighted-chi-square tail under a Gaussian
#'   approximation to the centred statistic, inverted by Lugannani-Rice) or
#'   `"saddlepoint"` (M2: the non-normal permutation saddlepoint, which carries the
#'   true skewness and higher cumulants of the permutation law into the tail).
#' @param nodes number of Gauss-Hermite nodes per dimension for the
#'   `"saddlepoint"` method (the Hubbard-Stratonovich integral is an `r`-dimensional
#'   Gaussian expectation evaluated by tensor quadrature). Ignored for `"gaussian"`.
#' @return a list with `p.value`, the observed `statistic` (`Q`), the `rank` of
#'   the metric, the chi-square `eigenvalues` (weights), the permutation `mean`
#'   of `Q`, the resolved `metric`, and the `method`. For `method = "saddlepoint"`
#'   it also carries `p.value.rstar`, the Barndorff-Nielsen `r*` higher-order
#'   refinement of `p.value`.
#' @details The `"saddlepoint"` method uses the exact permutation CGF of the
#'   centred statistic via the Hubbard-Stratonovich identity
#'   `M_Q(theta) = E_{W ~ N(0, I_r)}[exp(K_e(sqrt(2 theta) W))]`, evaluated by
#'   tensor Gauss-Hermite quadrature and inverted by Lugannani-Rice. The tensor
#'   rule costs `nodes^r`, so it is feasible only for small `r` (the number of
#'   representations); larger `r` needs the sparse-grid/QMC engine, which is not
#'   yet implemented. Below `E[Q]` the method falls back to the Gaussian-d tail,
#'   which is accurate there and not a rejection region.
#' @references Hansen, B. B. and Bowers, J. (2008). Covariate balance in simple,
#'   stratified and clustered comparative studies. Statistical Science 23,
#'   219--236.
#' @examples
#' scores <- cbind(c(1, 2, 3, 4, 5, 6), c(1, 1, 2, 2, 5, 9))
#' z <- c(0, 0, 0, 1, 1, 1)
#' # omnibus balance test (metric = the permutation covariance)
#' fastperm_spa_quadratic(scores, z, rep(1, 6))$p.value
#' # the non-normal permutation saddlepoint
#' fastperm_spa_quadratic(scores, z, rep(1, 6), method = "saddlepoint")$p.value
#' @export
fastperm_spa_quadratic <- function(scores, treatment, strata,
                                   metric = "cov",
                                   method = c("gaussian", "saddlepoint"),
                                   nodes = 32L) {
  method <- match.arg(method)
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
  ## weights are the eigenvalues of the rotated permutation covariance. Keep the
  ## eigenVECTORS too: when V_d (hence Sigma_y) is rank-deficient -- collinear
  ## representations -- the centred statistic lives in the r-dimensional range of
  ## Sigma_y, and rotating further onto that range (rotr, p x r) makes Q a sum of
  ## exactly r coordinates, which the M2 quadrature dimension must match. For a
  ## full-rank metric the extra rotation is orthogonal and leaves Q unchanged.
  rot     <- .quad_inv_sqrt(M)
  Sigma_y <- crossprod(rot, Vd %*% rot)                  # rot' V_d rot (p x p)
  eig     <- eigen((Sigma_y + t(Sigma_y)) / 2, symmetric = TRUE)
  keep    <- eig$values > sqrt(.Machine$double.eps) * max(eig$values, 1)
  lambda  <- eig$values[keep]
  r       <- length(lambda)
  rotr    <- rot %*% eig$vectors[, keep, drop = FALSE]   # p x r combined rotation

  t_obs <- as.numeric(crossprod(scores, as.numeric(treatment)))   # T at observed z
  d_y   <- as.numeric(crossprod(rotr, t_obs - mu))       # r range coordinates
  q_obs <- sum(d_y^2)

  out <- list(statistic = q_obs, rank = r, eigenvalues = lambda,
              mean = sum(lambda), metric = M, method = method)

  if (method == "gaussian") {
    out$p.value <- .quad_lr_upper(lambda, q_obs)
    return(out)
  }

  ## M2: tensor Gauss-Hermite over the r-dimensional Hubbard-Stratonovich integral
  if (nodes^r > 250000L)
    stop("rank r = ", r, " is too large for the tensor Gauss-Hermite engine ",
         "(", nodes, "^", r, " nodes); the sparse-grid/QMC engine for larger r is ",
         "not yet implemented. Use method = \"gaussian\", or reduce the number of ",
         "representations.", call. = FALSE)
  scoresY <- scores %*% rotr                             # n x r rotated scores
  muy     <- as.numeric(crossprod(rotr, mu))             # rotr' mu = E[T_y] in range
  idx     <- split(seq_len(nrow(scores)), as.factor(strata))
  tcounts <- as.numeric(treatment)
  mb      <- unname(vapply(idx, function(ix) sum(tcounts[ix]), numeric(1)))
  nb      <- unname(lengths(idx))
  KQ      <- .make_KQ_quad(scoresY, idx, mb, nb, muy, .tensor_normal(nodes, r))
  tail    <- .quad_spa_upper(KQ, q_obs, sum(lambda), lambda)
  out$p.value       <- unname(tail["lr"])
  out$p.value.rstar <- unname(tail["rstar"])
  out
}
