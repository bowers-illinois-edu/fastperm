## helper-quadratic-enumerate.R --- exact enumeration oracle for the permutation
## distribution of the quadratic form Q = (T - mu)' M^{-1} (T - mu). testthat
## auto-loads helper-*.R, so these are available to every test file.
##
## These functions are the ground truth for the quadratic-form saddlepoint: they
## walk the entire within-stratum permutation orbit, so the tail probabilities
## they return are exact (no Monte Carlo error) for the small designs the tests
## use. They are the same construction as dev/spike-quadratic-spa.R, kept here as
## reusable test infrastructure.

## M^{-1/2} via the spectral decomposition, dropping numerically-zero eigenvalues
## so a rank-deficient metric still yields a rotation onto its range (the same
## convention RItools' XtX_pseudoinv_sqrt uses).
.metric_sqrt_inv <- function(M) {
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  keep <- e$values > sqrt(.Machine$double.eps) * max(e$values, 1)
  e$vectors[, keep, drop = FALSE] %*%
    diag(1 / sqrt(e$values[keep]), sum(keep), sum(keep))
}

## every treatment vector sharing `treatment`'s per-stratum treated counts, as
## columns of an n x N matrix (N = prod_b choose(n_b, m_b)).
.enumerate_orbit <- function(treatment, strata) {
  idx <- split(seq_along(treatment), as.factor(strata))
  mb  <- vapply(idx, function(ix) sum(treatment[ix]), numeric(1))
  per <- lapply(seq_along(idx), function(b)
    utils::combn(length(idx[[b]]), mb[b], simplify = FALSE))
  grid <- expand.grid(lapply(per, seq_along))
  Z <- matrix(0, length(treatment), nrow(grid))
  for (k in seq_len(nrow(grid))) {
    z <- numeric(length(treatment))
    for (b in seq_along(idx))
      z[idx[[b]][per[[b]][[grid[k, b]]]]] <- 1
    Z[, k] <- z
  }
  Z
}

## exact permutation law of Q under the metric M. Returns Q computed two ways
## (directly as d' M^{-1} d and via the rotation as ||rot' d||^2, so a test can
## confirm the rotation reduction), the exact mean vector mu and permutation
## covariance Sigma of T, and the observed Q at the observed treatment vector.
enumerate_Q <- function(scores, treatment, strata, M) {
  scores <- as.matrix(scores)
  Z   <- .enumerate_orbit(treatment, strata)
  Tm  <- t(scores) %*% Z                       # p x N: T for each orbit member
  mu  <- rowMeans(Tm)                          # exact E[T] over the orbit
  D   <- Tm - mu
  A   <- chol2inv(chol((M + t(M)) / 2))        # M^{-1}
  rot <- .metric_sqrt_inv(M)
  d0  <- as.numeric(t(scores) %*% as.numeric(treatment)) - mu
  list(Q     = colSums(D * (A %*% D)),         # d' M^{-1} d, directly
       Q_rot = colSums((t(rot) %*% D)^2),      # ||rot' d||^2
       mu    = mu,
       Sigma = (D %*% t(D)) / ncol(D),         # exact permutation covariance of T
       q_obs = as.numeric(t(d0) %*% A %*% d0))
}

## exact upper-tail probability P(Q >= q) under the permutation null
enum_quadratic_tail <- function(en, q) mean(en$Q >= q - 1e-9)
