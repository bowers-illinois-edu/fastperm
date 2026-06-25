## cgf-stratified-mv.R --- the JOINT (multivariate) CGF of a vector of linear
## statistics under the within-stratum permutation null. The enabler for the
## quadratic-form and max-joint-tail saddlepoints.
##
## For a vector statistic T = (T_1, ..., T_d) with T_c = sum_i z_i score_{i,c},
## the joint moment generating function of the stratum-b contribution is
## E[exp(s' sum_{i in treated} v_i)] where v_i is unit i's score ROW. Writing the
## projected per-unit weight w_i = s' v_i, this is exactly the univariate
## elementary-symmetric MGF on the w_i. So the d-dimensional CGF K_T(s) is the
## 1-D CGF (cgf-stratified.R) evaluated on the projected scores `scores %*% s`,
## and its gradient and Hessian at s = 0 are the exact permutation mean vector and
## covariance. That projection is what lets the quadratic form and the max reuse
## the same per-block machinery instead of needing a new multivariate object.

#' Joint CGF of a vector of linear statistics under within-stratum permutation
#'
#' Builds the joint cumulant generating function `K_T(s)` of the vector statistic
#' `T_c = sum_i treatment_i * score_{i,c}`, referred to the permutation null that
#' fixes each stratum's treated count. `K_T(s)` is evaluated by projecting the
#' score matrix onto `s` and applying the univariate elementary-symmetric CGF, so
#' the cost is that of one 1-D CGF per evaluation. The exact permutation mean
#' vector and covariance are returned in closed form; `grad K_T(0)` and
#' `Hess K_T(0)` reproduce them.
#'
#' @param scores numeric matrix of per-unit scores (units x representations).
#' @param treatment 0/1 vector; its per-stratum sums fix the treated counts.
#' @param strata factor or vector of stratum labels.
#' @return an object of class `fastperm_cgf_mv`: a list with `cgf` (a closure of
#'   the length-d vector `s`), the exact `mean` vector and `cov` matrix, and `d`.
#' @export
fastperm_linear_cgf_mv <- function(scores, treatment, strata) {
  scores <- as.matrix(scores)
  treatment <- as.numeric(treatment); strata <- as.factor(strata)
  n <- nrow(scores); d <- ncol(scores)
  if (length(treatment) != n || length(strata) != n)
    stop("`scores` rows, `treatment`, and `strata` must describe the same units.",
         call. = FALSE)

  idx <- split(seq_len(n), strata)
  mb  <- unname(vapply(idx, function(ix) sum(treatment[ix]), numeric(1)))
  nb  <- unname(lengths(idx))

  ## exact mean vector and covariance (Strasser-Weber), summed over strata
  mu <- numeric(d); Sigma <- matrix(0, d, d)
  for (b in seq_along(idx)) {
    ix <- idx[[b]]; Sb <- scores[ix, , drop = FALSE]; m <- mb[b]; nbb <- nb[b]
    mu <- mu + (m / nbb) * colSums(Sb)
    if (nbb > 1 && m > 0 && m < nbb) {
      w <- m * (nbb - m) / (nbb * (nbb - 1))
      Sc <- sweep(Sb, 2, colMeans(Sb))
      Sigma <- Sigma + w * crossprod(Sc)
    }
  }
  dimnames(Sigma) <- list(colnames(scores), colnames(scores))
  names(mu) <- colnames(scores)

  ## K_T(s): project the scores onto s, then the univariate elementary-symmetric
  ## CGF of the projected weights
  cgf <- function(s) {
    w <- as.numeric(scores %*% s)
    acc <- 0
    for (b in seq_along(idx)) {
      ix <- idx[[b]]; m <- mb[b]; nbb <- nb[b]
      if (m == 0L) next
      if (m == nbb) { acc <- acc + sum(w[ix]); next }   # deterministic stratum
      acc <- acc + .log_esym(w[ix], m) - lchoose(nbb, m)
    }
    acc
  }

  structure(list(cgf = cgf, mean = mu, cov = Sigma, d = d),
            class = "fastperm_cgf_mv")
}

#' @export
print.fastperm_cgf_mv <- function(x, ...) {
  cat(sprintf("<fastperm_cgf_mv> joint CGF of %d linear statistics\n", x$d))
  invisible(x)
}
