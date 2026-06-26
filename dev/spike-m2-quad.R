## spike-m2-quad.R --- the M2 fix: evaluate the Hubbard-Stratonovich integral by
## Gauss-Hermite quadrature instead of Laplace. SCRATCH (dev/).
##
## The H-S identity is exactly an expectation against a standard Gaussian:
##   M_Q(theta) = E_{W ~ N(0, I_r)}[ exp(K_dy(sqrt(2 theta) W)) ].
## Leading-order Laplace around the (centred) mode w = 0 returns the Gaussian
## K_Q and discards all non-normality. Gauss-Hermite quadrature evaluates the
## expectation directly, carrying the full K_dy. This checks K_Q(GH) against the
## exact K_Q = log E[exp(theta Q)] from enumeration. Needs only the cgf closure --
## no derivatives.

## physicists' Gauss-Hermite nodes/weights: INT exp(-x^2) f(x) dx ~ sum w_k f(x_k)
gauss_hermite <- function(n) {
  i <- seq_len(n - 1); a <- sqrt(i / 2)
  J <- matrix(0, n, n); J[cbind(i, i + 1)] <- a; J[cbind(i + 1, i)] <- a
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = (e$vectors[1, o])^2 * sqrt(pi))
}

## M_Q(theta) = E_{N(0,I_2)}[exp(K_dy(sqrt(2 theta) W))] by tensor GH (r = 2).
## With physicists' nodes, W = sqrt(2) x and E_{N(0,1)}[g] = pi^{-1/2} sum w g(sqrt 2 x),
## so the CGF argument is s = sqrt(2 theta) * sqrt(2) * (x_k, x_l) = 2 sqrt(theta) (x_k, x_l).
KQ_gh <- function(theta, Kval, gh) {
  if (theta <= 0) return(0)
  b <- 2 * sqrt(theta); nq <- length(gh$x); acc <- -Inf
  for (k in seq_len(nq)) for (l in seq_len(nq)) {
    lt <- log(gh$w[k]) + log(gh$w[l]) + Kval(c(b * gh$x[k], b * gh$x[l]))
    acc <- max(acc, lt) + log1p(exp(min(acc, lt) - max(acc, lt)))  # logsumexp
  }
  acc - log(pi)                                                    # / pi^{r/2}, r = 2
}

metric_sqrt_inv <- function(M) { e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  keep <- e$values > sqrt(.Machine$double.eps) * max(e$values, 1)
  e$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(e$values[keep]), sum(keep), sum(keep)) }
enum_orbit <- function(treatment, strata) {
  idx <- split(seq_along(treatment), as.factor(strata)); mb <- vapply(idx, function(ix) sum(treatment[ix]), numeric(1))
  per <- lapply(seq_along(idx), function(b) utils::combn(length(idx[[b]]), mb[b], simplify = FALSE))
  grid <- expand.grid(lapply(per, seq_along)); Z <- matrix(0, length(treatment), nrow(grid))
  for (k in seq_len(nrow(grid))) { z <- numeric(length(treatment))
    for (b in seq_along(idx)) z[idx[[b]][per[[b]][[grid[k, b]]]]] <- 1; Z[, k] <- z }; Z }

sc <- cbind(rep(c(1, 2), 10) + rep(0:9, each = 2) * 0.3,
            c(rbind(rep(1, 10), c(2, 3, 5, 8, 13, 21, 34, 55, 89, 144))))
tr <- rep(c(0, 1), 10); st <- rep(1:10, each = 2)
gh <- gauss_hermite(40)

for (mlab in c("Vd", "Sigmax")) {
  M <- if (mlab == "Vd") fastperm_linear_cgf_mv(sc, tr, st)$cov else stats::cov(sc)
  rot <- metric_sqrt_inv(M); scy <- sc %*% rot
  cg  <- fastperm_linear_cgf_mv(scy, tr, st); muy <- cg$mean
  Kval <- function(s) cg$cgf(s) - sum(s * muy)
  lam  <- Re(eigen((cg$cov + t(cg$cov)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  lam  <- lam[lam > sqrt(.Machine$double.eps) * max(lam, 1)]
  Z <- enum_orbit(tr, st); Qen <- colSums((t(scy) %*% Z - muy)^2)
  KQ_exact <- function(th) { mx <- max(th * Qen); mx + log(mean(exp(th * Qen - mx))) }
  cap <- 1 / (2 * max(lam))
  cat(sprintf("\n=== 10 pairs | metric = %s | E[Q]=%.3f, lambda=(%s) ===\n",
              mlab, sum(lam), paste(sprintf("%.2f", lam), collapse = ", ")))
  cat(sprintf("%8s %14s %14s %14s\n", "theta", "K_Q exact", "K_Q GH", "K_Q Gauss"))
  for (th in cap * c(0.1, 0.3, 0.5, 0.7, 0.9)) {
    cat(sprintf("%8.4f %14.6f %14.6f %14.6f\n",
                th, KQ_exact(th), KQ_gh(th, Kval, gh), -0.5 * sum(log(1 - 2 * th * lam))))
  }
}
