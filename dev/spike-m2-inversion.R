## spike-m2-inversion.R --- M2 prototype: the full non-Gaussian permutation
## quadratic-form UPPER tail. K_Q(theta) by Gauss-Hermite of the Hubbard-
## Stratonovich integral (theta > 0 only -- the p-value region), theta-derivatives
## by finite difference, then Lugannani-Rice + Barndorff-Nielsen r*. Validated
## against EXACT enumeration of the orbit and compared to M1 (Gaussian d) and imhof.
## SCRATCH (dev/).
##
## Two questions, answered in order:
##  (1) how accurate is the GH evaluation of K_Q (self-convergence n=40/64, and vs
##      the exact log-mean-exp over the orbit)? this sets whether FD derivatives
##      are usable.
##  (2) does M2 (FD on GH K_Q + LR/r*) match the exact tail and beat M1?

gauss_hermite <- function(n) {
  i <- seq_len(n - 1); a <- sqrt(i / 2)
  J <- matrix(0, n, n); J[cbind(i, i + 1)] <- a; J[cbind(i + 1, i)] <- a
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = (e$vectors[1, o])^2 * sqrt(pi))
}
tensor_normal <- function(n, r) {
  gh <- gauss_hermite(n); z1 <- sqrt(2) * gh$x; lw1 <- log(gh$w) - 0.5 * log(pi)
  grid <- as.matrix(expand.grid(rep(list(seq_len(n)), r)))
  list(Z = matrix(z1[grid], nrow(grid), r),
       lom = rowSums(matrix(lw1[grid], nrow(grid), r)))
}
logsumexp <- function(l) { m <- max(l); m + log(sum(exp(l - m))) }

## K_Q(theta) = log E_W[exp(K_e(sqrt(2 theta) W))], theta > 0; NA otherwise.
make_KQ <- function(Kval, tns) function(theta) {
  if (theta <= 0) return(NA_real_)
  b <- sqrt(2 * theta)
  logsumexp(tns$lom + apply(tns$Z, 1, function(zz) Kval(b * zz)))
}
fd1 <- function(f, x, h) (f(x + h) - f(x - h)) / (2 * h)
fd2 <- function(f, x, h) (f(x + h) - 2 * f(x) + f(x - h)) / (h^2)

## upper-tail saddlepoint (q > E[Q] => theta > 0). h1, h2 separate FD steps.
lr_tail_upper <- function(KQ, q, m0, qmax, h1 = 1e-4, h2 = 2e-3) {
  KQp  <- function(th) fd1(KQ, th, h1)
  KQpp <- function(th) fd2(KQ, th, h2)
  if (q >= qmax) return(c(lr = 0, rstar = 0))
  lo <- h2 * 2; hi <- lo                                # keep theta - h2 > 0
  if (KQp(lo) >= q) return(c(lr = NA, rstar = NA))      # saddle below FD floor
  while (KQp(hi) < q && hi < 1e8) { lo <- hi; hi <- hi * 2 }
  th <- stats::uniroot(function(t) KQp(t) - q, c(lo, hi))$root
  w  <- sign(th) * sqrt(2 * (th * q - KQ(th)))
  u  <- th * sqrt(KQpp(th))
  lr <- stats::pnorm(w, lower.tail = FALSE) + stats::dnorm(w) * (1 / u - 1 / w)
  rs <- stats::pnorm(w + log(u / w) / w, lower.tail = FALSE)
  c(lr = min(1, max(0, lr)), rstar = min(1, max(0, rs)), theta = th)
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

suppressMessages(devtools::load_all(".", quiet = TRUE))
sc <- cbind(rep(c(1, 2), 10) + rep(0:9, each = 2) * 0.3,
            c(rbind(rep(1, 10), c(2, 3, 5, 8, 13, 21, 34, 55, 89, 144))))
tr <- rep(c(0, 1), 10); st <- rep(1:10, each = 2)

for (mlab in c("Vd", "Sigmax")) {
  M   <- if (mlab == "Vd") fastperm_linear_cgf_mv(sc, tr, st)$cov else stats::cov(sc)
  rot <- metric_sqrt_inv(M); scy <- sc %*% rot
  cg  <- fastperm_linear_cgf_mv(scy, tr, st); muy <- cg$mean
  Kval <- function(s) cg$cgf(s) - sum(s * muy)
  lam  <- Re(eigen((cg$cov + t(cg$cov)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  lam  <- lam[lam > sqrt(.Machine$double.eps) * max(lam, 1)]
  Z   <- enum_orbit(tr, st); Qen <- colSums((t(scy) %*% Z - muy)^2)
  m0  <- mean(Qen); qmax <- max(Qen)
  KQ40 <- make_KQ(Kval, tensor_normal(40, 2))
  KQ   <- make_KQ(Kval, tensor_normal(64, 2))
  KQex <- function(th) { mx <- max(th * Qen); mx + log(mean(exp(th * Qen - mx))) }

  cat(sprintf("\n=== metric = %s | E[Q]=%.4f qmax=%.4f lambda=(%s) ===\n",
              mlab, m0, qmax, paste(sprintf("%.3f", lam), collapse = ",")))
  cat("(1) K_Q accuracy: GH64 vs GH40 vs exact log-mean-exp\n")
  for (th in c(0.05, 0.2, 0.5, 1.0)) {
    cat(sprintf("    theta=%.2f  GH64=%.8f  GH40-GH64=%.1e  exact-GH64=%.1e\n",
                th, KQ(th), KQ40(th) - KQ(th), KQex(th) - KQ(th)))
  }
  cat("(2) upper-tail probabilities\n")
  cat(sprintf("    %9s %9s %9s %9s %9s %9s\n", "q", "exact", "M2_LR", "M2_r*", "M1_gauss", "imhof"))
  for (p in c(0.75, 0.9, 0.95, 0.99, 0.999)) {
    q  <- as.numeric(stats::quantile(Qen, p))
    ex <- mean(Qen >= q)
    m2 <- lr_tail_upper(KQ, q, m0, qmax)
    m1 <- .quad_lr_upper(lam, q)
    im <- if (requireNamespace("CompQuadForm", quietly = TRUE)) CompQuadForm::imhof(q, lambda = lam)$Qq else NA
    cat(sprintf("    %9.4f %9.5f %9.5f %9.5f %9.5f %9.5f\n", q, ex, m2["lr"], m2["rstar"], m1, im))
  }
}
