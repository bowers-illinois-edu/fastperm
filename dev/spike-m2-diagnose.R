## spike-m2-diagnose.R --- isolate where M2 fails: does the Hubbard-Stratonovich
## + Laplace K_Q(theta) match the EXACT K_Q(theta) = log E[exp(theta Q)] from
## enumeration? If H-S tracks the Gaussian curve but not the exact one, the inner
## Laplace is too crude to carry the non-normality; if H-S tracks exact, the fault
## is downstream in Lugannani-Rice. SCRATCH (dev/).

.loge_safe <- function(x, m) { if (m < 0 || m > length(x)) return(-Inf); if (m == 0) return(0); .log_esym(x, m) }
.incl_probs <- function(w, m) {
  n <- length(w); loge_m <- .loge_safe(w, m)
  p <- vapply(seq_len(n), function(i) exp(w[i] + .loge_safe(w[-i], m - 1) - loge_m), numeric(1))
  P <- matrix(0, n, n)
  if (n >= 2) for (i in 1:(n - 1)) for (j in (i + 1):n) {
    pij <- exp(w[i] + w[j] + .loge_safe(w[-c(i, j)], m - 2) - loge_m); P[i, j] <- pij; P[j, i] <- pij }
  diag(P) <- p; list(p = p, P = P)
}
cgf_mv_derivs <- function(scores, treatment, strata, s) {
  scores <- as.matrix(scores); d <- ncol(scores); idx <- split(seq_along(treatment), as.factor(strata))
  grad <- numeric(d); hess <- matrix(0, d, d)
  for (ix in idx) { Vb <- scores[ix, , drop = FALSE]; m <- sum(treatment[ix]); nb <- length(ix)
    if (m == 0) next; if (m == nb) { grad <- grad + colSums(Vb); next }
    ip <- .incl_probs(as.numeric(Vb %*% s), m)
    grad <- grad + as.numeric(crossprod(Vb, ip$p)); hess <- hess + crossprod(Vb, (ip$P - tcrossprod(ip$p)) %*% Vb) }
  list(grad = grad, hess = hess)
}
metric_sqrt_inv <- function(M) { e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  keep <- e$values > sqrt(.Machine$double.eps) * max(e$values, 1)
  e$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(e$values[keep]), sum(keep), sum(keep)) }
enum_orbit <- function(treatment, strata) {
  idx <- split(seq_along(treatment), as.factor(strata)); mb <- vapply(idx, function(ix) sum(treatment[ix]), numeric(1))
  per <- lapply(seq_along(idx), function(b) utils::combn(length(idx[[b]]), mb[b], simplify = FALSE))
  grid <- expand.grid(lapply(per, seq_along)); Z <- matrix(0, length(treatment), nrow(grid))
  for (k in seq_len(nrow(grid))) { z <- numeric(length(treatment))
    for (b in seq_along(idx)) z[idx[[b]][per[[b]][[grid[k, b]]]]] <- 1; Z[, k] <- z }; Z
}
KQ_hs <- function(theta, Kval, Kgrad, Khess, r) {
  if (theta <= 0) return(0); a <- sqrt(2 * theta)
  opt <- stats::optim(rep(0, r), fn = function(w) -(-0.5 * sum(w^2) + Kval(a * w)),
                      gr = function(w) w - a * Kgrad(a * w), method = "BFGS", control = list(reltol = 1e-12))
  ws <- opt$par; s_star <- a * ws; ld <- determinant(diag(r) - 2 * theta * Khess(s_star), logarithm = TRUE)
  if (ld$sign <= 0) return(NA_real_); (-0.5 * sum(ws^2) + Kval(s_star)) - 0.5 * as.numeric(ld$modulus)
}

## 10 matched pairs, Vd metric
sc <- cbind(rep(c(1, 2), 10) + rep(0:9, each = 2) * 0.3,
            c(rbind(rep(1, 10), c(2, 3, 5, 8, 13, 21, 34, 55, 89, 144))))
tr <- rep(c(0, 1), 10); st <- rep(1:10, each = 2)
M   <- fastperm_linear_cgf_mv(sc, tr, st)$cov
rot <- metric_sqrt_inv(M); scy <- sc %*% rot
cg  <- fastperm_linear_cgf_mv(scy, tr, st); muy <- cg$mean
Kval  <- function(s) cg$cgf(s) - sum(s * muy)
Kgrad <- function(s) cgf_mv_derivs(scy, tr, st, s)$grad - muy
Khess <- function(s) cgf_mv_derivs(scy, tr, st, s)$hess
lam <- Re(eigen((cg$cov + t(cg$cov)) / 2, symmetric = TRUE, only.values = TRUE)$values)
lam <- lam[lam > sqrt(.Machine$double.eps) * max(lam, 1)]

Z <- enum_orbit(tr, st); Qen <- colSums((t(scy) %*% Z - muy)^2)
KQ_exact <- function(theta) { mx <- max(theta * Qen); mx + log(mean(exp(theta * Qen - mx))) }

cap <- 1 / (2 * max(lam))
cat(sprintf("E[Q]=%.4f  lambda=(%.3f, %.3f)  pole(Gauss)=%.4f\n", mean(Qen), lam[1], lam[2], cap))
cat(sprintf("%8s %14s %14s %14s\n", "theta", "K_Q exact", "K_Q H-S", "K_Q Gauss"))
for (th in cap * c(0.1, 0.3, 0.5, 0.7, 0.9)) {
  ex <- KQ_exact(th); hs <- KQ_hs(th, Kval, Kgrad, Khess, 2); gs <- -0.5 * sum(log(1 - 2 * th * lam))
  cat(sprintf("%8.4f %14.6f %14.6f %14.6f\n", th, ex, hs, gs))
}
