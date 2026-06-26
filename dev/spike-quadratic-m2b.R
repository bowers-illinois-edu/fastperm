## spike-quadratic-m2b.R --- M2 re-prototype with ANALYTIC CGF derivatives, on a
## genuinely non-Gaussian design (matched pairs). SCRATCH (dev/).
##
## Fixes the noise that sank dev/spike-quadratic-m2.R: the inner Laplace now uses
## the analytic gradient of K_dy, and the determinant term the analytic Hessian
## (dev/spike-cgf-derivs.R), so K_Q(theta) is clean and its theta-derivatives for
## Lugannani-Rice are reliable. The test design is matched pairs, where d_y is far
## from Gaussian and M1 (Gaussian-d weighted chi-square) should visibly miss.

## --- analytic CGF derivatives (validated in spike-cgf-derivs.R) --------------
.loge_safe <- function(x, m) {
  if (m < 0 || m > length(x)) return(-Inf)
  if (m == 0) return(0)
  .log_esym(x, m)
}
.incl_probs <- function(w, m) {
  n <- length(w); loge_m <- .loge_safe(w, m)
  p <- vapply(seq_len(n),
              function(i) exp(w[i] + .loge_safe(w[-i], m - 1) - loge_m), numeric(1))
  P <- matrix(0, n, n)
  if (n >= 2) for (i in 1:(n - 1)) for (j in (i + 1):n) {
    pij <- exp(w[i] + w[j] + .loge_safe(w[-c(i, j)], m - 2) - loge_m)
    P[i, j] <- pij; P[j, i] <- pij
  }
  diag(P) <- p; list(p = p, P = P)
}
cgf_mv_derivs <- function(scores, treatment, strata, s) {
  scores <- as.matrix(scores); d <- ncol(scores)
  idx <- split(seq_along(treatment), as.factor(strata))
  grad <- numeric(d); hess <- matrix(0, d, d)
  for (ix in idx) {
    Vb <- scores[ix, , drop = FALSE]; m <- sum(treatment[ix]); nb <- length(ix)
    if (m == 0) next
    if (m == nb) { grad <- grad + colSums(Vb); next }
    ip <- .incl_probs(as.numeric(Vb %*% s), m)
    grad <- grad + as.numeric(crossprod(Vb, ip$p))
    hess <- hess + crossprod(Vb, (ip$P - tcrossprod(ip$p)) %*% Vb)
  }
  list(grad = grad, hess = hess)
}

## --- M2 inversion with analytic inner pieces ---------------------------------
## K_Q(theta) = g(w*) - 1/2 log det(I - 2 theta Hess K_dy(s*)),
## g(w) = -||w||^2/2 + K_dy(sqrt(2 theta) w), maximised with the analytic gradient.
make_KQ <- function(Kval, Kgrad, Khess, r) {
  function(theta) {
    if (theta <= 0) return(0)
    a   <- sqrt(2 * theta)
    opt <- stats::optim(rep(0, r),
                        fn = function(w) -(-0.5 * sum(w^2) + Kval(a * w)),
                        gr = function(w) w - a * Kgrad(a * w),
                        method = "BFGS", control = list(reltol = 1e-12))
    ws <- opt$par; s_star <- a * ws
    ld <- determinant(diag(r) - 2 * theta * Khess(s_star), logarithm = TRUE)
    if (ld$sign <= 0) return(NA_real_)
    (-0.5 * sum(ws^2) + Kval(s_star)) - 0.5 * as.numeric(ld$modulus)
  }
}

## Lugannani-Rice upper tail from a (now clean) K_Q; theta-derivatives by finite
## differences with a step tuned to the residual ~1e-8 noise in K_Q.
lr_upper <- function(KQ, lam, q, ht = 0.02) {
  m0 <- sum(lam); if (q <= m0) return(NA_real_)
  Kp  <- function(th) (KQ(th + ht) - KQ(th - ht)) / (2 * ht)
  Kpp <- function(th) (KQ(th + ht) - 2 * KQ(th) + KQ(th - ht)) / ht^2
  cap  <- 1 / (2 * max(lam))
  grid <- seq(cap * 0.02, cap * 0.98, length.out = 50)
  kp   <- vapply(grid, function(th) tryCatch(Kp(th), error = function(e) NA), numeric(1))
  ok   <- which(is.finite(kp)); cr <- ok[which(kp[ok] >= q)[1]]
  if (is.na(cr) || cr == ok[1]) return(NA_real_)
  th <- stats::uniroot(function(t) Kp(t) - q, c(grid[cr - 1], grid[cr]))$root
  w  <- sign(th) * sqrt(2 * (th * q - KQ(th))); u <- th * sqrt(Kpp(th))
  min(1, max(0, stats::pnorm(w, lower.tail = FALSE) + stats::dnorm(w) * (1/u - 1/w)))
}

## --- M1 and enumeration, as before -------------------------------------------
metric_sqrt_inv <- function(M) {
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  keep <- e$values > sqrt(.Machine$double.eps) * max(e$values, 1)
  e$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(e$values[keep]), sum(keep), sum(keep))
}
enum_orbit <- function(treatment, strata) {
  idx <- split(seq_along(treatment), as.factor(strata))
  mb  <- vapply(idx, function(ix) sum(treatment[ix]), numeric(1))
  per <- lapply(seq_along(idx), function(b) utils::combn(length(idx[[b]]), mb[b], simplify = FALSE))
  grid <- expand.grid(lapply(per, seq_along)); Z <- matrix(0, length(treatment), nrow(grid))
  for (k in seq_len(nrow(grid))) { z <- numeric(length(treatment))
    for (b in seq_along(idx)) z[idx[[b]][per[[b]][[grid[k, b]]]]] <- 1; Z[, k] <- z }
  Z
}
lr_upper_m1 <- function(lam, q) {
  m0 <- sum(lam); if (q <= m0) return(NA_real_)
  K <- function(th) -0.5 * sum(log(1 - 2 * th * lam))
  Kp <- function(th) sum(lam / (1 - 2 * th * lam)); Kpp <- function(th) sum(2 * lam^2 / (1 - 2 * th * lam)^2)
  pole <- 1 / (2 * max(lam)); th <- stats::uniroot(function(t) Kp(t) - q, c(1e-12, pole * (1 - 1e-12)))$root
  w <- sign(th) * sqrt(2 * (th * q - K(th))); u <- th * sqrt(Kpp(th))
  stats::pnorm(w, lower.tail = FALSE) + stats::dnorm(w) * (1/u - 1/w)
}
enum_tail <- function(Q, q) mean(Q >= q - 1e-9)

## --- 10 matched pairs, skewed second covariate (d_y far from Gaussian) -------
set.seed(1)
sc <- cbind(rep(c(1, 2), 10) + rep(0:9, each = 2) * 0.3,
            c(rbind(rep(1, 10), c(2, 3, 5, 8, 13, 21, 34, 55, 89, 144))))
tr <- rep(c(0, 1), 10)
st <- rep(1:10, each = 2)

for (mlab in c("Vd", "Sigmax")) {
  M <- if (mlab == "Vd") fastperm_linear_cgf_mv(sc, tr, st)$cov else stats::cov(sc)
  rot <- metric_sqrt_inv(M); scy <- sc %*% rot
  cg  <- fastperm_linear_cgf_mv(scy, tr, st); muy <- cg$mean
  Kval  <- function(s) cg$cgf(s) - sum(s * muy)
  Kgrad <- function(s) cgf_mv_derivs(scy, tr, st, s)$grad - muy
  Khess <- function(s) cgf_mv_derivs(scy, tr, st, s)$hess
  lam <- Re(eigen((cg$cov + t(cg$cov)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  lam <- lam[lam > sqrt(.Machine$double.eps) * max(lam, 1)]
  KQ  <- make_KQ(Kval, Kgrad, Khess, length(muy))

  Z <- enum_orbit(tr, st); Dy <- t(scy) %*% Z - muy; Qen <- colSums(Dy^2)
  cat(sprintf("\n=== 10 pairs | metric = %s | r=%d, E[Q]=%.3f ===\n", mlab, length(muy), sum(lam)))
  cat(sprintf("lambda = %s\n", paste(sprintf("%.3f", lam), collapse = ", ")))
  qs <- stats::quantile(Qen[Qen > sum(lam)], c(.5, .75, .9, .97), names = FALSE)
  cat(sprintf("%9s %12s %12s %12s %14s %14s\n", "q", "enum", "M1", "M2", "M1 rel.err", "M2 rel.err"))
  for (q in qs) {
    e <- enum_tail(Qen, q); m1 <- lr_upper_m1(lam, q); m2 <- lr_upper(KQ, lam, q)
    cat(sprintf("%9.4f %12.5f %12.5f %12.5f %13.1f%% %13.1f%%\n",
                q, e, m1, m2, 100 * (m1 - e) / e, 100 * (m2 - e) / e))
  }
}
