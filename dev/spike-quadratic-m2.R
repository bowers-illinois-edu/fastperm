## spike-quadratic-m2.R --- prototype the M2 non-normal permutation saddlepoint
## for Q = ||d_y||^2 and validate it by simulation. SCRATCH (dev/).
##
## Inversion (Hubbard-Stratonovich + Laplace). With K_dy the centred joint CGF of
## the rotated linear-statistic vector d_y (in R^r), and A = I,
##   E[exp(theta ||d_y||^2)]
##     = (2 pi)^{-r/2} INT exp(-||w||^2/2 + K_dy(sqrt(2 theta) w)) dw    (exact),
## then Laplace at the maximiser w* gives
##   K_Q(theta) = g(w*) - 1/2 log det(I_r - 2 theta Hess K_dy(s*)),  s* = sqrt(2 theta) w*,
##   g(w) = -||w||^2/2 + K_dy(sqrt(2 theta) w).
## Lugannani-Rice on K_Q gives the tail.
##
## Two checks:
##  (1) self-check: feed K_dy = (1/2) s' Sigma_y s. Then the integral is exact and
##      K_Q(theta) must equal -1/2 sum log(1 - 2 theta lambda_k) for all theta.
##  (2) real check: feed the actual permutation K_dy and compare the M2 tail to
##      exact enumeration and to M1 (the Gaussian-d weighted-chi-square).
##
## Run:
##   Rscript -e 'devtools::load_all("~/repos/fastperm-route-b"); \
##               source("~/repos/fastperm-route-b/dev/spike-quadratic-m2.R")'

## ---- small numerical-calculus helpers ---------------------------------------
num_grad <- function(f, x, h = 1e-5)
  vapply(seq_along(x), function(i) {
    e <- numeric(length(x)); e[i] <- h
    (f(x + e) - f(x - e)) / (2 * h)
  }, numeric(1))

num_hess <- function(f, x, h = 1e-4) {
  r <- length(x); H <- matrix(0, r, r)
  for (i in seq_len(r)) for (j in seq_len(r)) {
    ei <- numeric(r); ei[i] <- h; ej <- numeric(r); ej[j] <- h
    H[i, j] <- (f(x+ei+ej) - f(x+ei-ej) - f(x-ei+ej) + f(x-ei-ej)) / (4 * h * h)
  }
  (H + t(H)) / 2
}

## ---- the M2 inversion --------------------------------------------------------
## K_Q(theta) via H-S + Laplace. Returns NA outside the saddlepoint domain
## (where I - 2 theta Hess is not positive-definite).
KQ_hs <- function(theta, Kdy, r) {
  if (theta <= 0) return(0)                       # K_Q(0) = 0
  a <- sqrt(2 * theta)
  g  <- function(w) -0.5 * sum(w^2) + Kdy(a * w)
  opt <- stats::optim(rep(0, r), function(w) -g(w),
                      gr = function(w) w - a * num_grad(Kdy, a * w),
                      method = "BFGS")
  s_star <- a * opt$par
  H  <- num_hess(Kdy, s_star)
  ld <- determinant(diag(r) - 2 * theta * H, logarithm = TRUE)
  if (ld$sign <= 0) return(NA_real_)
  g(opt$par) - 0.5 * as.numeric(ld$modulus)
}

## Lugannani-Rice upper tail of Q from K_Q (numerical theta-derivatives). Scans a
## theta grid up to just below the Gaussian pole to bracket the saddlepoint root.
lr_upper_m2 <- function(Kdy, r, lam, q) {
  m0 <- sum(lam)
  if (q <= m0) return(NA_real_)
  K   <- function(th) KQ_hs(th, Kdy, r)
  ht  <- 1e-4
  Kp  <- function(th) (K(th + ht) - K(th - ht)) / (2 * ht)
  Kpp <- function(th) (K(th + ht) - 2 * K(th) + K(th - ht)) / ht^2
  cap  <- 1 / (2 * max(lam))                       # Gaussian-pole guess
  grid <- seq(cap * 1e-3, cap * 0.999, length.out = 60)
  kp   <- vapply(grid, function(th) tryCatch(Kp(th), error = function(e) NA), numeric(1))
  ok   <- which(is.finite(kp))
  cross <- ok[which(kp[ok] >= q)[1]]               # first theta with K_Q'(th) >= q
  if (is.na(cross) || cross == ok[1]) return(NA_real_)
  br <- c(grid[cross - 1], grid[cross])
  th <- stats::uniroot(function(t) Kp(t) - q, br)$root
  w  <- sign(th) * sqrt(2 * (th * q - K(th)))
  u  <- th * sqrt(Kpp(th))
  min(1, max(0, stats::pnorm(w, lower.tail = FALSE) + stats::dnorm(w) * (1/u - 1/w)))
}

## ---- shared bits: metric rotation, enumeration -------------------------------
metric_sqrt_inv <- function(M) {
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  keep <- e$values > sqrt(.Machine$double.eps) * max(e$values, 1)
  e$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(e$values[keep]), sum(keep), sum(keep))
}
enum_orbit <- function(treatment, strata) {
  idx <- split(seq_along(treatment), as.factor(strata))
  mb  <- vapply(idx, function(ix) sum(treatment[ix]), numeric(1))
  per <- lapply(seq_along(idx), function(b) utils::combn(length(idx[[b]]), mb[b], simplify = FALSE))
  grid <- expand.grid(lapply(per, seq_along))
  Z <- matrix(0, length(treatment), nrow(grid))
  for (k in seq_len(nrow(grid))) {
    z <- numeric(length(treatment))
    for (b in seq_along(idx)) z[idx[[b]][per[[b]][[grid[k, b]]]]] <- 1
    Z[, k] <- z
  }
  Z
}

## build the centred rotated CGF and the pieces M1/M2 need for one (data, metric)
build_case <- function(scores, treatment, strata, M) {
  scores   <- as.matrix(scores)
  rot      <- metric_sqrt_inv(M)
  scores_y <- scores %*% rot
  cg       <- fastperm_linear_cgf_mv(scores_y, treatment, strata)
  mu_y     <- cg$mean
  Ty_obs   <- as.numeric(t(scores_y) %*% as.numeric(treatment))
  lam      <- Re(eigen((cg$cov + t(cg$cov)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  lam      <- lam[lam > sqrt(.Machine$double.eps) * max(lam, 1)]
  ## exact enumerated upper tail of Q = ||d_y||^2
  Z   <- enum_orbit(treatment, strata)
  Dy  <- t(scores_y) %*% Z - mu_y
  Qen <- colSums(Dy^2)
  list(Kdy = function(s) cg$cgf(s) - sum(s * mu_y),
       Sigma_y = cg$cov, r = length(mu_y), lam = lam,
       q_obs = sum((Ty_obs - mu_y)^2), Qen = Qen)
}

K_wchisq  <- function(th, lam) -0.5 * sum(log(1 - 2 * th * lam))
lr_upper_m1 <- function(lam, q) {
  m0 <- sum(lam); if (q <= m0) return(NA_real_)
  Kp <- function(th) sum(lam / (1 - 2 * th * lam))
  Kpp<- function(th) sum(2 * lam^2 / (1 - 2 * th * lam)^2)
  pole <- 1 / (2 * max(lam))
  th <- stats::uniroot(function(t) Kp(t) - q, c(1e-12, pole * (1 - 1e-12)))$root
  w <- sign(th) * sqrt(2 * (th * q - K_wchisq(th, lam))); u <- th * sqrt(Kpp(th))
  stats::pnorm(w, lower.tail = FALSE) + stats::dnorm(w) * (1/u - 1/w)
}
enum_tail <- function(Q, q) mean(Q >= q - 1e-9)

## ============================================================================
## (1) SELF-CHECK: K_Q from H-S+Laplace must equal the Gaussian closed form when
## fed a Gaussian CGF. Use a non-trivial Sigma_y so unequal eigenvalues exercise
## the det term.
## ============================================================================
Sig <- matrix(c(2, 0.6, 0.6, 1), 2, 2)
lamS <- eigen(Sig, only.values = TRUE)$values
Kdy_gauss <- function(s) 0.5 * as.numeric(t(s) %*% Sig %*% s)
cat("=== (1) Gaussian self-check: K_Q(theta) H-S+Laplace vs -1/2 sum log(1-2 theta lambda) ===\n")
cat(sprintf("%8s %14s %14s %12s\n", "theta", "KQ_hs", "KQ_closed", "abs.diff"))
for (th in c(0.02, 0.05, 0.10, 0.15)) {
  hs <- KQ_hs(th, Kdy_gauss, 2)
  cf <- -0.5 * sum(log(1 - 2 * th * lamS))
  cat(sprintf("%8.3f %14.8f %14.8f %12.2e\n", th, hs, cf, abs(hs - cf)))
}

## ============================================================================
## (2) REAL CHECK: a skewed design where M1 (Gaussian-d) should be visibly off.
## n = 12, m = 6, one stratum; col2 highly skewed.
## ============================================================================
sc <- cbind(1:12, c(1, 1, 1, 2, 2, 3, 3, 5, 8, 13, 21, 34))
tr <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
st <- rep(1, 12)

for (mlab in c("Vd", "Sigmax")) {
  M  <- if (mlab == "Vd") build_case(sc, tr, st, diag(2))$Sigma_y else stats::cov(sc)
  ## Vd as a metric: rebuild with M = the permutation covariance on raw scores
  if (mlab == "Vd") {
    cg0 <- fastperm_linear_cgf_mv(sc, tr, st); M <- cg0$cov
  }
  bc <- build_case(sc, tr, st, M)
  cat(sprintf("\n=== (2) skewed n=12 m=6 | metric = %s | r=%d, E[Q]=%.3f, q_obs=%.3f ===\n",
              mlab, bc$r, sum(bc$lam), bc$q_obs))
  cat(sprintf("lambda = %s\n", paste(sprintf("%.3f", bc$lam), collapse = ", ")))
  qs <- sort(unique(c(bc$q_obs,
                      stats::quantile(bc$Qen[bc$Qen > sum(bc$lam)], c(.5, .8), names = FALSE))))
  cat(sprintf("%9s %12s %12s %12s\n", "q", "enum", "M1", "M2"))
  for (q in qs)
    cat(sprintf("%9.4f %12.5f %12.5f %12.5f\n",
                q, enum_tail(bc$Qen, q), lr_upper_m1(bc$lam, q),
                lr_upper_m2(bc$Kdy, bc$r, bc$lam, q)))
}
