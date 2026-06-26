## spike-quadratic-spa.R --- prototype + simulation check for the Route B
## quadratic-form saddlepoint. SCRATCH (dev/): validate the math empirically
## before any package code or vignette prose. Jake's rule: the simulation shows
## when the math is confused.
##
## What this checks, in order:
##  (A) the rotation reduction   Q = (T-mu)' M^{-1} (T-mu) = ||d_y||^2   exactly,
##  (B) M1, the Gaussian-d weighted-chi-square tail of Q via Lugannani-Rice,
##      against (i) exact enumeration of the permutation null and (ii) imhof,
##  for two metrics M (M = V_d, the permutation covariance, and M = cov(scores),
##  the RItools Sigma_x default) and two designs (one stratum; two strata).
##
## Run:
##   Rscript -e 'devtools::load_all("~/repos/fastperm-route-b"); \
##               source("~/repos/fastperm-route-b/dev/spike-quadratic-spa.R")'

## ---- metric^{-1/2} via the spectral decomposition (pseudoinverse sqrt) -------
## rot rot' = M^{-1}, dropping numerically-zero eigenvalues so a rank-deficient
## metric still yields a valid rotation onto its range (RItools does the same).
metric_sqrt_inv <- function(M) {
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  keep <- e$values > sqrt(.Machine$double.eps) * max(e$values, 1)
  e$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(e$values[keep]),
                                            sum(keep), sum(keep))
}

## ---- enumerate the within-stratum permutation orbit --------------------------
## Every treatment vector with the same per-stratum treated counts as `treatment`,
## as columns of an n x N matrix. N = prod_b choose(n_b, m_b).
enumerate_orbit <- function(treatment, strata) {
  idx <- split(seq_along(treatment), as.factor(strata))
  mb  <- vapply(idx, function(ix) sum(treatment[ix]), numeric(1))
  per <- lapply(seq_along(idx), function(b)
    utils::combn(length(idx[[b]]), mb[b], simplify = FALSE))   # treated positions
  grid <- expand.grid(lapply(per, seq_along))
  Z <- matrix(0, length(treatment), nrow(grid))
  for (k in seq_len(nrow(grid))) {
    z <- numeric(length(treatment))
    for (b in seq_along(idx)) {
      ix <- idx[[b]]
      z[ix[per[[b]][[grid[k, b]]]]] <- 1
    }
    Z[, k] <- z
  }
  Z
}

## exact permutation distribution of Q = (T-mu)' A (T-mu), A = M^{-1}, uniform
## over the orbit. Returns the Q values (one per orbit member) and the rotated
## scores used, so the caller can check the rotation reduction.
enumerate_Q <- function(scores, treatment, strata, M) {
  scores <- as.matrix(scores)
  Z  <- enumerate_orbit(treatment, strata)
  Tm <- t(scores) %*% Z                      # p x N: T for each orbit member
  mu <- rowMeans(Tm)                         # exact E[T] over the orbit
  D  <- Tm - mu                              # centered, p x N
  A  <- chol2inv(chol((M + t(M)) / 2))       # M^{-1}
  Q_direct <- colSums(D * (A %*% D))         # d' A d per column
  rot <- metric_sqrt_inv(M)
  Dy  <- t(rot) %*% D                        # rotated centered vector
  Q_rot <- colSums(Dy^2)                     # ||d_y||^2
  list(Q = Q_direct, Q_rot = Q_rot, mu = mu, D = D, Sigma = (D %*% t(D)) / ncol(D))
}

## upper-tail share P(Q >= q) from the enumerated distribution (each orbit member
## equally likely); "observed-included" convention is irrelevant here since we
## evaluate the law itself, not a one-sample p-value.
enum_tail <- function(Qvals, q) mean(Qvals >= q - 1e-9)

## ---- M1: Gaussian-d weighted-chi-square tail via Lugannani-Rice --------------
## If d_y ~ N(0, Sigma_y), then Q = ||d_y||^2 = sum_k lambda_k chi^2_1 with
## lambda = eigen(Sigma_y). Its CGF is K(theta) = -1/2 sum log(1 - 2 theta lambda),
## with closed-form K' and K''. We invert the upper tail with Lugannani-Rice.
K_wchisq  <- function(theta, lam) -0.5 * sum(log(1 - 2 * theta * lam))
Kp_wchisq <- function(theta, lam) sum(lam / (1 - 2 * theta * lam))
Kpp_wchisq<- function(theta, lam) sum(2 * lam^2 / (1 - 2 * theta * lam)^2)

lr_upper_wchisq <- function(q, lam) {
  m0 <- sum(lam)                                   # E[Q]
  if (q <= m0) return(NA_real_)                    # prototype: upper tail only
  th_max <- 1 / (2 * max(lam))
  root <- stats::uniroot(function(th) Kp_wchisq(th, lam) - q,
                         c(1e-10, th_max * (1 - 1e-8)))$root
  w <- sign(root) * sqrt(2 * (root * q - K_wchisq(root, lam)))
  u <- root * sqrt(Kpp_wchisq(root, lam))
  1 - stats::pnorm(w) + stats::dnorm(w) * (1 / u - 1 / w)
}

imhof_tail <- function(q, lam) {
  if (!requireNamespace("CompQuadForm", quietly = TRUE)) return(NA_real_)
  CompQuadForm::imhof(q, lambda = lam)$Qq
}

## ---- run one (design, metric) case and print the comparison ------------------
run_case <- function(label, scores, treatment, strata, metric = c("Vd", "Sigmax")) {
  metric <- match.arg(metric)
  scores <- as.matrix(scores)
  ## the metric: V_d is the permutation covariance (from the orbit); Sigma_x is
  ## the within-stratum-pooled sample covariance, here the unstratified cov().
  enum_for_cov <- enumerate_Q(scores, treatment, strata,
                              M = diag(ncol(scores)))            # cheap: get Sigma
  Vd <- enum_for_cov$Sigma
  M  <- if (metric == "Vd") Vd else stats::cov(scores)
  en <- enumerate_Q(scores, treatment, strata, M = M)

  ## (A) rotation reduction
  rot_err <- max(abs(en$Q - en$Q_rot))

  ## rotated covariance and its eigenvalues drive M1
  rot     <- metric_sqrt_inv(M)
  Sigma_y <- t(rot) %*% Vd %*% rot
  lam     <- Re(eigen((Sigma_y + t(Sigma_y)) / 2, symmetric = TRUE,
                       only.values = TRUE)$values)
  lam     <- lam[lam > sqrt(.Machine$double.eps) * max(lam, 1)]

  ## observed Q at the observed treatment vector
  Tobs <- as.numeric(t(scores) %*% as.numeric(treatment))
  qobs <- as.numeric(t(Tobs - en$mu) %*% chol2inv(chol((M + t(M)) / 2)) %*% (Tobs - en$mu))

  cat(sprintf("\n=== %s | metric = %s ===\n", label, metric))
  cat(sprintf("orbit size N = %d, rank r = %d, E[Q] = %.4f, q_obs = %.4f\n",
              length(en$Q), length(lam), sum(lam), qobs))
  cat(sprintf("(A) rotation reduction max|Q_direct - Q_rot| = %.2e\n", rot_err))
  cat(sprintf("lambda = %s\n", paste(sprintf("%.4f", lam), collapse = ", ")))

  ## (B) compare tails over a grid of q above the mean, plus q_obs
  qs <- sort(unique(c(qobs, stats::quantile(en$Q[en$Q > sum(lam)],
                                            c(.5, .75, .9), names = FALSE))))
  cat(sprintf("%8s %12s %12s %12s\n", "q", "enum", "M1(LR)", "imhof"))
  for (q in qs)
    cat(sprintf("%8.4f %12.5f %12.5f %12.5f\n",
                q, enum_tail(en$Q, q), lr_upper_wchisq(q, lam), imhof_tail(q, lam)))
  invisible(NULL)
}

## ---- two small designs -------------------------------------------------------
## col1 a smooth rank-like score; col2 deliberately skewed so d is non-Gaussian
## and M1 (Gaussian-d) should visibly differ from enumeration -- the gap M2 must
## later close.
sc1 <- cbind(c(1, 2, 3, 4, 5, 6),
             c(1, 1, 2, 2, 5, 9))
tr1 <- c(0, 0, 0, 1, 1, 1)                 # treated = the three high units
st1 <- rep(1, 6)

sc2 <- cbind(c(1, 2, 3, 4, 2, 4, 6, 8),
             c(1, 1, 3, 5, 1, 2, 2, 7))
tr2 <- c(0, 0, 1, 1, 0, 0, 1, 1)
st2 <- c(1, 1, 1, 1, 2, 2, 2, 2)

run_case("1 stratum, n=6 m=3", sc1, tr1, st1, "Vd")
run_case("1 stratum, n=6 m=3", sc1, tr1, st1, "Sigmax")
run_case("2 strata, n=8",      sc2, tr2, st2, "Vd")
run_case("2 strata, n=8",      sc2, tr2, st2, "Sigmax")
