## experiment-tilt-mechanism.R --- pick the right continuity correction by looking
## at the mechanism. At the 0.05-tail saddle of a B=12 orbit:
##  (1) is the TILTED law of Q max-dominated by a few extreme orbit points (then a
##      lattice/span correction cannot fix it), or smooth-but-discrete (then it can)?
##  (2) does the SCALAR sum saddlepoint on the same orbit achieve low bias (Robinson
##      O(1/B) for permutation sums) while the quadratic Q does not -- i.e. is the
##      residual specific to the quadratic-form route?
## SCRATCH (dev/).

suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(20260625)

orbit_T <- function(scoresY, idx, muy) {            # full orbit of d = T - mu (rows)
  Ty <- scoresY[idx[[1]], , drop = FALSE]
  for (b in 2:length(idx)) {
    opts <- scoresY[idx[[b]], , drop = FALSE]; M <- nrow(Ty); k <- nrow(opts)
    Ty <- Ty[rep(seq_len(M), each = k), , drop = FALSE] +
          opts[rep(seq_len(k), times = M), , drop = FALSE]
  }
  sweep(Ty, 2, muy)
}

## scalar saddlepoint upper tail from an enumerated 1-D orbit X at threshold x
scalar_lr <- function(X, x) {
  N <- length(X)
  K   <- function(th) { a <- max(th * X); a + log(sum(exp(th * X - a))) - log(N) }
  Kp  <- function(th) { a <- max(th * X); w <- exp(th * X - a); sum(w * X) / sum(w) }
  Kpp <- function(th) { a <- max(th * X); w <- exp(th * X - a); w <- w / sum(w)
                        sum(w * X^2) - (sum(w * X))^2 }
  hi <- 1e-3; while (Kp(hi) < x && hi < 1e4) hi <- hi * 2
  th <- stats::uniroot(function(t) Kp(t) - x, c(1e-9, hi))$root
  w  <- sign(th) * sqrt(2 * (th * x - K(th))); u <- th * sqrt(Kpp(th))
  stats::pnorm(w, lower.tail = FALSE) + stats::dnorm(w) * (1 / u - 1 / w)
}

emp_lr_full <- function(Q, q) {                     # returns tail + saddle diagnostics
  N <- length(Q)
  K   <- function(th) { a <- max(th * Q); a + log(sum(exp(th * Q - a))) - log(N) }
  Kp  <- function(th) { a <- max(th * Q); w <- exp(th * Q - a); sum(w * Q) / sum(w) }
  Kpp <- function(th) { a <- max(th * Q); w <- exp(th * Q - a); w <- w / sum(w)
                        sum(w * Q^2) - (sum(w * Q))^2 }
  hi <- 1e-3; while (Kp(hi) < q && hi < 1e4) hi <- hi * 2
  th <- stats::uniroot(function(t) Kp(t) - q, c(1e-9, hi))$root
  a <- max(th * Q); wt <- exp(th * Q - a); wt <- wt / sum(wt)   # tilted weights
  m <- sum(wt * Q); v <- sum(wt * (Q - m)^2)
  sk <- sum(wt * (Q - m)^3) / v^1.5; ku <- sum(wt * (Q - m)^4) / v^2 - 3
  ess <- 1 / sum(wt^2)                                          # tilt effective size
  topw <- sort(wt, decreasing = TRUE)
  w <- sign(th) * sqrt(2 * (th * q - K(th))); u <- th * sqrt(Kpp(th))
  lr <- stats::pnorm(w, lower.tail = FALSE) + stats::dnorm(w) * (1 / u - 1 / w)
  list(lr = lr, th = th, tilt_skew = sk, tilt_kurt = ku, ess = ess, N = N,
       top1 = topw[1], top10 = sum(topw[1:10]), top100 = sum(topw[1:100]))
}

cat("B=12 orbit, 0.05 tail. Tilted-law shape at the saddle and scalar-vs-quadratic.\n\n")
for (i in 1:5) {
  B <- 12L
  sc <- matrix(stats::rexp(3 * B * 2), 3 * B, 2)
  st <- rep(seq_len(B), each = 3); tr <- rep(c(1, 0, 0), B)
  cg  <- fastperm_linear_cgf_mv(sc, tr, st)
  rot <- .quad_inv_sqrt(cg$cov); scoresY <- sc %*% rot
  muy <- as.numeric(crossprod(rot, cg$mean))
  idx <- split(seq_along(tr), as.factor(st))
  D   <- orbit_T(scoresY, idx, muy)                 # N x 2 centred orbit
  Q   <- rowSums(D^2); N <- length(Q)
  q   <- sort(Q, decreasing = TRUE)[round(0.05 * N)]
  exq <- mean(Q >= q - 1e-9)
  d   <- emp_lr_full(Q, q)
  ## scalar control: first rotated coordinate's orbit, matched 0.05 tail
  X   <- D[, 1]; x <- sort(X, decreasing = TRUE)[round(0.05 * N)]
  exx <- mean(X >= x - 1e-9); sx <- scalar_lr(X, x)
  cat(sprintf("design %d:\n", i))
  cat(sprintf("  Q  : bias=%+.3f  saddle th=%.3f  tilt skew=%.2f kurt=%.2f  ess=%.0f (%.1f%% of N)\n",
              (d$lr - exq) / exq, d$th, d$tilt_skew, d$tilt_kurt, d$ess, 100 * d$ess / N))
  cat(sprintf("       tilt top-1 wt=%.2e  top-10=%.3f  top-100=%.3f\n",
              d$top1, d$top10, d$top100))
  cat(sprintf("  T1 : bias=%+.3f  (scalar sum saddlepoint, same orbit)\n",
              (sx - exx) / exx))
}
