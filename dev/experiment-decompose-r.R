## experiment-decompose-r.R --- split M2's r-dependent bias into the INTRINSIC
## continuous chi^2_r saddlepoint error vs the FINITE-set/discreteness error.
## For r = 1, 2, 3 at p = 0.05:
##   (A) analytic continuous chi^2_r LR error: .quad_lr_upper(rep(1,r), qchisq) vs .05
##   (B) finite-set error: emp-LR on n iid EXACT chi^2_r draws (n = 59049 ~ 3^10),
##       a CONTINUOUS law represented by a finite set -- pure finite-set effect
##   (C) orbit error (from experiment-curvature-vs-r.R): +0.155, +0.067, +0.017
## If (A) carries the r-trend, the residual is the chi-square shape (a continuous
## saddlepoint property M1 shares); if (B) does, it is discreteness. SCRATCH (dev/).

suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(20260625)

emp_lr <- function(Q, q) {
  N <- length(Q)
  K   <- function(th) { a <- max(th * Q); a + log(sum(exp(th * Q - a))) - log(N) }
  Kp  <- function(th) { a <- max(th * Q); w <- exp(th * Q - a); sum(w * Q) / sum(w) }
  Kpp <- function(th) { a <- max(th * Q); w <- exp(th * Q - a); w <- w / sum(w)
                        sum(w * Q^2) - (sum(w * Q))^2 }
  hi <- 1e-3; while (Kp(hi) < q && hi < 1e4) hi <- hi * 2
  th <- stats::uniroot(function(t) Kp(t) - q, c(1e-9, hi))$root
  w  <- sign(th) * sqrt(2 * (th * q - K(th))); u <- th * sqrt(Kpp(th))
  stats::pnorm(w, lower.tail = FALSE) + stats::dnorm(w) * (1 / u - 1 / w)
}

orbit_bias <- c(0.155, 0.067, 0.017)   # from experiment-curvature-vs-r.R, B=10
n <- 59049L; reps <- 40L
cat("Decomposition of M2's r-dependent bias at p = 0.05:\n\n")
cat("  r   (A) analytic chi^2_r LR   (B) finite chi^2_r set   (C) orbit (B=10)\n")
for (r in 1:3) {
  qr <- stats::qchisq(0.95, df = r)
  A  <- (.quad_lr_upper(rep(1, r), qr) - 0.05) / 0.05
  bB <- vapply(seq_len(reps), function(i) {
    Q <- stats::rchisq(n, df = r)
    q <- sort(Q, decreasing = TRUE)[round(0.05 * n)]
    (emp_lr(Q, q) - mean(Q >= q - 1e-12)) / mean(Q >= q - 1e-12)
  }, numeric(1))
  cat(sprintf("  %d        %+.4f               %+.4f (se %.3f)        %+.3f\n",
              r, A, mean(bB), stats::sd(bB) / sqrt(reps), orbit_bias[r]))
}
