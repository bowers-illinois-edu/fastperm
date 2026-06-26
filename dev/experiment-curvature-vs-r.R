## experiment-curvature-vs-r.R --- test the curved-boundary discreteness mechanism
## and its practical scaling. If M2's residual is the count of discrete orbit points
## against the CURVED level set {Q = q} (a Gauss-circle-type discrepancy), then:
##   r = 1: boundary {D1^2 >= q} is two points -- NO curvature -> bias ~ 0
##   r = 2: circle boundary -> bias ~ +7%
##   r = 3: sphere boundary, more surface -> bias larger still
## This both confirms the mechanism and answers whether M2 degrades as the covariate
## count r grows -- the RItools many-covariate regime. SCRATCH (dev/).

suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(20260625)

orbit_T <- function(scoresY, idx, muy) {
  Ty <- scoresY[idx[[1]], , drop = FALSE]
  for (b in 2:length(idx)) {
    opts <- scoresY[idx[[b]], , drop = FALSE]; M <- nrow(Ty); k <- nrow(opts)
    Ty <- Ty[rep(seq_len(M), each = k), , drop = FALSE] +
          opts[rep(seq_len(k), times = M), , drop = FALSE]
  }
  sweep(Ty, 2, muy)
}
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

one_design <- function(B, r, ptarget = 0.05) {
  sc <- matrix(stats::rexp(3 * B * r), 3 * B, r)
  st <- rep(seq_len(B), each = 3); tr <- rep(c(1, 0, 0), B)
  cg  <- fastperm_linear_cgf_mv(sc, tr, st)
  rot <- .quad_inv_sqrt(cg$cov); scoresY <- sc %*% rot
  muy <- as.numeric(crossprod(rot, cg$mean))
  idx <- split(seq_along(tr), as.factor(st))
  D   <- orbit_T(scoresY, idx, muy); Q <- rowSums(D^2); N <- length(Q)
  q   <- sort(Q, decreasing = TRUE)[max(1L, round(ptarget * N))]
  exact <- mean(Q >= q - 1e-9)
  (emp_lr(Q, q) - exact) / exact
}

B <- 10L; nd <- 80L
cat(sprintf("emp-LR (= M2) signed bias at p=0.05, B=%d, %d designs, by dimension r:\n", B, nd))
cat("(prediction if curved-boundary discreteness: r=1 flat ~0, bias grows with r)\n\n")
cat("  r   mean_bias      se\n")
for (r in 1:3) {
  b <- vapply(seq_len(nd), function(i) one_design(B, r), numeric(1))
  cat(sprintf("  %d   %+.4f     %.4f\n", r, mean(b), stats::sd(b) / sqrt(nd)))
}
