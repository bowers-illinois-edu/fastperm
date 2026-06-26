## experiment-bias-trajectory.R --- does M2's tail-probability bias VANISH as B
## grows, and at what rate? The GH/HS pipeline reproduces the exact orbit CGF to
## 1e-14 (experiment-empirical-cgf-check.R), so the cheap empirical-CGF LR (no GH)
## IS M2's tail to 4 decimals. Trace its signed bias to B = 14 to see whether it
## converges to the chi^2_2 floor (+0.5%) or stalls. SCRATCH (dev/).

suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(20260625)

orbit_Q <- function(scoresY, idx, muy) {
  Ty <- scoresY[idx[[1]], , drop = FALSE]
  for (b in 2:length(idx)) {
    opts <- scoresY[idx[[b]], , drop = FALSE]; M <- nrow(Ty); k <- nrow(opts)
    Ty <- Ty[rep(seq_len(M), each = k), , drop = FALSE] +
          opts[rep(seq_len(k), times = M), , drop = FALSE]
  }
  rowSums(sweep(Ty, 2, muy)^2)
}

## empirical-CGF Lugannani-Rice tail for an enumerated orbit Q at threshold q
emp_lr <- function(Q, q) {
  N <- length(Q)
  Kemp <- function(th) { a <- max(th * Q); a + log(sum(exp(th * Q - a))) - log(N) }
  Kp   <- function(th) { a <- max(th * Q); w <- exp(th * Q - a); sum(w * Q) / sum(w) }
  Kpp  <- function(th) { a <- max(th * Q); w <- exp(th * Q - a); w <- w / sum(w)
                         sum(w * Q^2) - (sum(w * Q))^2 }
  hi <- 1e-3; while (Kp(hi) < q && hi < 1e4) hi <- hi * 2
  th <- stats::uniroot(function(t) Kp(t) - q, c(1e-9, hi))$root
  w  <- sign(th) * sqrt(2 * (th * q - Kemp(th))); u <- th * sqrt(Kpp(th))
  stats::pnorm(w, lower.tail = FALSE) + stats::dnorm(w) * (1 / u - 1 / w)
}

one_design <- function(B, ptarget = 0.05) {
  sc <- matrix(stats::rexp(3 * B * 2), 3 * B, 2)
  st <- rep(seq_len(B), each = 3); tr <- rep(c(1, 0, 0), B)
  cg  <- fastperm_linear_cgf_mv(sc, tr, st)
  rot <- .quad_inv_sqrt(cg$cov); scoresY <- sc %*% rot
  muy <- as.numeric(crossprod(rot, cg$mean))
  idx <- split(seq_along(tr), as.factor(st))
  Q   <- orbit_Q(scoresY, idx, muy); N <- length(Q)
  q   <- sort(Q, decreasing = TRUE)[max(1L, round(ptarget * N))]
  exact <- mean(Q >= q - 1e-9)
  m1  <- .quad_lr_upper(c(1, 1), q)            # chi^2_2 LR (lambda = (1,1))
  e2  <- emp_lr(Q, q)                          # = M2 tail to 4 decimals
  c(m1 = (m1 - exact) / exact, m2 = (e2 - exact) / exact)
}

plan <- list(`8` = 120, `9` = 120, `10` = 80, `11` = 50, `12` = 25,
             `13` = 12, `14` = 6)
res <- data.frame()
for (Bs in names(plan)) {
  B <- as.integer(Bs); nd <- plan[[Bs]]
  E <- vapply(seq_len(nd), function(i) one_design(B), numeric(2))
  res <- rbind(res, data.frame(B = B, ndesign = nd, N = 3^B,
    bias_m1 = mean(E["m1", ]), bias_m2 = mean(E["m2", ]),
    se_m2 = stats::sd(E["m2", ]) / sqrt(nd)))
}
print(res, digits = 4)

cat("\nlog-log slope fits (bias magnitude vs B):\n")
for (col in c("bias_m1", "bias_m2")) {
  fit <- stats::lm(log(abs(res[[col]])) ~ log(res$B))
  cat(sprintf("  %-8s  slope = %+.3f  (over B = 8..14)\n", col,
              stats::coef(fit)[2]))
}
cat("\nchi^2_2 LR floor at q = qchisq(0.95, 2):",
    sprintf("%+.4f\n", (.quad_lr_upper(c(1, 1), stats::qchisq(0.95, 2)) -
                        0.05) / 0.05))
