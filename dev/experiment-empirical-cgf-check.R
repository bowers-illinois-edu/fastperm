## experiment-empirical-cgf-check.R --- is M2's +6% large-B tail bias in the CGF
## PIPELINE (Hubbard-Stratonovich + GH-24 + finite-difference derivatives) or in
## the LR inversion of the discrete orbit Q itself?
##
## Reference engine: the EXACT empirical CGF of the enumerated orbit,
## K_emp(theta) = log mean_i exp(theta Q_i), with EXACT tilted moments
## K_emp' = E_tilt[Q], K_emp'' = Var_tilt[Q] (no GH, no HS, no finite differences),
## inverted by the SAME Lugannani-Rice formula as .quad_spa_upper. Comparison:
##   exact     : k/N, the truth
##   M2        : .quad_spa_upper on the GH/HS CGF (the shipped engine)
##   emp-LR    : LR on the exact empirical CGF (isolates CGF pipeline vs inversion)
## Also prints K_Q^GH(theta_hat) vs K_emp(theta_hat) at the saddle. SCRATCH (dev/).

suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(20260625)

prep <- function(scores, treatment, strata, nodes = 24L) {
  cg  <- fastperm_linear_cgf_mv(scores, treatment, strata)
  rot <- .quad_inv_sqrt(cg$cov); scoresY <- scores %*% rot
  muy <- as.numeric(crossprod(rot, cg$mean))
  Sy  <- crossprod(rot, cg$cov %*% rot)
  lam <- Re(eigen((Sy + t(Sy)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  lam <- lam[lam > sqrt(.Machine$double.eps) * max(lam, 1)]
  idx <- split(seq_along(treatment), as.factor(strata))
  mb  <- vapply(idx, function(ix) sum(treatment[ix]), numeric(1)); nb <- lengths(idx)
  KQ  <- .make_KQ_quad(scoresY, idx, mb, nb, muy, .tensor_normal(nodes, length(lam)))
  list(scoresY = scoresY, idx = idx, muy = muy, lam = lam, KQ = KQ, m0 = sum(lam))
}
orbit_Q <- function(scoresY, idx, muy) {
  Ty <- scoresY[idx[[1]], , drop = FALSE]
  for (b in 2:length(idx)) {
    opts <- scoresY[idx[[b]], , drop = FALSE]; M <- nrow(Ty); k <- nrow(opts)
    Ty <- Ty[rep(seq_len(M), each = k), , drop = FALSE] +
          opts[rep(seq_len(k), times = M), , drop = FALSE]
  }
  rowSums(sweep(Ty, 2, muy)^2)
}

## exact empirical-CGF saddlepoint tail for an enumerated orbit Q at threshold q
emp_lr <- function(Q, q) {
  N <- length(Q)
  Kemp <- function(th) { a <- max(th * Q); a + log(sum(exp(th * Q - a))) - log(N) }
  ## exact tilted mean and variance (= K', K'')
  Kp  <- function(th) { a <- max(th * Q); w <- exp(th * Q - a); sum(w * Q) / sum(w) }
  Kpp <- function(th) { a <- max(th * Q); w <- exp(th * Q - a); w <- w / sum(w)
                        m <- sum(w * Q); sum(w * Q^2) - m^2 }
  hi <- 1e-3; while (Kp(hi) < q && hi < 1e4) hi <- hi * 2
  th <- stats::uniroot(function(t) Kp(t) - q, c(1e-9, hi))$root
  w  <- sign(th) * sqrt(2 * (th * q - Kemp(th))); u <- th * sqrt(Kpp(th))
  lr <- stats::pnorm(w, lower.tail = FALSE) + stats::dnorm(w) * (1 / u - 1 / w)
  c(lr = lr, th = th, Kemp = Kemp(th))
}

ndesign <- 40L; ptarget <- 0.05
for (B in c(8L, 10L, 12L)) {
  agg <- matrix(NA_real_, ndesign, 4,
                dimnames = list(NULL, c("m2", "emp", "dKrel", "exact")))
  for (i in seq_len(ndesign)) {
    sc <- matrix(stats::rexp(3 * B * 2), 3 * B, 2)
    st <- rep(seq_len(B), each = 3); tr <- rep(c(1, 0, 0), B)
    p  <- prep(sc, tr, st)
    Q  <- orbit_Q(p$scoresY, p$idx, p$muy); N <- length(Q)
    q  <- sort(Q, decreasing = TRUE)[max(1L, round(ptarget * N))]
    exact <- mean(Q >= q - 1e-9)
    m2 <- suppressWarnings(unname(.quad_spa_upper(p$KQ, q, p$m0, p$lam)["lr"]))
    e  <- emp_lr(Q, q)
    agg[i, ] <- c((m2 - exact) / exact, (e["lr"] - exact) / exact,
                  (p$KQ(e["th"]) - e["Kemp"]) / abs(e["Kemp"]), exact)
  }
  cat(sprintf("\nB=%2d (N=%d, %d designs): mean signed rel error\n", B, 3^B, ndesign))
  cat(sprintf("  M2 (GH/HS pipeline)        : %+.4f\n", mean(agg[, "m2"])))
  cat(sprintf("  emp-LR (exact orbit CGF)   : %+.4f\n", mean(agg[, "emp"])))
  cat(sprintf("  K_Q^GH vs K_emp at saddle  : %+.2e (rel)\n", mean(agg[, "dKrel"])))
}
