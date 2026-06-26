## experiment-rate-vs-B-diagnostic.R --- decompose the rate-vs-B result.
##
## The averaged experiment (experiment-rate-vs-B.R) showed M2's relative error
## DROP from 0.21 (B=4) to ~0.05 (B=8) and then PLATEAU. A single power-law slope
## (-1.49) is misleading: the curve is a decrease to a floor, contaminated at small
## B by orbit granularity. This script decomposes the two confounds.
##
## (a) FLOOR c_0: the Lugannani-Rice error for a TRUE chi^2_2 (= the metric-V_d
##     limit of Q as B -> infinity). Both M1 and M2 must converge to this; it is
##     not removable by a better CGF. Measured directly, no permutation.
## (b) GRANULARITY: the orbit has only 3^B points, so the empirical tail at depth p
##     moves in steps of ~1/(p*3^B) -- 25% relative at B=4, <1% by B=7. This
##     inflates the small-B errors and the fitted slope.
## (c) BIAS vs SPREAD: averaging |rel error| conflates the systematic Edgeworth
##     bias (what the rate is about) with design-to-design spread. Report the
##     SIGNED mean and the SD separately, and fit the rate on B >= 6 only.
## SCRATCH (dev/).

suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(20260625)

## ---- (a) the chi^2_2 LR floor c_0, across the tail depths the experiment hits ----
cat("(a) Lugannani-Rice error for a TRUE chi^2_2 (the B -> infinity floor):\n")
cat("    q      p_true     LR        rel_err\n")
for (q in c(4, 4.6, 5.0, 5.5, 5.991, 6.5, 8)) {
  ptrue <- stats::pchisq(q, df = 2, lower.tail = FALSE)
  lr    <- .quad_lr_upper(c(1, 1), q)
  cat(sprintf("  %5.3f  %.5f   %.5f   %+.4f\n", q, ptrue, lr,
              (lr - ptrue) / ptrue))
}

## ---- (b)+(c) signed bias and spread vs B, granularity quantified ----
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

## SIGNED relative error of M1 and M2-LR for one random design
one_design <- function(B, ptarget = 0.05, nodes = 24L) {
  sc <- matrix(stats::rexp(3 * B * 2), 3 * B, 2)
  st <- rep(seq_len(B), each = 3); tr <- rep(c(1, 0, 0), B)
  p  <- prep(sc, tr, st, nodes)
  Q  <- orbit_Q(p$scoresY, p$idx, p$muy); N <- length(Q)
  q  <- sort(Q, decreasing = TRUE)[max(1L, round(ptarget * N))]
  exact <- mean(Q >= q - 1e-9)
  m1 <- .quad_lr_upper(p$lam, q)
  s2 <- suppressWarnings(.quad_spa_upper(p$KQ, q, p$m0, p$lam))
  c(m1 = (m1 - exact) / exact, lr = (unname(s2["lr"]) - exact) / exact)
}

Bvals <- 4:12; ndesign <- 200L; ptarget <- 0.05
res <- data.frame()
for (B in Bvals) {
  E <- vapply(seq_len(ndesign), function(i) one_design(B, ptarget), numeric(2))
  res <- rbind(res, data.frame(
    B = B, N = 3^B,
    gran = 1 / (ptarget * 3^B),               # relative tail granularity
    bias_m1 = mean(E["m1", ]), sd_m1 = stats::sd(E["m1", ]),
    bias_lr = mean(E["lr", ]), sd_lr = stats::sd(E["lr", ]),
    mae_lr  = mean(abs(E["lr", ]))))
}
cat("\n(b)+(c) signed bias, spread, and granularity vs B (",
    ndesign, "designs):\n", sep = "")
print(res, digits = 3)

cat("\nslopes fit on B >= 6 (granularity < 3%, so saddlepoint-limited):\n")
sub <- res[res$B >= 6, ]
for (col in c("bias_m1", "bias_lr", "mae_lr", "sd_lr")) {
  y <- abs(sub[[col]])
  fit <- stats::lm(log(y) ~ log(sub$B))
  cat(sprintf("  |%-8s| slope = %+.3f\n", col, stats::coef(fit)[2]))
}
