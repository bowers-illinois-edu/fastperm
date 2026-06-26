## experiment-rate-vs-B.R --- Route B step (1): the relative-error-vs-B experiment.
##
## Claim under test: the M2 saddlepoint tail's relative error against the EXACT
## permutation distribution decreases like O(1/B) (slope -1 on log-log) as the
## number of strata B grows -- i.e. it ESCAPES the Osipov O(n^{-1/4}) barrier
## (slope -1/4) that made Kolassa-Robinson abandon the quadratic form -- and beats
## the Gaussian M1 (expected slope -1/2 for skewed d). The Barndorff-Nielsen r*
## should be steeper still (toward -3/2).
##
## Design: matched TRIPLES (n_b = 3, m = 1) so the centred statistic d is skewed
## (matched pairs give symmetric d and would hide the skewness separation). Metric
## M = V_d, so the rotated covariance is I, lambda = (1, 1), and M1 is exactly the
## chi^2_2 limit -- the gap to exact is pure non-normality. Truth is full
## enumeration of the 3^B orbit (feasible to B = 12). Evaluated at a fixed tail
## depth p ~ 0.05 across B. SCRATCH (dev/).

suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(20260625)

Bmax <- 12
## per-unit scores ~ exponential (skewed), so within-triple centring leaves d
## skewed; r = 2 representations. Nested across B (use the first 3B units).
scores_all <- matrix(stats::rexp(3 * Bmax * 2), 3 * Bmax, 2)
strata_all <- rep(seq_len(Bmax), each = 3)
trt_all    <- rep(c(1, 0, 0), Bmax)              # one treated per triple

## build the rotation, weights, and M2 CGF once per design
prep <- function(scores, treatment, strata, nodes = 32L) {
  cg  <- fastperm_linear_cgf_mv(scores, treatment, strata)
  rot <- .quad_inv_sqrt(cg$cov); scoresY <- scores %*% rot
  muy <- as.numeric(crossprod(rot, cg$mean))
  Sy  <- crossprod(rot, cg$cov %*% rot)          # = I for metric V_d
  lam <- Re(eigen((Sy + t(Sy)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  lam <- lam[lam > sqrt(.Machine$double.eps) * max(lam, 1)]
  idx <- split(seq_along(treatment), as.factor(strata))
  mb  <- vapply(idx, function(ix) sum(treatment[ix]), numeric(1)); nb <- lengths(idx)
  KQ  <- .make_KQ_quad(scoresY, idx, mb, nb, muy, .tensor_normal(nodes, length(lam)))
  list(scoresY = scoresY, idx = idx, muy = muy, lam = lam, KQ = KQ, m0 = sum(lam))
}

## exact orbit Q for an m = 1 design: Cartesian sum of one chosen rotated score
## per stratum (vectorised, no R loop over the 3^B members).
orbit_Q <- function(scoresY, idx, muy) {
  Ty <- scoresY[idx[[1]], , drop = FALSE]
  for (b in 2:length(idx)) {
    opts <- scoresY[idx[[b]], , drop = FALSE]; M <- nrow(Ty); k <- nrow(opts)
    Ty <- Ty[rep(seq_len(M), each = k), , drop = FALSE] +
          opts[rep(seq_len(k), times = M), , drop = FALSE]
  }
  rowSums(sweep(Ty, 2, muy)^2)
}

## one random design at stratum count B: relative error of M1 / M2-LR / M2-r*
## against the exact tail at depth ptarget. mean(Q >= q) is exact (Q is
## real-valued, so no ties / no continuity correction needed).
one_design <- function(B, ptarget = 0.05, nodes = 24L) {
  sc <- matrix(stats::rexp(3 * B * 2), 3 * B, 2)
  st <- rep(seq_len(B), each = 3); tr <- rep(c(1, 0, 0), B)
  p  <- prep(sc, tr, st, nodes)
  Q  <- orbit_Q(p$scoresY, p$idx, p$muy); N <- length(Q)
  q  <- sort(Q, decreasing = TRUE)[max(1L, round(ptarget * N))]
  exact <- mean(Q >= q - 1e-9)
  m1 <- .quad_lr_upper(p$lam, q)
  s2 <- suppressWarnings(.quad_spa_upper(p$KQ, q, p$m0, p$lam))
  c(m1 = abs(m1 - exact) / exact,
    lr = abs(unname(s2["lr"]) - exact) / exact,
    rstar = abs(unname(s2["rstar"]) - exact) / exact)
}

## average the relative error over many random designs at each B -- the asymptotic
## rate is a statement about the EXPECTED error; single designs are too noisy.
Bvals <- 4:10
ndesign <- 150L
res <- data.frame()
for (B in Bvals) {
  errs <- vapply(seq_len(ndesign), function(i) one_design(B), numeric(3))
  res <- rbind(res, data.frame(B = B,
    rel_m1 = mean(errs["m1", ]), rel_lr = mean(errs["lr", ]),
    rel_rstar = mean(errs["rstar", ])))
}
print(res, digits = 4)
cat(sprintf("\nlog-log slopes of MEAN relative error vs B (%d designs each;\n",
            ndesign))
cat("target: M2 ~ -1 (escapes Osipov -1/4); M1 ~ -1/2):\n")
for (col in c("rel_m1", "rel_lr", "rel_rstar")) {
  fit <- stats::lm(log(res[[col]]) ~ log(res$B))
  cat(sprintf("  %-10s slope = %+.3f\n", col, stats::coef(fit)[2]))
}
