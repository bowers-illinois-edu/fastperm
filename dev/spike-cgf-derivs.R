## spike-cgf-derivs.R --- prototype + validate the ANALYTIC gradient and Hessian
## of the joint within-stratum permutation CGF. SCRATCH (dev/).
##
## The per-stratum CGF K_b(s) = log E[exp(s' T_b)] for T_b = sum_{i in S} v_i over
## a uniform m-subset S is a log-partition function, so
##   grad K_b(s) = E_tilt[T_b] = V_b' p,         p_i  = Pr_tilt(i in S),
##   Hess K_b(s) = Cov_tilt[T_b] = V_b' (P - p p') V_b,
##                                               P_ii = p_i, P_ij = Pr_tilt(i,j in S).
## The tilt weights are exp(s' v_i). Inclusion probabilities come from the same
## elementary-symmetric polynomials the CGF already uses:
##   p_i  = y_i e_{m-1}(y_{-i})   / e_m(y),
##   p_ij = y_i y_j e_{m-2}(y_{-i,-j}) / e_m(y),   y_i = exp(s' v_i),
## all evaluated in log space (sums of positive terms, no cancellation).
##
## Checks: (1) grad/Hess match numerical differentiation of the existing cgf
## closure at several tilts; (2) at s = 0 they recover the closed-form mean and
## covariance (cg$mean, cg$cov), i.e. p_i = m/n and p_ij = m(m-1)/(n(n-1)).
##
## Run:
##   Rscript -e 'devtools::load_all("~/repos/fastperm-route-b"); \
##               source("~/repos/fastperm-route-b/dev/spike-cgf-derivs.R")'

## log e_m of exp(x), guarding the degenerate degrees the recursion mishandles
.loge_safe <- function(x, m) {
  if (m < 0 || m > length(x)) return(-Inf)   # e_m = 0
  if (m == 0) return(0)                       # e_0 = 1
  .log_esym(x, m)                             # from R/cgf-stratified.R
}

## per-stratum inclusion probabilities at projected weights w_i = s' v_i
.incl_probs <- function(w, m) {
  n <- length(w)
  loge_m <- .loge_safe(w, m)
  p <- vapply(seq_len(n),
              function(i) exp(w[i] + .loge_safe(w[-i], m - 1) - loge_m),
              numeric(1))
  P <- matrix(0, n, n)
  if (n >= 2) for (i in 1:(n - 1)) for (j in (i + 1):n) {
    pij <- exp(w[i] + w[j] + .loge_safe(w[-c(i, j)], m - 2) - loge_m)
    P[i, j] <- pij; P[j, i] <- pij
  }
  diag(P) <- p
  list(p = p, P = P)
}

## analytic grad and Hess of K_T(s) summed over the independent strata
cgf_mv_derivs <- function(scores, treatment, strata, s) {
  scores <- as.matrix(scores); d <- ncol(scores)
  idx <- split(seq_along(treatment), as.factor(strata))
  grad <- numeric(d); hess <- matrix(0, d, d)
  for (ix in idx) {
    Vb <- scores[ix, , drop = FALSE]
    m  <- sum(treatment[ix]); nb <- length(ix)
    if (m == 0) next
    if (m == nb) { grad <- grad + colSums(Vb); next }   # deterministic stratum
    ip   <- .incl_probs(as.numeric(Vb %*% s), m)
    grad <- grad + as.numeric(crossprod(Vb, ip$p))
    hess <- hess + crossprod(Vb, (ip$P - tcrossprod(ip$p)) %*% Vb)
  }
  list(grad = grad, hess = hess)
}

## numerical gradient/Hessian for the cross-check
num_grad <- function(f, x, h = 1e-5)
  vapply(seq_along(x), function(i) { e <- numeric(length(x)); e[i] <- h
    (f(x + e) - f(x - e)) / (2 * h) }, numeric(1))
num_hess <- function(f, x, h = 1e-4) {
  r <- length(x); H <- matrix(0, r, r)
  for (i in seq_len(r)) for (j in seq_len(r)) {
    ei <- numeric(r); ei[i] <- h; ej <- numeric(r); ej[j] <- h
    H[i, j] <- (f(x+ei+ej) - f(x+ei-ej) - f(x-ei+ej) + f(x-ei-ej)) / (4*h*h)
  }
  (H + t(H)) / 2
}

## --- designs: matched pairs and triples (the small-stratum, non-Gaussian
## --- regime where M2 will matter), plus one larger stratum ------------------
mk <- function(scores, treatment, strata, label, svals) {
  cg <- fastperm_linear_cgf_mv(scores, treatment, strata)
  cat(sprintf("\n=== %s ===\n", label))
  ## (1) match numerical differentiation at several tilts
  for (s in svals) {
    an <- cgf_mv_derivs(scores, treatment, strata, s)
    ng <- num_grad(cg$cgf, s); nh <- num_hess(cg$cgf, s)
    cat(sprintf("s=(%s): max|grad-num|=%.2e  max|hess-num|=%.2e\n",
                paste(sprintf("%.2f", s), collapse = ","),
                max(abs(an$grad - ng)), max(abs(an$hess - nh))))
  }
  ## (2) s = 0 recovers the closed-form mean and covariance
  a0 <- cgf_mv_derivs(scores, treatment, strata, c(0, 0))
  cat(sprintf("s=0: max|grad-mean|=%.2e  max|hess-cov|=%.2e\n",
              max(abs(a0$grad - cg$mean)), max(abs(a0$hess - cg$cov))))
}

## 6 matched pairs (n_b = 2, m = 1), skewed second covariate
pr_sc <- cbind(c(1, 3, 2, 5, 1, 4, 2, 6, 3, 8, 1, 9),
               c(1, 2, 1, 8, 2, 3, 1, 13, 2, 21, 1, 34))
pr_tr <- rep(c(0, 1), 6)
pr_st <- rep(1:6, each = 2)
mk(pr_sc, pr_tr, pr_st, "6 matched pairs", list(c(0,0), c(.1,-.05), c(.3,.2)))

## 4 matched triples (n_b = 3, m = 1)
tr_sc <- cbind(c(1,2,3, 2,4,6, 1,5,2, 3,1,8),
               c(1,1,2, 1,2,5, 2,1,1, 1,3,7))
tr_tr <- rep(c(1,0,0), 4)
tr_st <- rep(1:4, each = 3)
mk(tr_sc, tr_tr, tr_st, "4 matched triples", list(c(0,0), c(.2,-.1), c(-.15,.25)))

## one larger stratum n = 8, m = 4
lg_sc <- cbind(1:8, c(1,1,2,3,5,8,13,21))
lg_tr <- c(0,0,0,0,1,1,1,1)
lg_st <- rep(1, 8)
mk(lg_sc, lg_tr, lg_st, "1 stratum n=8 m=4", list(c(0,0), c(.05,.02), c(.1,-.03)))
