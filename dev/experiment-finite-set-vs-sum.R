## experiment-finite-set-vs-sum.R --- WHY does the exact-CGF (emp-LR) tail bias
## GROW with B for the permutation orbit, instead of falling to the chi^2_2 LR
## floor (+0.5%)?
##
## The orbit is the uniform distribution on 3^B points. Two candidate causes:
##  (i)  generic finite-SET artifact: the empirical CGF of any finite set is tail-
##       dominated by its largest member, biasing the saddlepoint. If so, a finite
##       set of n iid EXACT chi^2_2 draws (a set whose shape is fixed, not a sum)
##       would show the same bias -- and it should SHRINK as n grows (better CGF).
##  (ii) orbit SUM structure: Q = ||sum of B bounded vectors||^2, so the support and
##       Q_max GROW with B. If the bias is specific to this, the chi^2_2 finite-set
##       control will NOT reproduce the growing-with-n bias.
##
## Test A: uniform on n iid chi^2_2 draws, emp-LR at its own 0.05 quantile, vs 0.05.
##         n matched to 3^B for B = 8, 10, 12. Many reps. SCRATCH (dev/).

suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(20260625)

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

## bias of emp-LR for the uniform distribution on n iid chi^2_2 draws,
## measured against the SET's own 0.05 quantile (exactly analogous to the orbit)
chisq_set_bias <- function(n, ptarget = 0.05) {
  Q <- stats::rchisq(n, df = 2)
  q <- sort(Q, decreasing = TRUE)[max(1L, round(ptarget * n))]
  exact <- mean(Q >= q - 1e-12)
  (emp_lr(Q, q) - exact) / exact
}

cat("Test A: emp-LR bias for the uniform law on n iid EXACT chi^2_2 draws\n")
cat("(if the orbit bias is a generic finite-set artifact, this matches it and\n")
cat(" shrinks as n grows; the chi^2_2 LR floor against the TRUE law is +0.5%)\n\n")
cat("matched_B      n   reps   mean_bias       se\n")
for (bn in list(c(8, 6561, 200), c(10, 59049, 80), c(12, 531441, 25))) {
  B <- bn[1]; n <- bn[2]; reps <- bn[3]
  b <- vapply(seq_len(reps), function(i) chisq_set_bias(n), numeric(1))
  cat(sprintf("  B=%2d  %7d   %4d   %+.4f   %.4f\n",
              B, n, reps, mean(b), stats::sd(b) / sqrt(reps)))
}
