# Kuonen (1999) saddlepoint formulas for the Gaussian quadratic-form tail (M1 reference)

Port-ready reference for the Gaussian-d (weighted-chi-square) quadratic-form tail:
the M1 closed form our HS+GH route (M2) generalizes by replacing the closed-form
chi-square-mixture CGF with the exact non-Gaussian K_Q. Formulas relayed by the
LitSurvey teammate 2026-06-25 and verified here by direct differentiation; cite
Kuonen, D. (1999), "Saddlepoint approximations for distributions of quadratic
forms in normal variables," Biometrika 86(4):929-935, DOI 10.1093/biomet/86.4.929,
and Daniels (1987), Int. Stat. Rev. 55(1):37-48 for the limiting form. Companion:
the design note `route-b-quadratic-inversion.md` and `route-b-literature.md`.

## Statistic and reduction

Q = d' A d with d ~ N(mu_d, Sigma_d). Reduce to a weighted sum of independent
noncentral chi-squares:

    Q =_d sum_{i=1}^{m} lambda_i X_i,   X_i ~ chi^2_{h_i}(delta_i)  independent.

Reduction (compute once): B = Sigma_d^{1/2} A Sigma_d^{1/2}, eigen-decompose
B = P diag(.) P'. The DISTINCT eigenvalues are lambda_i with multiplicity h_i (the
chi-square degrees of freedom). Set nu = P' Sigma_d^{-1/2} mu_d; the noncentrality
for eigenvalue lambda_i is delta_i = sum_{j : eigval_j = lambda_i} nu_j^2. Central
case mu_d = 0 gives all delta_i = 0.

## Our application is always CENTRAL (delta_i = 0)

Under the within-stratum permutation null, d = T - mu is centered, so mu_d = 0 and
every delta_i = 0. Our shipped M1 (`R/saddlepoint-quadratic.R`, `.quad_lr_upper`)
lists each eigenvalue separately (h_i = 1 per entry, repeats = multiplicity) and
uses the central CGF below. The noncentral terms are recorded here only for
completeness and any future non-null / power calculation; they are NOT needed for
the null tail. Two metric checks:
 - M = V_d (unadjusted): A = V_d^{-1}, d ~ N(0, V_d), so lambda_i = 1, h_i = 1,
   delta_i = 0 and Q ~ chi^2_r exactly; K(zeta) collapses to -(r/2) log(1-2 zeta).
 - M = Sigma_x: distinct lambda_i = eigenvalues of Sigma_x^{-1/2} V_d Sigma_x^{-1/2}
   (genuine weighted sum of central chi-squares).

## CGF and derivatives

Single noncentral chi-square: K_X(t) = -(h/2) log(1-2t) + delta t/(1-2t). Hence

    K(zeta)    = sum_i [ -(h_i/2) log(1 - 2 zeta lambda_i)
                         + delta_i lambda_i zeta / (1 - 2 zeta lambda_i) ]
    K'(zeta)   = sum_i [  h_i lambda_i / (1 - 2 zeta lambda_i)
                         + delta_i lambda_i / (1 - 2 zeta lambda_i)^2 ]
    K''(zeta)  = sum_i [ 2 h_i lambda_i^2 / (1 - 2 zeta lambda_i)^2
                         + 4 delta_i lambda_i^2 / (1 - 2 zeta lambda_i)^3 ]
    K'''(zeta) = sum_i [ 8 h_i lambda_i^3 / (1 - 2 zeta lambda_i)^3
                         + 24 delta_i lambda_i^3 / (1 - 2 zeta lambda_i)^4 ]

Cumulants at zeta = 0 (verified against the standard weighted-noncentral-chi-square
moments): K'(0) = sum lambda_i (h_i + delta_i) = E[Q];
K''(0) = 2 sum lambda_i^2 (h_i + 2 delta_i);
K'''(0) = 8 sum lambda_i^3 (h_i + 3 delta_i).

## Domain and saddlepoint equation

All lambda_i > 0 (A positive definite): K is analytic and steep on
zeta in (-inf, 1/(2 lambda_max)), lambda_max = max_i lambda_i. K'(zeta) rises from
0+ to +inf, so any observed x in (0, inf) has a unique saddle zeta-hat solving
K'(zeta-hat) = x, with zeta-hat > 0 iff x > E[Q]. Solve by Newton (derivative K'')
or bisection on (lower, 1/(2 lambda_max)); clamp strictly below the pole, e.g.
(1 - 1e-10)/(2 lambda_max), to keep all s_i = 1 - 2 zeta lambda_i > 0. (Indefinite
A would use (1/(2 lambda_min), 1/(2 lambda_max)); not needed for the PD metric.)

## Lugannani-Rice upper tail (x != E[Q])

    w-hat = sign(zeta-hat) * sqrt( 2 ( zeta-hat * x - K(zeta-hat) ) )   (sqrt >= 0)
    u-hat = zeta-hat * sqrt( K''(zeta-hat) )
    P(Q >= x) = 1 - Phi(w-hat) + phi(w-hat) * ( 1/u-hat - 1/w-hat )

Phi, phi the standard normal CDF and pdf; the bracket is (1/u-hat - 1/w-hat). Both
terms diverge as x -> E[Q] but the difference has a finite limit (next section).

## Removable singularity at x = E[Q] (Daniels limiting form)

    P(Q >= E[Q]) = 1/2 - K'''(0) / ( 6 * sqrt(2 pi) * K''(0)^{3/2} )

Switch to this when |zeta-hat| (or |x - E[Q]|) is below a small tolerance, to avoid
0/0. For a weighted sum of chi-squares K'''(0) > 0, so the value at the mean is
BELOW 1/2.

### KNOWN M1 GAP (found 2026-06-25): the mean-band fallback drops the skewness term
Our shipped `.quad_lr_upper` falls back to the plain normal tail
`pnorm((q - m0)/sd0)` when |th| < 1e-6 (and at q = m0 returns exactly 1/2). That
DISCARDS the -K'''(0)/(6 sqrt(2 pi) K''(0)^{3/2}) skewness correction above. The
error is confined to a narrow band around p ~ 0.5 and is irrelevant for tail
p-values (where the full LR formula runs), so M1's tail behaviour and all tests
are unaffected. The cheap fix when M1 is next touched: replace the |th|-small
fallback with the Daniels limiting form (exact at the mean) rather than the normal
tail. Track with the M2 build.

## r* (Barndorff-Nielsen) higher-order variant

Same w-hat, u-hat:

    r-star = w-hat + (1/w-hat) * log( u-hat / w-hat )
    P(Q >= x) = 1 - Phi(r-star)

Relative error O(n^{-3/2}) vs Lugannani-Rice's O(n^{-1}); a drop-in once w-hat and
u-hat are computed. (Same root as the M2 r* option; see the design note.)

## Numerical notes for the port

 - Accumulate s_i = 1 - 2 zeta lambda_i in one pass and reuse across K, K', K''.
 - Clamp zeta-hat strictly below 1/(2 lambda_max).
 - The HS+GH M2 engine reuses this EXACT Lugannani-Rice / r* inversion; only the
   CGF source changes (closed-form chi-square mixture here -> exact non-Gaussian
   K_Q = log E_W[exp(K_d(sqrt(2 theta) W))] there). That swap is the M1 -> M2 step.
