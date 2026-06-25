# Route B: the metric-parameterized quadratic-form saddlepoint, and the RItools backend

Status: design note, work not started. The enabler is built and validated
(`R/cgf-stratified-mv.R`, `fastperm_linear_cgf_mv`): the joint within-stratum
permutation CGF of the representation/covariate vector d, with gradient and
Hessian at 0 reproducing mu and Sigma. The inversion below is the unbuilt
piece. Full prior context: `~/repos/riposte/FASTPERM_INTEGRATION.md` and the
memories `fastperm-route-b-resume`, `fastperm-riposte-api-contract`.

## The object, and why one piece of code serves three front ends

Three different tests reduce to the same computation: the tail of a quadratic
form in a stratified, design-based vector of linear statistics.

    Q = d' M^{-1} d

where d is the vector of within-stratum linear statistics (treated-minus-control
sums of scores or covariates, centered within stratum), distributed under the
within-stratum permutation of the treatment indicator, and M is a fixed
symmetric positive-definite metric. The three call sites differ only in M:

- riposte quadratic combination: M = the pooled S-W covariance Sigma.
- RItools `balanceTest` omnibus d^2 (Hansen-Bowers 2008): M = V_d, the
  randomization covariance of d. Then Q is asymptotically chi^2 with df =
  rank(V_d). (RItools `R/Design.R:993-1052`, p-value via pchisq in
  `R/balanceTest.R:374-376`.)
- RItools `sigma_x` omnibus (devel-sigma-x-omnibus branch): M = Sigma_x, the
  within-stratum-pooled sample covariance of the covariates. Then Q is
  asymptotically a weighted sum of chi-squares, sum_k lambda_k chi^2_1, with
  lambda the eigenvalues of Sigma_x^{-1/2} V_d Sigma_x^{-1/2}. (RItools
  `R/sigma_x_test.R`.)

The point that makes this one function, not three: M is constant across the
permutation. Both V_d and Sigma_x are functions of the covariates and the
per-stratum sizes, not of the treatment vector z (see RItools
`randomization_cov_d` and `default_sigma_x`, both z-free). So M is a parameter
of the inversion, supplied by the caller. Do not hard-wire riposte's Sigma.

## Why a saddlepoint, given what RItools already has

The sigma_x branch is not naive. It already gives Q four reference
distributions (`sigma_x_test.R:264`): a 2-moment Satterthwaite match using
asymptotic moments; the same match using the exact randomization E[Q] and
Var[Q] (the default, `satterthwaite_finite`, from a Finucan fourth-moment port
of Mark Fredrickson's i113-highermoments work); Imhof/Davies inversion of the
weighted chi-square; and direct Monte Carlo of the permutation null.

Each leaves a gap, and the gaps share a shape:

- Imhof/Davies are analytic and exact for the weighted chi-square, but assume d
  is Gaussian. In few-strata, small-n, many-covariate designs -- the regime
  balance testing lives in -- d is skewed and discrete, and that is exactly
  where the Gaussian-d tail drifts.
- `satterthwaite_finite` uses the true first two moments but forces a scaled
  chi-square shape, so it ignores the third and higher cumulants that set a
  tail probability.
- `simulate` is the honest reference, but it is Monte Carlo: O(1/sqrt(B)) noise,
  and resolving a small tail probability needs very large B.

The saddlepoint fills the empty cell: analytic and fast like Imhof, but built
on the exact permutation CGF of d, so it carries the true skewness and kurtosis
into the tail rather than assuming them away. It does not replace `simulate` as
ground truth; it reaches simulate's target without the noise and at a fraction
of the cost for small p-values. For M = V_d (the chi-square case) the gain is
smaller, since the chi-square shape is already asymptotically right; it refines
the small-sample tail. For M = Sigma_x the gain is real, because both the
non-normality of d and the weighted-chi-square shape bite at once.

## The inversion: Hubbard-Stratonovich, then Laplace, then Lugannani-Rice

We have the joint CGF K_d(t) = log E[exp(t' d)] (built). We want
K_Q(theta) = log E[exp(theta * d' A d)] with A = M^{-1}, then a tail from it.

Linearize the quadratic by a Gaussian auxiliary integral (Hubbard-Stratonovich).
For theta > 0 and A positive-definite, with A^{1/2} = M^{-1/2},

    exp(theta d'A d) = (2 pi)^{-p/2} INT exp(-||w||^2/2 + sqrt(2 theta) d' A^{1/2} w) dw.

Take the expectation over d inside the integral; the d-dependence is now linear,
so it becomes K_d evaluated at sqrt(2 theta) A^{1/2} w:

    E[exp(theta d'A d)] = (2 pi)^{-p/2} INT exp(-||w||^2/2 + K_d(sqrt(2 theta) A^{1/2} w)) dw.

Because d has finite support (a finite permutation set), K_d is entire and grows
at most linearly, so the -||w||^2/2 term dominates and this p-dimensional
integral converges for all theta in the relevant range. Approximate it by
Laplace's method: find the maximizer w* of the exponent, expand to second order,
and read off K_Q(theta) (the metric enters only through A^{1/2}, confirming M is
a clean parameter). Then apply Lugannani-Rice to K_Q and its theta-derivatives
for P(Q > q). The single-statistic tail already uses log-stable Lugannani-Rice
(`R/saddlepoint.R`); reuse it.

This is the Kuonen (2003) saddlepoint for quadratic forms, generalized from
normal d to the permutation CGF of d via the H-S representation.

## Validation (same discipline as the rest of the package)

1. Exact enumeration of Q on a small stratified case (full within-stratum
   enumeration of z), the primary oracle -- as for the single-rep tail and the
   enabler.
2. The RItools `simulate` backend at large B as a second oracle for the
   Sigma_x case; the saddlepoint tail should match it within Monte Carlo error
   and resolve further into the tail than B allows.
Check all three metrics: M = I, M = V_d (recover chi^2_rank in the large-stratum
limit), M = Sigma_x (recover the weighted-chi-square tail when d is near-Gaussian).

## API sketch

A single entry that takes the raw per-stratum inputs and the metric, and returns
both the tail and the cumulants, so callers can build cheaper approximations from
the same source:

    fastperm_quadratic_spa(scores, treatment, strata,
                           metric = c("cov", <matrix>),  # M; "cov" = V_d
                           observed = NULL)
    -> list(p_upper, statistic, K_Q, cumulants = c(k1, k2, k3, k4), saddle, rank)

Returning the cumulants matters for the consolidation below: the first two feed
a Satterthwaite match, all four feed higher-order corrections, and the same CGF
feeds the saddlepoint. One validated source replaces three hand-derived ones.

## RItools integration (later; do not start until Route B is built and the branch settles)

The integration is small and non-invasive: add a fifth option,
`null = "saddlepoint"`, to `sigma_x_test()` that calls `fastperm_quadratic_spa`
with M = Sigma_x. Nothing else in RItools changes. The branch currently has live
vignette edits, so stay clear of it until then.

The larger payoff is consolidation, not just a shared engine. RItools is
computing exact randomization cumulants in at least three places --
`sigma_x_T_moments` here, the i113-highermoments branch it ports from, and (for
the first two cumulants) the HB08 covariance. fastperm's stratified CGF yields
all of them from one validated recursion. The endgame is RItools consuming the
fastperm cumulants/CGF for both `satterthwaite_finite` and `saddlepoint`,
retiring the hand-ported Finucan formulas.

Constraint: RItools is on CRAN and fastperm is not. The first integration must
keep fastperm in Suggests with the existing backends as the default fallback
(mirror riposte's permute-default / saddlepoint-optional pattern). A hard
dependency waits until fastperm itself is on CRAN.

## The max combination is a separate, later problem

riposte's T_max needs a joint-tail probability, not a quadratic-form tail --
either a multivariate saddlepoint or pmvnorm on the exact S-W correlation. Out
of scope here; sequence it after the quadratic inversion.
