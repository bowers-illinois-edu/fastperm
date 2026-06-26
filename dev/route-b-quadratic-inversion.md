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

**SUPERSEDED 2026-06-25: the Laplace step is a proven dead end -- see the
"Update" section at the end of this file. The Hubbard-Stratonovich
representation is correct, but it must be evaluated by quadrature, not Laplace.
The original reasoning is kept below as the record.**

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

## Update 2026-06-25: M1 shipped, Laplace dead end, the quadrature fix, M3 (sparse grids)

State after a build-and-simulate session. Scratch (all reproducible):
`dev/spike-quadratic-spa.R`, `dev/spike-cgf-derivs.R`, `dev/spike-quadratic-m2.R`,
`dev/spike-quadratic-m2b.R`, `dev/spike-m2-quad.R`, `dev/spike-m2-diagnose.R`.

### M1 is done and committed
`fastperm_spa_quadratic(scores, treatment, strata, metric = "cov" | <matrix>,
method = "gaussian")` ships the Gaussian-d weighted-chi-square tail via
Lugannani-Rice (commit 3482d3f on route-b; R CMD check 0/0/0; enumeration +
imhof oracle in tests). The rotation reduction Q = ||rot' d||^2 with
rot rot' = M^{-1} is exact (verified vs enumeration to 1e-9), so the metric is
preprocessing and the core works on a sum of squares. `method = "saddlepoint"`
(M2) currently errors with a "not yet implemented" message.

### The Laplace step (the superseded section above) is wrong -- proven
Leading-order Laplace of the H-S integral around its mode returns EXACTLY the
Gaussian K_Q, carrying none of the non-normality. Reason: K_d is centered, so
grad K_d(0) = 0, hence w = 0 is always a critical point of the inner objective,
and the second-order expansion there gives -1/2 log det(I - 2 theta Sigma) --
the Gaussian CGF. The non-normality lives only in higher-order Laplace terms.
Confirmed numerically: H-S K_Q equalled the Gaussian K_Q to six digits and both
diverged from the exact (enumerated) K_Q in the tail. The analytic
gradient/Hessian built to clean up the Laplace numerics (`spike-cgf-derivs.R`:
grad K = V'p, Hess K = V'(P - pp')V via tilted inclusion probabilities;
validated to recover the closed-form mean and covariance exactly) therefore do
NOT rescue M2. They remain useful -- clean K_Q' for Lugannani-Rice, and the
cumulant infrastructure RItools wants -- but they are not the fix.

### M2 SHIPPED (method = "saddlepoint"), validated against enumeration
The H-S representation is an exact Gaussian expectation:
    M_Q(theta) = E_{W ~ N(0, I_r)}[ exp(K_e(sqrt(2 theta) W)) ],
evaluated by tensor Gauss-Hermite (not Laplace). K_Q(theta) = log M_Q is the EXACT
permutation CGF of Q (Q has bounded support, so K_Q is entire -- no pole, unlike
the Gaussian weighted-chi-square). The inversion is the same Lugannani-Rice as M1
plus the Barndorff-Nielsen r*; only the CGF source changes.

Implementation (R/saddlepoint-quadratic.R): vendored Golub-Welsch GH (.gh_nodes),
tensor N(0,I_r) rule (.tensor_normal), a VECTORISED multi-node elementary-symmetric
evaluation (.log_esym_cols -- the key speedup; the per-row apply prototype took
4 min, the vectorised package call is ~0.1s), the K_Q closure (.make_KQ_quad), and
the upper-tail inverter (.quad_spa_upper: FD theta-derivatives -- GH is machine
precision, so FD is fine -- saddle solve, LR + r*). The public function gains
method = "saddlepoint", returns p.value (LR) and p.value.rstar, and guards
nodes^r > 250000 with an error pointing at the unbuilt M3.

Validation (dev/spike-m2-inversion.R and the tests). GH self-converges to ~1e-15
(GH-40 vs GH-64) for moderate theta, ~1e-6 in the deep tail. On the dense 10-pair
orbit with skewed (fibonacci) scores and M = Sigma_x, M2 tracks the exact tail to
~3-5% where M1 (Gaussian) is off by 15-340%:
    q=25.4: exact 0.0508, M2 0.0524 (3%), M1 0.0884 (74%)
    q=32.4: exact 0.0117, M2 0.0112 (5%), M1 0.0516 (340%)
LR and r* are nearly identical here. Below E[Q], and at the support edge / coarse
lattices where FD goes unreliable, M2 falls back to the Gaussian-d tail (finite,
with a warning) -- never NaN, never a hard error.

KNOWN LIMITATION (next M2 refinement): coarse lattices. For a single or double
stratum, V_d is proportional to cov(scores) (so the metrics coincide) and Q's
support is a coarse, atom-heavy lattice; the continuous saddlepoint then under/
overshoots the discrete tail and needs a continuity correction (lit note sec 5;
Daniels 1987 CC2, 2 sinh(s/2)). M2 is accurate on the DENSE multi-stratum support
that the balance-test regime actually has; the coarse-lattice continuity
correction is the next refinement, not yet implemented.

Port-ready formulas (CGF and derivatives, saddle equation, LR, r*, the Daniels mean
limit) are in `dev/kuonen-spa-formulas.md`, which also records the M1 mean-band fix
(commit 33849f1): the old fallback returned the plain normal tail near the mean,
dropping the skewness term and (over a wider band than its 1e-6 guard) giving
non-monotone values from 1/u - 1/w cancellation; now the Daniels limit plus a
Taylor series of the correction over |w| < 0.1. Tails unchanged.

### M3: do not forget the sparse-grid / QMC work
Tensor Gauss-Hermite costs n^r, feasible only for small r (the number of
representations / covariates): trivial at r = 2-4, infeasible at riposte's ~6,
impossible at RItools' covariate counts (10-50). The H-S expectation is the same
for any r, so the dimension problem is purely the quadrature rule. M3 = a
quadrature scheme matched to r:
 - tensor Gauss-Hermite for small r,
 - sparse-grid (Smolyak) Gauss-Hermite for moderate r,
 - (quasi-)Monte Carlo of the same expectation for large r (scales to any r;
   needs variance control / importance sampling, since exp(K_d) can be heavy).
This is what RItools' many-covariate balance test will actually require, so M3
is not optional polish -- it is the path to the RItools backend at realistic
covariate counts. Two background surveys (academic methods; existing R/C++ for
quadratic-form tails, saddlepoints, and sparse-grid/QMC quadrature) were
launched 2026-06-25 to find better-known methods or reusable code before we
build M3 (and possibly to simplify M2).

### Higher-order inversion (r*) and the Osipov rate question
After reading Kolassa-Robinson (2011) in full (see route-b-literature.md sec 3),
two points bear on the inversion step.

1. Higher-order asymptotics apply to ONE step. The HS identity is exact and the
   GH quadrature error is numerical (node-count-controlled, not asymptotic-in-n),
   so the only asymptotic-in-n error in the pipeline is the final scalar tail
   inversion. Lugannani-Rice there is O(1/n). The scalar Barndorff-Nielsen r*
   (modified signed root w* = w + (1/w) log(v/w), v = theta_hat sqrt(K_Q''))
   upgrades it to O(n^{-3/2}) as a near-free drop-in: same saddle theta_hat, same
   K_Q, only the root changes. K_Q'' comes from the same GH expectation with an
   extra factor, or by finite-differencing K_Q. Implement as
   correction = c("none", "r-star"). The scalar reduction here beats K-R: their
   multivariate r* (eq 9) stays O(1/n); 1-D r* gets the extra half-order off the
   shelf.

2. The Osipov barrier: mechanism RESOLVED, but it points to a different residual
   threat. Full analysis in route-b-literature.md sec 3 (after reading K-R 2011 and
   RHHQ 1990 in full). Short version: the n^{-1/4} is RHHQ's surface-area penalty
   (Cor 2.1: error = (surface * n^{1/2} + volume) * O(...); the boundary term carries
   an extra n^{1/2} for a curved level set, none for a flat one). It is an artifact
   of the Edgeworth-OVER-A-CURVED-REGION method, not a property of Q. Our HS+GH
   integrates over all W directions -- the angular integral over the degenerate
   sphere of boundary points -- EXACTLY, so the surface penalty never arises. Two
   honest caveats: (a) escaping curvature is NOT proving scalar O(1/n) -- the real
   residual threat is the DISCRETENESS of Q's permutation law at small B (ties to
   eigenvalue commensurability, sec 5 of the lit note); (b) the honest competitor is
   Imhof, not Osipov -- as B grows the inner Laplace becomes exact and Imhof is
   itself O(1/B), so our gain over Imhof is a better CONSTANT at small B (carrying
   true skew/kurtosis of d), not a better rate. Decisive experiment has two readings:
   our error vs B (rate; slope -1 = scalar O(1/B), shallower = discreteness not
   Osipov), and our error vs Imhof at fixed small B (constant; favors us, widening
   with skewness). The Table-3 selling point stands either way: K-R's "Quadratic"
   row (D^2 by chi-square) is badly off in the tail and they call accurate D^2 tails
   "not available"; an exact-CGF saddlepoint is a better calibration of the same
   statistic that practitioners already report. (Referee Q "why not Lambda?": we
   keep the deployed Hansen-Bowers d^2 and have the exact CGF K-R lacked.)

### Update 2026-06-26: the decisive experiment ran -- result and mechanism
Full write-up + 8 scripts in dev/route-b-rate-experiment-findings.md and
dev/experiment-*.R; summary in vignette sec 6-7. The "error vs B" reading came back
SHALLOWER than slope -1, so caveat (a) wins over the optimistic O(1/B) reading: M2's
relative error at fixed tail depth p=0.05 does NOT fall like 1/B. It sits at 3-7%
over B=8..12 and grows toward a finite floor, while the crude chi^2 M1 falls toward
its own +0.5% limit (would overtake M2 near B~45). So M2's win is a better CONSTANT
(4-10x smaller p-value error than M1 at realistic B), NOT a better rate -- exactly
the better-constant-vs-Imhof framing of caveat (b).

The residual is the DISCRETENESS of caveat (a), and four controls pin its nature and
kill the obvious fixes. (i) Pipeline exact: GH/HS K_Q matches the empirical orbit
CGF to 1e-14 at the saddle, so M2 == LR-of-exact-orbit-CGF; the gap is genuine LR
error inverting the exact CGF of a FINITE support. (ii) NOT curvature: bias FALLS
with dimension r (+15.5% r=1, +6.7% r=2, +1.7% r=3 at B=10) -- a count of discrete
points against the curved level set {Q=q} (Gauss-circle) would GROW with r. So the
RHHQ surface mechanism, escaped in the smooth/Edgeworth sense, does NOT reappear
through the finite orbit. (iii) NOT the continuous saddlepoint: analytic chi^2_r LR
error is +0.5% for r=1,2,3 alike. (iv) NOT generic permutation discreteness: the
SCALAR sum saddlepoint on the same orbit is accurate <1% (flat boundary fine;
Robinson O(1/B) holds) -- the residual is specific to the quadratic form. What
remains: the gap between continuous LR and the exact tail of a finite-discrete Q;
a uniform law on n iid exact chi^2_r draws shows the same bias, and the orbit
(structured grid) beats random scatter at every r.

Consequences. (1) No standard continuity correction targets this (not a regular
lattice -- local span ~1/N -> 0; not max-domination -- tilt ess ~6% of N; not
curvature; not the chi^2 shape), so the "correction = none/r-star" plan in point 1
stays and we add NO lattice/CC term. r* and LR agree here, so r* does not move it
either. (2) PRACTICAL: the omnibus is used at r>=2 (r=1 is a scalar test), where the
bias is small and FALLS with r -- ~negligible at RItools' r=10-50, i.e. M2 IMPROVES
as covariates are added -- and it is CONSERVATIVE (p slightly too large, controls
Type I). Decision (user-approved 2026-06-26): do NOT build a continuity correction;
wire M2 into riposte instead. A finite-discrete Euler-Maclaurin boundary term is the
only thing that could touch it, is research-level, and is not worth the complexity.
