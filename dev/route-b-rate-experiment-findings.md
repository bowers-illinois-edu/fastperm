# Route B step (1): the relative-error-vs-B experiment --- findings

Date: 2026-06-25. Scripts: `dev/experiment-rate-vs-B.R`,
`dev/experiment-rate-vs-B-diagnostic.R`, `dev/experiment-empirical-cgf-check.R`,
`dev/experiment-bias-trajectory.R`, `dev/experiment-finite-set-vs-sum.R`.

## The claim under test

The design note framed M2 (the saddlepoint on the exact permutation CGF `K_Q`) as
ESCAPING the Osipov O(n^{-1/4}) barrier that made Kolassa-Robinson set the
quadratic-form tail aside, with the optimistic reading that, because M2 does the
angular integral exactly (Hubbard-Stratonovich + Gauss-Hermite), its relative
error against the exact permutation tail should fall like O(1/B) (log-log slope
-1) as the number of strata B grows, beating the Gaussian M1 (slope -1/2 for
skewed d).

The experiment: matched TRIPLES (n_b = 3, m = 1, so the centred statistic d is
skewed), metric M = V_d (rotated covariance = I, lambda = (1,1), so M1 is exactly
the chi^2_2 limit and its gap to exact is pure non-normality), exact tail by full
3^B orbit enumeration, evaluated at fixed tail depth p ~ 0.05.

## What the simulations actually show

The O(1/B) reading is WRONG. M2's relative error at fixed tail depth does NOT fall
like 1/B; it sits at roughly 3-7% over B = 8..12 and, if anything, grows toward a
finite-orbit floor. The cause is a discreteness effect, isolated below. This is
the kind of correction the user's standing guidance anticipated: "simulations will
show when the math is confused or wrong."

### 1. Single designs are too noisy; average over designs

The first run used one nested random design per B and gave flat/positive slopes
(rel_m1 +0.06, rel_lr +0.09): each B had a different design, so the specific
non-normality and quantile varied. Averaging the relative error over 150 random
designs per B (the rate is a statement about the EXPECTED error) gave clean,
monotone curves. Single realizations cannot show the rate.

### 2. The averaged curves are not power laws

Mean |rel error|, 150 designs, B = 4..10:

```
 B  rel_m1  rel_lr
 4  0.668   0.212
 6  0.466   0.082
 8  0.367   0.050
10  0.289   0.059   <- M2 has stopped falling and ticked UP
```

M2 drops fast to ~0.05 by B = 8 and then plateaus; M1 keeps decreasing. Fitting a
single slope through the plateauing M2 curve returns a spuriously steep -1.49 --- a
fit artifact, not a rate. Two confounds sit at the two ends: small-B orbit
granularity (the 3^B orbit's tail moves in steps of ~1/(0.05 * 3^B): 25% relative
at B = 4, < 1% by B = 7) inflates the steep early drop, and a floor flattens the
tail.

### 3. Decompose: signed bias vs spread; measure the chi^2_2 floor

Averaging |error| conflates the systematic Edgeworth bias (what the rate is about)
with design-to-design spread. Splitting them (200 designs, B = 4..12):

- The chi^2_2 LR floor c_0 (Lugannani-Rice for a TRUE chi^2_2, the B -> infinity
  limit of Q) is TINY: +0.5% at p = 0.05 (+0.2% to +0.8% across the tail). So the
  ~7% plateau is NOT the chi^2_2-LR floor.
- M2's SIGNED bias crosses zero near B = 6 (captured skewness exactly cancels LR
  truncation there) and then climbs: -0.10, -0.03, +0.001, +0.018, +0.047, +0.041,
  +0.053, +0.057, +0.063 for B = 4..12. It does NOT converge to c_0.
- M2's design SPREAD sd_lr falls cleanly (0.22 -> 0.033, slope -1.81). So the floor
  is a systematic positive BIAS, not realization scatter.

### 4. The pipeline is exact; the bias is genuine LR truncation

Comparing M2 (GH/HS CGF + finite-difference derivatives) against an independent
reference --- the EXACT empirical CGF of the enumerated orbit,
K_emp(theta) = log mean_i exp(theta Q_i), with exact tilted moments and the same LR
formula (no GH, no HS, no FD):

```
B    M2 (GH/HS)   emp-LR (exact orbit CGF)   K_Q^GH vs K_emp at saddle
 8   +0.0316      +0.0316                     -3e-14 (rel)
10   +0.0698      +0.0698                     -5e-14
12   +0.0682      +0.0682                     -1e-13
```

M2 and emp-LR agree to four decimals; the GH/HS CGF matches the exact orbit CGF to
1e-14 at the saddle. The +7% is NOT a quadrature / finite-difference / pipeline
artifact. It is the genuine Lugannani-Rice error inverting the exact CGF of the
discrete orbit.

### 5. Trace the bias to B = 14: it grows, it does not vanish

emp-LR == M2 to 4 decimals, and is cheap (no GH), so it traces the bias to larger
B. Signed bias, B = 8..14:

```
 B   bias_m1   bias_m2
 8   +0.340    +0.0345
10   +0.280    +0.0564
12   +0.235    +0.0626
14   +0.198    +0.0806 (noisy, 6 designs)
```

M1 (the CRUDE chi^2_2 CGF) DECREASES toward the +0.5% floor (slope -0.89). M2 (the
EXACT orbit CGF) INCREASES toward ~6-8%. Extrapolating, M1 would overtake M2 near
B ~ 45. M2's advantage is real but confined to small/moderate B --- the regime the
balance-test application lives in, not the asymptotic-in-B limit.

### 6. The cause: finite-orbit discreteness, not curvature

The orbit is the uniform law on 3^B points --- a FINITE set. Inverting the exact
CGF of a finite set with continuous LR overestimates the tail, because the tilt
needed to reach p = 0.05 is influenced by the finite extreme points (the empirical
CGF is max-dominated in the tail, and Q_max grows with B). Control: the uniform law
on n iid EXACT chi^2_2 draws (shape fixed, not a sum) shows the same kind of bias:

```
matched_B   n         emp-LR bias
 8          6561       +0.161
10          59049      +0.100
12          531441     +0.079
```

The chi^2_2-set bias is large and POSITIVE against a law whose true LR error is
+0.5%, and it falls only slowly with n (~ n^{-0.16}). The orbit (structured
Minkowski-sum grid) is better-behaved than random scatter --- +3.5% vs +16% at
matched size B = 8 --- and rises toward the finite-set curve as B grows; they meet
near +7% at N ~ 5e5. So the residual is the finite-support / discreteness gap,
slowly vanishing, NOT a curvature / surface-area effect. This matches the
pre-registered honest residual "(a) discreteness of Q at small B" and contradicts
the curvature framing.

## Revised, simulation-supported claim

- M2 dramatically improves the CONSTANT over M1 and over the crude D^2 tail: 3-7%
  vs 25-67% p-value error at p = 0.05 over B = 8..12, a 4-10x reduction. This is
  the practical win, and it is the sense in which M2 escapes the useless Osipov-rate
  D^2 tail that Kolassa-Robinson set aside.
- M2 does NOT achieve a better asymptotic RATE for the tail probability at fixed
  depth. Its residual is a finite-orbit discreteness floor that decreases only
  slowly (the chi^2_2-set control falls like ~ n^{-0.16}) and can grow with B in the
  moderate range before meeting that floor.
- The honest competitor framing from the literature note stands: vs Imhof, M2's
  gain is a better CONSTANT at small B, not a better rate.

## 7. Mechanism: what the residual is NOT, and what it is

Four further experiments pinned the mechanism, and ruled out the first three
candidate corrections (`dev/experiment-tilt-mechanism.R`,
`dev/experiment-curvature-vs-r.R`, `dev/experiment-decompose-r.R`).

- NOT max-domination of the empirical CGF. At the 0.05 saddle of a B = 12 orbit the
  tilt's effective size is ~6% of N (~32,000 points); the top 100 points carry ~3%
  of the mass. The tilted law is smooth, not dominated by a few extremes. So a
  correction aimed at the tail-max does not apply.
- NOT generic permutation discreteness. The SCALAR sum saddlepoint on the same
  orbit (first rotated coordinate, one-sided tail) is accurate to < 1% (bias +0.029,
  -0.004, +0.001, -0.002, -0.005 across five designs). Robinson's O(1/B) for scalar
  permutation sums holds; the residual is specific to the quadratic form.
- NOT boundary curvature. If the cause were counting discrete points against the
  curved ellipsoid {Q = q} (a Gauss-circle discrepancy), the bias would GROW with
  the dimension r. It FALLS: +15.5% (r = 1), +6.7% (r = 2), +1.7% (r = 3) at B = 10.
  r = 1 (Q = D1^2, the two-sided fold) is the worst, not the curved high-r case.
- NOT the intrinsic chi^2_r saddlepoint. The analytic continuous chi^2_r LR error is
  +0.5% for r = 1, 2, 3 alike; the method handles every r. So a 1D lattice / Daniels
  span correction (built for flat-boundary lattices) is the wrong tool, and r* (the
  smooth higher-order fix) does not move it either.

What it IS: the finite-discrete-representation gap. Representing the same continuous
chi^2_r law by a FINITE set of n points reintroduces the bias the continuous CGF
does not have -- a uniform law on n = 59049 i.i.d. exact chi^2_r draws shows +19.8%
(r=1), +9.7% (r=2), +8.2% (r=3), against an analytic LR error of +0.5%. The orbit is
a structured Minkowski-sum grid, so it is better than random scatter at every r
(+15.5 vs +19.8 at r=1, +1.7 vs +8.2 at r=3) and improves much faster with r.

## 8. Practical reading and recommendation

- The Mahalanobis omnibus is used at r >= 2 (a single covariate is a scalar test,
  not a quadratic form), and there the bias is small and falling: ~7% at r = 2, 1.7%
  at r = 3, and by extrapolation negligible at the r = 10-50 RItools targets. M2 gets
  MORE accurate as covariates are added -- the opposite of a problem for the many-
  covariate backend. The alarming +15% is the r = 1 case, which is not a use case.
- The bias is CONSERVATIVE: M2 reports p-values slightly too large (overstates the
  tail), so it controls Type I error and costs only a little power.
- No standard continuity correction targets this residual: it is not a lattice
  (no fixed span; local span -> 0 would predict a vanishing bias the data contradict),
  not max-domination, not curvature, not the chi-square shape. A genuine finite-
  discrete / Euler-Maclaurin-type correction is research-level with uncertain payoff.

Recommendation: do NOT implement a continuity correction. It is the wrong tool for
the diagnosed mechanism, and unnecessary for the regime the method serves -- M2 is
already 4-10x better than the shipped M1, conservative, accurate to a few percent at
r >= 2, and improving with r. Document the residual (done here and in the vignette)
and move to wiring M2 into riposte. Keep the coarse-lattice (single/double stratum)
Gaussian fallback that already ships.
