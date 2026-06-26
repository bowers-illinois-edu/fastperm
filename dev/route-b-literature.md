# Route B: literature and positioning (survey synthesis, 2026-06-25)

Synthesis of five background literature/code surveys run 2026-06-25, on the tail
of a quadratic form Q = (T - mu)' M^{-1} (T - mu) in a NON-Gaussian (lattice
permutation) vector, by saddlepoint. References carry a verification status; the
flagged ones must be confirmed directly before they enter a manuscript (per the
bibliography rules in `~/repos/ai_workflow/CLAUDE_BIB.md`). Companion: the design
note `route-b-quadratic-inversion.md` (the method + the M1/M2/M3 plan).

## 1. The gap: every standard quadratic-form tail method assumes Gaussian d

Imhof (1961), Davies (1980), Farebrother (1984, Ruben's method), Liu-Tang-Zhang
(2009), and the saddlepoint version Kuonen (1999) all reduce d'Ad to a weighted
sum of chi-squares via the spectral decomposition of A. That reduction is a
Gaussian-only fact (squared standard normals are chi-squares; the diagonalizing
rotation preserves independence only for normals). `CompQuadForm`'s title is
literally "Distribution Function of Quadratic Forms in Normal Variables." So for
our non-Gaussian lattice d there is NO off-the-shelf method -- the four oracles
we use (imhof etc.) are Gaussian-only by construction. Saddlepoint inversion
needs the CGF of Q, K_Q(theta); for non-Gaussian d there is generally no closed
form for it, and that missing K_Q is the whole obstacle.

A concrete confirmation that the deployed quadratic-form tail tools are
Gaussian-only: the genetics SKAT score test (Wu, Lee, Cai, Li, Boehnke & Lin
2011, Am. J. Hum. Genet. 89(1):82-93) is exactly a quadratic-form score
statistic, and that literature inverts its tail with Davies (exact), Liu-Tang-
Zhang (moment-match), and Kuonen (saddlepoint) -- all assuming the score vector
is asymptotically Gaussian. The same Gaussian-d assumption, in a high-volume
applied setting. (Surfaced by the LitSurvey teammate, 2026-06-25.)

## 2. Our contribution: get K_Q for non-Gaussian d, and it looks novel

The Hubbard-Stratonovich identity makes K_Q available:
    M_Q(theta) = E_{W ~ N(0, I_r)}[ exp(K_d(sqrt(2 theta) W)) ],
exact for theta > 0, evaluated by Gauss-Hermite quadrature (the standard tool
for a Gaussian expectation), then inverted by Lugannani-Rice. Each piece is
textbook: the identity is Hubbard (1959) / Stratonovich (1957), known in ML as
the "Gaussian integral trick"; GH quadrature of a Gaussian expectation is
routine; Lugannani-Rice (1980) is the standard scalar inversion.

The SPECIFIC pipeline -- HS representation of the MGF of a quadratic form in
NON-Gaussian d, GH-quadratured, then LR-inverted -- was not found published in
targeted searching across statistics, econometrics, genetics, and ML. Calibrate
honestly: absence of evidence, not proof of novelty; a careful prior-art check is
needed before claiming it. But it appears to be a genuine recombination filling a
real gap, and it is validated at the CGF level (GH recovered the exact enumerated
K_Q to six digits for M = V_d and M = Sigma_x; see `spike-m2-quad.R`).

## 3. The one true precedent, and how ours differs

Kolassa & Robinson (2011, Annals of Statistics 39(6):3357-3368, arXiv:1203.3106)
is the nearest developed precedent for saddlepoint inference on a multivariate
PERMUTATION test (k-sample one-way design, two-sample multivariate). Read in full
2026-06-25; the mechanism is NOT what the secondary sources implied, and the
difference sharpens our position rather than blurring it.

What they actually do: they do NOT approximate the tail of a quadratic form. They
replace the quadratic form with the likelihood-ratio-like statistic
  Lambda(x-bar) = tau_hat' x-bar - kappa(tau_hat),  tau_hat solving kappa'(tau_hat) = x-bar,
the convex-conjugate (saddlepoint/MLE) transform of the mean. They state plainly
(p. 2) that Lambda "can be approximated to first order by a quadratic form in the
means" but that "it does not seem to be possible to approximate tail probabilities
for quadratic forms with relative errors of order n^{-1}." The cited reason is
Osipov (1981): the Cramer large-deviation tail of a quadratic form has relative
error at best O(n^{-1/4}), because integrating the multivariate Edgeworth
correction over the CURVED (spherical) level set leaves an n^{-1/4} curvature term.
They switch to Lambda to escape that barrier. Their Theorem 1 then approximates
P(Lambda >= lambda) by integrating the formal multivariate saddlepoint density
over the level set {Lambda >= lambda} (a likelihood-ratio level set, NOT an
ellipsoid {Q >= q}), with O(1/n) relative error even when no density exists.
Theorem 2 reduces that multivariate integral to a radial chi-square survival
Q-bar_{d1} times an angular factor G(u), an integral over the unit sphere S_{d1}
evaluated by small Monte Carlo (M = 10 to 1000 sphere points) or Genz (2003)
spherical cubature. They give two equivalent O(1/n) output forms: a multivariate
Lugannani-Rice (eq 8, chi-square + correction) and a multivariate Barndorff-
Nielsen r* (eq 9, chi-square at a modified argument u* = u - log(G(u))/(nu)).
Both are O(1/n); the multivariate r* does NOT buy an extra half-order.

Their Table 3 is the direct evidence for the gap we fill: the classical
Mahalanobis D^2 tail calibrated by chi-square ("Quadratic" row) is 0.0135 vs a
Monte-Carlo 0.0276 at u = 0.5 and 0.0001 vs 0.0006 at u = 0.7 -- badly off -- and
they close "such approximations are not available for the classical quadratic
form." That classical D^2 is exactly what riposte / RItools balanceTest /
Hansen-Bowers (2008) report. So K-R AVOID the quadratic form; we TARGET it. We are
not competing on their statistic; we are filling the cell they declare empty.

How ours differs, precisely: we keep the quadratic form Q as the reported
statistic and get its tail by collapsing Q to a SCALAR exact CGF K_Q (Hubbard-
Stratonovich + Gauss-Hermite, exact identity + numerically exact quadrature),
then a ONE-dimensional Lugannani-Rice. No curved-region multivariate integration
appears, so the Osipov n^{-1/4} barrier does not arise the way it does for K-R.
Two genuinely distinct dimension-killers: K-R reduce to radial x spherical and do
cheap sphere-MC; we Gaussian-smooth to a scalar and do 1-D inversion.

Where the n^{-1/4} comes from, now confirmed from the primary source (RHHQ 1990,
read in full; the technical foundation K-R cite as [11]). RHHQ Corollary 2.1 bounds
the tilted-Edgeworth error over a shifted region B as
  error = (sur(B) n^{1/2} + vol(B)) * O(n^{-d/2 - (s-2)/2}),
splitting a VOLUME (interior) term from a SURFACE (boundary) term, and the surface
term carries an extra n^{1/2}. A FLAT boundary (half-space, the scalar one-sided
tail) dodges the penalty (RHHQ Remark 2, via Chebyshev); a CLOSED CURVED boundary
(the ellipsoid {Q >= q}) keeps it. That surface penalty, specialized to a quadratic
form and worked through Osipov's (1981) constants, is the n^{-1/4}. So the barrier
is a property of the multivariate-Edgeworth-OVER-A-CURVED-REGION method, not of the
scalar random variable Q. K-R's own wording is a hedge -- "does not seem to be
possible" -- not a proven lower bound. Our route never forms that surface integral:
the HS W-integral runs over ALL directions of W, which IS the angular integration
over the degenerate sphere of equally-likely boundary points, and Gauss-Hermite
does it exactly. The curvature is integrated exactly, before any asymptotics. This
also re-explains the inner-Laplace collapse: a leading Laplace of the W-integral
re-imposes the degenerate surface approximation and discards the exact angular
integration; GH is what keeps it exact. (Mechanism cross-checked independently by
the LitSurvey teammate and the RHHQ-primary-source agent, 2026-06-25.)

That the barrier is a method artifact, not intrinsic to Q, is confirmed on two
independent grounds (Osipov-primary-source agent, 2026-06-25; Osipov 1981 itself
was paywalled, so its internal remainder was not read, but the verdict does not
depend on it):
 (i) Same event, two rates. K-R state their LR-like statistic IS a quadratic form
     in the means to first order; Osipov's direct quadratic-form method gives
     n^{-1/4} on it, K-R's re-tilted per-point saddlepoint integration gives O(1/n)
     on it. A rate one method beats by an order on the identical event is the
     method's, not the event's.
 (ii) The exact K-R wording (p.2): the relative error is "of order sqrt(n*lambda)
     n^{-1/4} ... at best of order n^{-1/4}." Because Q is of order 1/n, the
     chi-square-scale thresholds give sqrt(n*lambda) = O(1); and the error DEGRADES
     as the deviation grows -- the opposite of a sound tail method, a tell of a
     method limit.
Subtlety to keep: the degenerate dominating set (rate function constant over the
whole boundary sphere, so the minimizer is the (d-1)-sphere, not a unique point --
Bahadur-Rao/Ney dominating-point theory assumes a unique point) is a real property
of the EVENT, but it is NOT a barrier: K-R pass through exactly this sphere (their
polar reduction integrates over S_{d1}) and still reach O(1/n). So neither curvature
nor the degenerate sphere forces n^{-1/4}; both only defeat the naive direct method.
A related reference if we develop this: a Lithuanian Math. J. paper on large
deviations of spherically distributed vectors when the dominating point degenerates
asymptotically (find exact cite before use). NET: curvature is cleared; the lattice/
discreteness residual below is the honest open question.

Two honest caveats against overclaiming O(1/n):

(1) Escaping the curvature artifact is NOT the same as proving scalar O(1/n). The
scalar LR of Q is O(1/n) only if Q's own law satisfies the scalar saddlepoint
conditions (RHHQ's S.1-S.4 applied to the scalar functional, not the multivariate
region): K_Q analytic and steep up to the pole theta = 1/(2 lambda_max) (holds),
a non-degenerate saddle K_Q'' > 0 (holds), and a Cramer/non-lattice smoothness
condition on the finite-B PERMUTATION law of Q. That last is the real residual
threat -- not curvature, but the DISCRETENESS of Q in very small strata (matched
pairs/triples). It ties to eigenvalue commensurability (sec 5): distinct/
incommensurate eigenvalues give dense Q support -> effectively continuous -> LR
fine; commensurate weights with tiny B give few atoms -> enumerate exactly. A rate
proof would verify S.1-S.4 for the finite-B lattice permutation law of Q.

(2) The honest competitor is Imhof, not Osipov's method. The HS exponent
K_d(sqrt(2 theta) W) = sum_b K_b(...) scales like B, so a leading Laplace of the
W-integral has relative error O(1/B) -- and that leading Laplace IS the Gaussian
weighted-chi-square (Imhof). Therefore Imhof is itself O(1/B) accurate as B grows,
and our GH correction is an O(1/B) refinement that vanishes in the large-B limit.
Against Imhof the RATE in B is the same order; our gain is a better CONSTANT at
small B, because we carry the true skewness and kurtosis of d that Imhof assumes
away. The gain is largest exactly where d is most non-Gaussian: few strata, skewed
scores -- the balance-test regime. State it as "better constant at small B," not
"better rate."

The decisive experiment, refined into two readings:
 - Rate: our error vs B. Plot |p_hat - p_enum|/p_enum against B (or total n) on
   log-log axes. Slope -1 confirms scalar O(1/B); a SHALLOWER slope flags the
   discreteness threat of caveat (1), NOT the Osipov surface penalty (which we do
   not incur). Overlay Lugannani-Rice vs Barndorff-Nielsen r*.
 - Constant: our error vs Imhof's error at fixed small B, across score
   distributions from near-Gaussian to strongly skewed. Should favor us, by a
   margin that widens with the skewness of d.

Strategic answer to "why not just use Lambda, like K-R?": K-R switched statistics
for two reasons -- the direct multivariate saddlepoint on the quadratic form is
degenerate (-> n^{-1/4}), and they had no closed-form scalar CGF for a general
multivariate permutation quadratic form. We have what they lacked: the exact
elementary-symmetric-polynomial CGF K_d plus the HS+GH reduction, so the scalar
route is feasible in our setting where it was not in theirs. Our contribution is to
KEEP the deployed Hansen-Bowers d^2 that practitioners actually report and fix its
computation, rather than asking them to adopt Lambda. A referee will ask why not
Lambda; "backward compatibility with the reported statistic, and the exact CGF is
available here" is the answer to have ready. (Strategic framing from LitSurvey.)

Higher-order asymptotics. Our only asymptotic-in-n error is the final scalar
inversion (HS is exact; GH error is numerical, controlled by node count). Upgrading
that step from Lugannani-Rice (O(1/n)) to the scalar Barndorff-Nielsen r* (the
modified signed root w* = w + (1/w) log(v/w), v = theta_hat sqrt(K_Q''); O(n^{-3/2}))
is a near-free drop-in: same saddle, same K_Q. This is where the scalar reduction
beats K-R: scalar r* gets the extra half-order off-the-shelf, while their
multivariate r* stays at O(1/n). Caveat: in the few-strata regime n is small, so
order labels are guides; trust the enumeration slope. The r* test and the Osipov
test are the same experiment (scale n, overlay LR vs r*).

Two different r*s, do not conflate them. (i) The TAIL r* we would use: the
Barndorff-Nielsen form of the saddlepoint CDF for a statistic with a known CGF,
w* = w + (1/w) log(u/w) read straight off K_Q -- a pure inversion refinement, no
nuisance parameters. This is what K-R generalize to the multivariate case (their
eq 9) and what Brazzale-Davison-Reid (2007) document practically. (ii) The
nuisance-parameter r* (Brazzale's core program, the `hoa` package): the same
algebraic shape but for a parameter of interest in a PARAMETRIC likelihood with
nuisance parameters, where the hard part is the information/nuisance correction
(tangent exponential models, modified profile likelihood, Skovgaard). Our null is
a permutation distribution with no likelihood and no nuisance parameter, so (ii)
does not apply; we use (i), with Brazzale-Davison-Reid (2007) as the how-to and a
small-sample-accuracy precedent. Validate r* against enumeration, not by invoking
the asymptotic order.

Hansen & Bowers (2008) calibrate their omnibus d^2 by chi-square (CLT), never by
saddlepoint, so this genuinely extends their balance test.

DONE 2026-06-25: read Kolassa-Robinson (2011) in full (notes above). STILL TODO:
read Osipov (1981) and Kuonen (1999) directly (DOIs in the download list below).

## 4. Saddlepoint vs Edgeworth: saddlepoint wins in the tails

Edgeworth is a center expansion; in the tails it can go negative and its relative
error can exceed 100%. Saddlepoint (Lugannani-Rice) has bounded relative error
uniformly into the tails. Monti (1993) shows Edgeworth is a re-expanded, degraded
saddlepoint; Robinson (1982) documents saddlepoint beating Edgeworth in the
extreme tails for permutation tests specifically. This justifies the saddlepoint
over the Satterthwaite/Edgeworth/moment-match route (which is what RItools'
`satterthwaite_finite` does).

## 5. Lattice / continuity correction: usually NOT needed here

A continuity correction is driven by the lattice span of Q itself, not by the
discreteness of d's components. Three cases:
 - Distinct/incommensurate eigenvalue weights (the M = Sigma_x case): support is
   dense and irregular, Q is effectively continuous -> continuous Lugannani-Rice,
   NO correction. Kolassa-Robinson treat exactly this and prove O(1/n).
 - Equal/commensurate weights (e.g. M = V_d gives lambda all 1, Q ~ sum of
   squares on a coarser lattice): a correction can measurably help.
 - Rationally-related distinct weights: a fine lattice, correction negligible.
If a correction is wanted, use Daniels' (1987) second correction CC2: half-integer
offset, with u replaced by 2 sinh(s/2) sqrt(K''(s)) (-> the continuous u as the
span -> 0). Verify the exact CC1/CC2 forms against Butler (2007) Sec. 1.2.5
before publishing equations. Practical plan: continuous LR by default; add CC2
only if a coarse-lattice case shows a tail discrepancy against enumeration.

UPDATE 2026-06-26 (rate experiment, dev/route-b-rate-experiment-findings.md): the
framing above (CC driven by Q's lattice span) is NOT the residual we actually found.
Even on dense, incommensurate support (M = Sigma_x-like, the first case above) M2
carries a 3-7% positive tail bias at p=0.05 that a span-driven CC2 cannot touch:
the local span shrinks like 1/N, which would predict a vanishing bias the data
contradict. Controls show it is not lattice span, not max-domination, not curvature
(it FALLS with dimension r), and not the continuous chi^2_r saddlepoint (analytic
chi^2_r LR is +0.5% for all r). It is the gap between continuous LR and the exact
tail of a FINITE-discrete Q, specific to the quadratic form (the scalar sum
saddlepoint on the same orbit is accurate <1%). So CC2 is still right for the genuine
coarse-lattice edge case (single/double stratum), but it is NOT a fix for the
dense-support residual. That residual is small and falling at r>=2 (negligible at
RItools' r=10-50), conservative, and we ship NO correction for it -- decision and
mechanism in the design memo sec "Update 2026-06-26".

## 6. Software (from the code survey)

No package does our non-Gaussian quadratic-form saddlepoint; we write the LR
inversion (~20 lines) and the M_Q quadrature ourselves. Reuse:
 - Quadrature: `mvQuad` (tensor product AND Smolyak sparse Gauss-Hermite, N(0,1)
   pre-weighted) covers M2 (tensor, small r) and M3 (sparse, moderate r);
   `statmod::gauss.quad.prob` for the 1-D rule; for r toward 50, QMC via
   `spacefillr` (MIT, Owen-scrambled Sobol) or `randtoolbox` (BSD-3).
 - Cross-check oracles (tests only): `CompQuadForm::imhof` and
   `survey::pchisqsum(method="saddlepoint")` (Kuonen, Gaussian case);
   `coin::independence_test(teststat="quadratic", distribution=approximate(...))`
   (Monte Carlo permutation p-value on the non-Gaussian case -- a cross-check
   beyond small-orbit enumeration); `RItools::balanceTest` (domain check). A
   further independent cross-check (LitSurvey, 2026-06-25): Imhof-style numerical
   inversion of the CHARACTERISTIC function of Q using the SAME E_W representation
   with imaginary argument, M_Q(it) = E_W[exp(K_d(sqrt(2 it) W))]. It hits the
   same K_Q by a different inversion than Lugannani-Rice, so a tail match between
   the two validates the inversion independently of enumeration. LR stays the
   primary (tail-accurate); CF inversion is a test oracle only.
 - M3 dimension reduction (LitSurvey, 2026-06-25): before quadrature, work in the
   eigenbasis of M (we already rotate to Q = ||rot' d||^2) AND drop near-zero
   eigendirections, which lowers the effective integration dimension r the inner
   Gaussian average runs over. This is the cheapest first move against the tensor-
   GH n^r cost, applied before sparse-grid (Smolyak) or QMC. Does not change the
   pipeline, only the cost.
License: fastperm targets MIT; mvQuad/statmod/CompQuadForm/survey/coin/RItools are
GPL, so keep them in Suggests. For the MIT runtime, vendor a ~30-line
Golub-Welsch Gauss-Hermite for small-r M2; defer the mvQuad-vs-QMC choice to M3.
`spacefillr` (MIT) and `randtoolbox` (BSD) are the runtime-safe QMC options.

## 7. References (status: [V] verified existence/metadata; [F] flagged, verify before citing)

Gaussian-only quadratic-form tails:
 - [V] Imhof (1961), Biometrika 48(3-4):419-426. CF inversion.
 - [V] Davies (1980), JRSS-C 29(3):323-333. AS155, CF inversion.
 - [V] Farebrother (1984), JRSS-C 33(3):332-339. AS204, Ruben series.
 - [V] Liu, Tang & Zhang (2009), CSDA 53(4):853-856. 4-cumulant chi-square match.
 - [V] Wood (1989), Comm. Stat. Sim. Comp. 18(4):1439-1456. F approximation.
Saddlepoint for quadratic forms (Gaussian):
 - [V] Kuonen (1999), Biometrika 86(4):929-935. DOI 10.1093/biomet/86.4.929. OPEN
   ACCESS (EPFL: infoscience.epfl.ch/record/84834/files/860929.pdf). The standard
   saddlepoint-for-quadratic-forms-in-normal-variables cite; our M1 reference.
   NOTE: the printed CGF has a typo (missing zeta in the 2nd-term numerator); use
   the Imhof (1961) eq. 2.3 corrected form. Formulas captured + verified against
   shotGroups/GMMAT, survey/FREGAT, and arXiv:2201.11762 in dev/kuonen-spa-formulas.md.
 - [V] Marsh (1998), Econometric Theory 14(5):539-559. Ratios, normal.
 - [V] Butler & Paolella (2008), Bernoulli 14(1):140-154. Ratios, normal.
Non-Gaussian quadratic forms:
 - [V/F] Khuri & Good (1977), JRSS-B 39(2):217-221; Good (1968) [F: venue
   unverified, cited via Khuri-Good]. CF of d'Ad from the CF of non-normal d.
 - [V] Mathai & Provost (1992), Quadratic Forms in Random Variables, Dekker.
   Mostly Gaussian/elliptical. [F] "Mathai 2005" elliptical extension unverified.
Saddlepoint for permutation tests:
 - [V] Robinson (1982), JRSS-B 44(1):91-101. Scalar; saddlepoint > Edgeworth.
 - [V] Davison & Hinkley (1988), Biometrika 75(3):417-431. Empirical-CGF, scalar.
 - [V] Booth & Butler (1990), Biometrika 77(4):787-796. Conditional double SPA.
 - [V] Abd-Elfattah & Butler (2007), Biometrika 94(3):543-551. Weighted log-rank.
 - [V] Kolassa (2003), Annals of Statistics 31(1):274-293. DOI 10.1214/aos/1046294465.
   Multivariate SPA tail. Kolassa solo; NOT in the K-R 2011 reference list.
 - [V] Kolassa & Robinson (2011), Annals of Statistics 39(6):3357-3368
   (arXiv:1203.3106). DOI 10.1214/11-AOS945. THE precedent. READ IN FULL 2026-06-25:
   they approximate the LR-LIKE statistic Lambda, NOT the quadratic form (which they
   declare intractable-to-O(1/n), citing Osipov 1981). See section 3.
Large-deviation barrier and technical core (from K-R's own reference list, [V]):
 - [V] Osipov (1981), J. Multivariate Analysis 11(1):115-126. DOI
   10.1016/0047-259X(81)90131-1 (MR0618780). Quadratic-form tail O(n^{-1/4}). READ.
 - [V] Robinson, Hoglund, Holst & Quine (1990), Ann. Probab. 18(2):727-753. DOI
   10.1214/aop/1176990856 (MR1055431). READ IN FULL 2026-06-25. The tilted-Edgeworth
   error bound (Thm 1); conditions S.1-S.4 (saddlepoint) / E.1-E.3 (Edgeworth) that
   K-R relabel A1-A4; e2 = their order-2 Edgeworth e_{s-3}, Q1 (n^{-1/2}), Q2 (n^{-1}).
   Cor 2.1: error = (surface * n^{1/2} + volume) * O(...), the curved-boundary penalty.
   (Note: K-R 2011 cites this DOI as ...853; the correct Project Euclid DOI is ...856.)
 - [V] Robinson, Ronchetti & Young (2003), Ann. Statist. 31(4):1154-1169. DOI
   10.1214/aos/1059655909 (MR2001646). Multivariate M-estimate SPA; the LR form.
 - [V] Jing, Feuerverger & Robinson (1994), Biometrika 81(1):211-215. DOI
   10.1093/biomet/81.1.211 (MR1279668). Bootstrap SPA; the r* analog (K-R eq 9).
 - [V] Genz (2003), J. Comput. Appl. Math. 157(1):187-195. DOI
   10.1016/S0377-0427(03)00413-8 (MR1996475). Spherical cubature (K-R G(u)); M3.
Inversion and corrections:
 - [V] Lugannani & Rice (1980), Adv. Appl. Prob. 12(2):475-490. The LR formula. O(1/n).
Higher-order inversion (r*, the modified signed root; O(n^{-3/2}) scalar):
 - [V] Barndorff-Nielsen (1986), Biometrika 73(2):307-322. DOI 10.1093/biomet/73.2.307
   (JSTOR 10.2307/2336207; MR0855891). The scalar r* origin (K-R ref [2]).
 - [V] Barndorff-Nielsen (1991), Biometrika 78(3):557-563. DOI 10.1093/biomet/78.3.557.
   Modified signed log likelihood ratio.
 - [V] Brazzale, Davison & Reid (2007), Applied Asymptotics: Case Studies in
   Small-Sample Statistics, Cambridge. The practical synthesis of the r* / higher-
   order asymptotics program; the how-to reference (with the `hoa` R bundle: cond,
   marg, nlreg) for computing r*. RELEVANCE: relevant for the SCALAR tail r* we
   would use (w* = w + (1/w) log(u/w) off K_Q) and as a small-sample-accuracy
   precedent; NOT for its core content, which is the nuisance-parameter r* for
   PARAMETRIC likelihoods (tangent exponential models, modified profile likelihood,
   Skovgaard) -- a different problem than our randomization-statistic tail (no
   likelihood, no nuisance parameter). Cite as the r* how-to; do not import the
   nuisance machinery. See the two-r* distinction in sec 3.
 - [V/F] Brazzale & Davison (2008), "Accurate parametric inference for small
   samples," Statistical Science 23(4):465-484 [F: volume/pages/DOI from memory,
   verify -- likely DOI 10.1214/08-STS273]. Review of the same program.
 - [V/F] Daniels (1987), Int. Stat. Rev. 55(1):37-48. Lattice CC1/CC2; exact
   equation forms from secondary sources, verify vs Butler (2007) Sec. 1.2.5.
 - [V] Monti (1993), Stat. Prob. Lett. 17(2):131-140. Edgeworth = degraded SPA.
 - [V] Butler (2007), Saddlepoint Approximations with Applications, CUP. Lattice
   corrections in Ch. 1; [F] no standalone "Booth & Butler" lattice article found.
Edgeworth / finite population:
 - [V] Babu & Singh (1985), JMVA 17(3):261-278. Edgeworth, finite-pop, linear.
 - [V] Bhattacharya & Rao (1976/2010). Normal approximation & asymptotic expansions.
Foundations / domain:
 - [V] Hubbard (1959), PRL 3(2):77-78; Stratonovich (1957). The HS transform.
 - [V] Hansen & Bowers (2008), Statistical Science 23(2):219-236
   (arXiv:0808.3857). Omnibus balance d^2, chi-square reference, NOT saddlepoint.
 - [V/F] Wu, Lee, Cai, Li, Boehnke & Lin (2011), Am. J. Hum. Genet. 89(1):82-93.
   SKAT: a quadratic-form score test whose tail is inverted by Davies/Liu/Kuonen,
   all Gaussian -- a concrete high-volume deployment of the Gaussian-only toolkit.
   [F] page range from secondary source (LitSurvey); verify before citing.

## 8. Open items

 - DONE: read Kolassa-Robinson (2011) AND Robinson-Hoglund-Holst-Quine (1990) in
   full (section 3). Curvature mechanism resolved: the n^{-1/4} is RHHQ's surface-
   area penalty (Cor 2.1, extra n^{1/2} on the boundary term), a property of the
   Edgeworth-over-curved-region METHOD, not of Q. Our HS+GH does that surface
   integral exactly, so the penalty does not arise. The real residual threat is the
   DISCRETENESS of Q at small B, not curvature.
 - THE decisive experiment, two readings (section 3): (i) RATE -- our error vs B on
   log-log; slope -1 confirms scalar O(1/B), a shallower slope flags the small-B
   discreteness threat (NOT Osipov, which we do not incur); overlay LR vs r*.
   (ii) CONSTANT -- our error vs Imhof's at fixed small B across near-Gaussian-to-
   skewed scores; should favor us, widening with skewness. Run before any rate claim.
 - Read Osipov (1981) in full for the exact n^{-1/4} constants (RHHQ gives the
   surface mechanism; Osipov specializes it to quadratic forms). Kuonen (1999) for
   the Gaussian M1 CGF/LR formulas. RHHQ S.1-S.4 are the conditions a scalar-Q rate
   proof must verify for the finite-B lattice permutation law.
 - Implement the scalar r* option in the M2 inverter (drop-in: same saddle, same
   K_Q, modified root w* = w + (1/w) log(v/w), v = theta_hat sqrt(K_Q'')).
 - Verify the flagged [F] citations before any manuscript.
 - If this becomes a methods paper: frame as assembling HS + GH + Lugannani-Rice
   (optionally r*) to obtain a non-Gaussian permutation QUADRATIC-FORM tail -- the
   object K-R declare intractable -- citing Kuonen (1999), Osipov (1981), and
   Kolassa-Robinson (2011) as nearest neighbors and extending HB08. The headline
   claim is an accurate tail for the classical Mahalanobis D^2 itself, contingent on
   the relative-error-vs-n experiment confirming the rate.
