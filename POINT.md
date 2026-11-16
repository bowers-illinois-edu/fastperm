# The point of fastperm

fastperm is an R package for researchers who test the sharp null hypothesis
that treatment changed no one's outcome. They compare a test statistic with
its distribution over re-randomizations of treatment. They do this in
randomized experiments, and in observational studies in which, after matching
or stratifying, they treat assignment as random within the matched sets or
strata. The package is also for the authors of other packages that need such
p-values quickly: riposte already calls fastperm's saddlepoint functions.

A researcher should be able to choose the test statistic and the design for
scientific reasons, not because some software happens to cover them. The
statistic can be any function of the data, including statistics that are not
a sum of one score per matched set, such as a variance or a quadratic form.
The design can be matched, stratified, or cluster-randomized.

In an observational study the researcher also chooses Gamma: an unmeasured
covariate could make one unit's odds of treatment as much as Gamma times
another's within the same matched set. fastperm is meant to report the
largest p-value that hidden bias of that size could produce, as a curve over
Gamma. Rosenbaum's packages (sensitivitymv, senstrat, and others) compute
that bound from a normal approximation, one Gamma at a time, and only for
statistics that are a sum of one score per matched set or stratum. fastperm
is meant to compute it from a saddlepoint approximation instead. When the
statistic takes many distinct values, the saddlepoint's error stays small as
a fraction of the p-value even when the p-value is very small. The paper aims
to show that the saddlepoint also gives the whole curve as one formula in
Gamma. Speed would follow from that formula and is not the aim in itself. It
still matters in practice: riposte's combined test, computed from 1,000
brute-force re-randomizations, takes about 75 seconds per p-value at a
million units.

As of November 2026, only part of this exists. The package returns
saddlepoint p-values for linear and quadratic statistics when treatment is
permuted within strata. The Gamma sensitivity curve, the importance-sampling
and Monte Carlo methods for other statistics, and `perm_test()`, the one
function that would choose among these methods for the user, are planned and
not built.

The package is the software for the companion paper in
`~/repos/fastperm-paper`, working title "A Saddlepoint Sensitivity Curve for
Randomization Inference in Observational Studies," aimed at Biometrika. A
reader of the paper should be able to compute its curves on their own data.

fastperm does not replace Rosenbaum's packages where his normal approximation
is accurate, such as M-statistics in large matched-pair samples. His packages
are also faster there. fastperm does not claim the saddlepoint is accurate
when the statistic takes few distinct values, as with small strata or heavily
tied outcomes. There the approximation can understate the exact p-value. No
exact p-value can fall below the smallest one the permutation distribution
attains, yet the saddlepoint can report values below it. fastperm tests
hypotheses. It does not estimate the size of effects.
