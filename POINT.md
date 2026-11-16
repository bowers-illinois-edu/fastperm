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
a sum of one score per matched set or stratum, such as a variance or a
quadratic form. The design can be matched, stratified, or cluster-randomized.

fastperm computes a p-value from a saddlepoint approximation to the
distribution of the statistic over re-randomizations, rather than by drawing
re-randomizations. When the statistic takes many distinct values, the
saddlepoint's error stays small as a fraction of the p-value even when the
p-value is very small. Speed is not the aim in itself, but it matters in
practice: riposte's combined test, computed from 1,000 brute-force
re-randomizations, takes about 75 seconds per p-value at a million units.

As of November 2026, the package returns saddlepoint p-values for linear and
quadratic statistics when treatment is permuted within strata. The
importance-sampling and Monte Carlo methods for other statistics, and
`perm_test()`, the one function that would choose among these methods for the
user, are planned and not built.

The package is the software for the companion paper in
`~/repos/fastperm-paper`, working title "A Saddlepoint Sensitivity Curve for
Randomization Inference in Observational Studies," aimed at Biometrika. A
reader of the paper should be able to apply its methods to their own data.

Rosenbaum's packages for matched designs (sensitivitymv, senstrat, and
others) compute p-values from a normal approximation. fastperm does not
replace them where that approximation is accurate, such as for M-statistics
in large matched-pair samples. His packages are also faster there. fastperm
does not claim the saddlepoint is accurate when the statistic takes few
distinct values, as with small strata or heavily tied outcomes. There the
approximation can understate the exact p-value. No exact p-value can fall
below the smallest one the permutation distribution attains, yet the
saddlepoint can report values below it. fastperm tests hypotheses. It does
not estimate the size of effects.

## Future work

Sensitivity analysis is an aim for later work. In an observational study the
researcher would also choose Gamma: an unmeasured covariate could make one
unit's odds of treatment as much as Gamma times another's within the same
matched set. fastperm would report the largest p-value that hidden bias of
that size could produce, as a curve over Gamma. Rosenbaum's packages compute
that bound from a normal approximation, one Gamma at a time, and only for
statistics that are a sum of one score per matched set or stratum. fastperm
would compute it from a saddlepoint approximation instead. We aim to prove in
the paper that the saddlepoint gives the whole curve as one formula in Gamma.
`PLAN.md` in the paper repository lists that result as a theorem not yet
proved.

Factorial, stepped-wedge, and regression-discontinuity designs are also
future work.
