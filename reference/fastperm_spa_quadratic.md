# Saddlepoint tail probability of a quadratic form under the permutation null

Computes the within-stratum permutation tail probability of the
quadratic form `Q = (T - mu)' M^{-1} (T - mu)`, where
`T_c = sum_i treatment_i * score_{i,c}` is a vector of linear
statistics, `mu` is its exact permutation mean, and `M` is a metric.
This is the omnibus covariate-balance statistic of Hansen and Bowers
(2008) when `M` is the permutation covariance of `T`. The permutation
null fixes each stratum's treated count.

## Usage

``` r
fastperm_spa_quadratic(
  scores,
  treatment,
  strata,
  metric = "cov",
  method = c("gaussian", "saddlepoint"),
  nodes = 32L
)
```

## Arguments

- scores:

  numeric matrix of per-unit scores (units x representations).

- treatment:

  0/1 vector; its per-stratum sums fix the treated counts.

- strata:

  factor or vector of stratum labels, the same length as the rows of
  `scores`.

- metric:

  either the string `"cov"` (use the exact permutation covariance of
  `T`, giving the Hansen-Bowers omnibus statistic) or a symmetric
  positive-definite `p x p` matrix (for example a covariate covariance).

- method:

  `"gaussian"` (M1: the weighted-chi-square tail under a Gaussian
  approximation to the centred statistic, inverted by Lugannani-Rice) or
  `"saddlepoint"` (M2: the non-normal permutation saddlepoint, which
  carries the true skewness and higher cumulants of the permutation law
  into the tail).

- nodes:

  number of Gauss-Hermite nodes per dimension for the `"saddlepoint"`
  method (the Hubbard-Stratonovich integral is an `r`-dimensional
  Gaussian expectation evaluated by tensor quadrature). Ignored for
  `"gaussian"`.

## Value

a list with `p.value`, the observed `statistic` (`Q`), the `rank` of the
metric, the chi-square `eigenvalues` (weights), the permutation `mean`
of `Q`, the resolved `metric`, and the `method`. For
`method = "saddlepoint"` it also carries `p.value.rstar`, the
Barndorff-Nielsen `r*` higher-order refinement of `p.value`.

## Details

The metric enters only through the rotation `R` with `R R' = M^{-1}`:
exactly, `Q = || R'(T - mu) ||^2`, a sum of squares whose
Gaussian-approximation weights are the eigenvalues of the rotated
permutation covariance.

The `"saddlepoint"` method uses the exact permutation CGF of the centred
statistic via the Hubbard-Stratonovich identity
`M_Q(theta) = E_{W ~ N(0, I_r)}[exp(K_e(sqrt(2 theta) W))]`, evaluated
by tensor Gauss-Hermite quadrature and inverted by Lugannani-Rice. The
tensor rule costs `nodes^r`, so it is feasible only for small `r` (the
number of representations); larger `r` needs the sparse-grid/QMC engine,
which is not yet implemented. Below `E[Q]` the method falls back to the
Gaussian-d tail, which is accurate there and not a rejection region.

## References

Hansen, B. B. and Bowers, J. (2008). Covariate balance in simple,
stratified and clustered comparative studies. Statistical Science 23,
219–236.

## Examples

``` r
scores <- cbind(c(1, 2, 3, 4, 5, 6), c(1, 1, 2, 2, 5, 9))
z <- c(0, 0, 0, 1, 1, 1)
# omnibus balance test (metric = the permutation covariance)
fastperm_spa_quadratic(scores, z, rep(1, 6))$p.value
#> [1] 0.1352993
# the non-normal permutation saddlepoint
fastperm_spa_quadratic(scores, z, rep(1, 6), method = "saddlepoint")$p.value
#> Warning: observed Q is at or beyond the saddlepoint's reliable range; falling back to the Gaussian-d tail.
#> [1] 0.1352993
```
