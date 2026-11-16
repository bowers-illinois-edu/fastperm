# Exact CGF of a linear statistic under the within-stratum permutation null

Builds the cumulant generating function of the linear statistic
`T = sum_i treatment_i * score_i` referred to the permutation null that
holds each stratum's treated count fixed (a simple random sample of
treated units within each stratum). The returned object carries the
exact permutation mean, variance, and support of `T`, plus closures for
`K(theta)` and its first two derivatives, which
[`saddlepoint_tail()`](https://bowers-illinois-edu.github.io/fastperm/reference/saddlepoint_tail.md)
inverts into a tail probability.

## Usage

``` r
fastperm_linear_cgf(score, treatment, strata)
```

## Arguments

- score:

  numeric vector of per-unit scores `c_i` (for a centred representation,
  the within-stratum-centred values).

- treatment:

  0/1 vector; its per-stratum sums give the treated counts that define
  the permutation null.

- strata:

  factor or vector of stratum labels, the same length as `score`.

## Value

an object of class `fastperm_cgf`: a list with `cgf`, `d1`, `d2`
(closures in `theta`), and the scalars `mean`, `variance`, `supmin`,
`supmax`.

## Examples

``` r
score <- c(1, 2, 3, 4)
cg <- fastperm_linear_cgf(score, c(0, 1, 0, 1), rep(1, 4))
cg$mean
#> [1] 5
```
