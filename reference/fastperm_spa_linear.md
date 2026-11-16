# Saddlepoint permutation p-value for a linear statistic

Computes the within-stratum permutation p-value of the linear statistic
`T = sum_i treatment_i * score_i` by saddlepoint approximation –
closed-form, with no re-randomization. The permutation null fixes each
stratum's treated count.

## Usage

``` r
fastperm_spa_linear(
  score,
  treatment,
  strata,
  alternative = c("two.sided", "greater", "less")
)
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

- alternative:

  `"two.sided"` (default), `"greater"` (`P(T >= t_obs)`), or `"less"`
  (`P(T <= t_obs)`).

## Value

a list with `p.value`, the observed `statistic`, the permutation `mean`
and `sd`, the `alternative`, and the underlying `cgf` object.

## Examples

``` r
score <- c(1, 2, 3, 4)
fastperm_spa_linear(score, c(0, 1, 0, 1), rep(1, 4), alternative = "greater")$p.value
#> [1] 0.2448314
```
