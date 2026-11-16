# Saddlepoint tail probability of a linear statistic

Inverts a
[`fastperm_linear_cgf()`](https://bowers-illinois-edu.github.io/fastperm/reference/fastperm_linear_cgf.md)
object into a one-sided permutation tail probability `P(T >= t)`
(default) or `P(T <= t)`, using the Lugannani-Rice saddlepoint
approximation. Probabilities outside the statistic's support are
returned exactly (0 or 1). The small tail is computed in log space, so
`log.p = TRUE` returns accurate values far past the reach of simulation.

## Usage

``` r
saddlepoint_tail(cgf, t, lower.tail = FALSE, log.p = FALSE)
```

## Arguments

- cgf:

  a `fastperm_cgf` object from
  [`fastperm_linear_cgf()`](https://bowers-illinois-edu.github.io/fastperm/reference/fastperm_linear_cgf.md).

- t:

  the value at which to evaluate the tail.

- lower.tail:

  if `FALSE` (default) return `P(T >= t)`; if `TRUE`, `P(T <= t)`.

- log.p:

  if `TRUE`, return the natural log of the probability.

## Value

a probability in `[0, 1]` (or its log if `log.p = TRUE`).

## Examples

``` r
cg <- fastperm_linear_cgf(c(1, 2, 3, 4), c(0, 1, 0, 1), rep(1, 4))
saddlepoint_tail(cg, 6)              # P(T >= 6)
#> [1] 0.2448314
```
