# Joint CGF of a vector of linear statistics under within-stratum permutation

Builds the joint cumulant generating function `K_T(s)` of the vector
statistic `T_c = sum_i treatment_i * score_{i,c}`, referred to the
permutation null that fixes each stratum's treated count. `K_T(s)` is
evaluated by projecting the score matrix onto `s` and applying the
univariate elementary-symmetric CGF, so the cost is that of one 1-D CGF
per evaluation. The exact permutation mean vector and covariance are
returned in closed form; `grad K_T(0)` and `Hess K_T(0)` reproduce them.

## Usage

``` r
fastperm_linear_cgf_mv(scores, treatment, strata)
```

## Arguments

- scores:

  numeric matrix of per-unit scores (units x representations).

- treatment:

  0/1 vector; its per-stratum sums fix the treated counts.

- strata:

  factor or vector of stratum labels.

## Value

an object of class `fastperm_cgf_mv`: a list with `cgf` (a closure of
the length-d vector `s`), the exact `mean` vector and `cov` matrix, and
`d`.
