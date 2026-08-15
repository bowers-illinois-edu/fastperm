# fastperm

Fast, general-purpose permutation tests for arbitrary test statistics and
as-if-randomized designs. The long-term aim includes Rosenbaum-style
sensitivity analysis over the Gamma bias parameter; see Status for what is
built today.

## Status

Version 0.0.0.9001. **Tier 1 (saddlepoint) is implemented for linear and
quadratic statistics under within-stratum permutation. Tiers 0, 2, and 3 are
not, and neither is any Gamma / sensitivity-analysis code.** The package
installs, exports five functions, and carries six test files and a
mathematical vignette.

Implemented and exported:

| Function | What it does |
|---|---|
| `fastperm_linear_cgf()` | Exact CGF of a linear statistic under the within-stratum permutation null |
| `fastperm_linear_cgf_mv()` | Joint CGF of a vector of linear statistics under the same null |
| `saddlepoint_tail()` | Saddlepoint tail probability of a linear statistic |
| `fastperm_spa_linear()` | Saddlepoint permutation p-value for a linear statistic |
| `fastperm_spa_quadratic()` | Saddlepoint tail probability of a quadratic form |

Not implemented: Tier 0 (deferral to Rosenbaum's packages), Tier 2
(cross-entropy importance sampling), Tier 3 (Monte Carlo with Besag-Clifford
stopping), the Gamma tilt, and the warm-started sensitivity curve. The
exponential tilting that does appear in `R/` is the saddlepoint's own
mechanism, not Rosenbaum's Gamma tilt.

The methodology and remaining API live in `docs/memo.md`; the companion paper
repository is `~/repos/fastperm-paper/`.

## Motivation

Rosenbaum's sensitivity-analysis software (`sensitivitymv`, `sensitivitymw`,
`sensitivityfull`) is powerful but locked into specific designs and test
statistics. This package generalizes that machinery: a user supplies a test
statistic T, a design (permutation null), and optionally a Gamma, and gets a
p-value without being restricted to Rosenbaum's templates.

## Usage

```r
# Saddlepoint permutation p-value for a linear statistic, 8 matched pairs.
y      <- c(14.1, 16.0, 7.4, 15.3, 4.3, 11.5, 7.6, 10.9,
            9.5, 15.7, 7.5, 11.9, 12.2, 7.9, 13.4, 13.4)
strata <- rep(1:8, each = 2)
trt    <- rep(c(0, 1), 8)

fastperm_spa_linear(score = y, treatment = trt, strata = strata)$p.value
#> [1] 0.06052
```

Expect `NaN` with a `sqrt(cg$d2(th))` warning when the observed assignment sits
at or near the edge of the permutation support; see Known limitations.

## The four-tier dispatch (the plan)

The intended design routes each call to the fastest method that is accurate
for the user's combination of design, statistic, and Gamma. **Only Tier 1
exists today.**

- **Tier 0** --- defer to Rosenbaum's own software when the design and
  statistic match a case his asymptotics cover. *Not implemented.*
- **Tier 1** --- saddlepoint approximation for linear or asymptotically-linear
  statistics. Closed-form and fast. **Implemented** for linear and quadratic
  statistics.
- **Tier 2** --- cross-entropy adaptive importance sampling for nonlinear
  statistics on exchangeable or stratified designs. *Not implemented.*
- **Tier 3** --- Monte Carlo with Besag-Clifford adaptive stopping, as a
  fallback. *Not implemented.*

Sensitivity analyses are intended to warm-start the tilt across the Gamma
grid rather than recompute at each Gamma. *Not implemented.*

## Known limitations

**Coarse lattices.** The saddlepoint is a continuous approximation to a
discrete permutation distribution. When the statistic takes few distinct
values -- small strata, heavily tied or zero-inflated outcomes -- it is
unreliable in the tail, in the anti-conservative direction, and it can fail
outright at the extremes. On a worked example (10 matched pairs, a
zero-inflated integer outcome, 19 distinct statistic values among 1024
assignments), evaluated at every attainable assignment:

- `fastperm_spa_linear()` returns `NaN` at the support maximum, which is the
  assignment a real effect in a small stratum produces;
- elsewhere in the tail it understates the exact p-value, returning 0.124
  against an exact 0.188 and 0.174 against an exact 0.250.

There is a hard limit underneath this that no approximation can cross. A
discrete permutation distribution has a smallest attainable p-value -- 0.0625
in that example -- and no p-value below it is achievable at any computational
budget. A continuous approximation will nonetheless report values below the
floor on some designs. Where multiplicity corrections need thresholds near or
under that floor, use exact enumeration instead of the saddlepoint.

The quadratic module discusses the same effect: `vignettes/
mathematical-foundations.Rmd` treats the residual at small stratum counts as a
finite-discrete-representation effect rather than curvature, and a
lattice-aware correction is noted there as deliberately out of scope for now.
A dispatcher that chooses Tier 1 automatically will need a support-size check,
not only a check on the CGF's convergence domain: on a coarse lattice the CGF
is perfectly well behaved and the approximation is still wrong.

## Installation

```r
remotes::install_github("bowers-illinois-edu/fastperm")
```

Installs from source with no compiled code.

## Related

- Companion paper repository (local sibling): `~/repos/fastperm-paper/`
- Design memo: `docs/memo.md`
- Mathematical vignette: `vignettes/mathematical-foundations.Rmd`
- Intellectual lineage: Robinson (1982, JRSS-B); Mehta, Patel, and
  Senchaudhuri (1988, JASA); Besag and Clifford (1991); Rosenbaum's
  sensitivity-analysis line of work; Gao et al. (2024) on the conditional
  saddlepoint.

## License

MIT. See `LICENSE` and `LICENSE.md`.
