# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`fastperm` is an R package for fast, general-purpose permutation testing with optional Rosenbaum-style sensitivity analysis over the Gamma bias parameter. The design generalizes Rosenbaum's sensitivity-analysis framework: a user supplies a test statistic, a design (permutation null), and optionally a Gamma, and the package returns a p-value without being locked into Rosenbaum's specific test-statistic/design templates.

Internally, the package is organized around a four-tier dispatch that routes each call to the fastest method that is accurate for the user's combination of design, test statistic, and Gamma. The package is in early development --- no methodology has been implemented. The full design document lives at `docs/memo.md`, and the companion paper/analysis repo is at `~/repos/fastperm-paper/`.

## Companion Paper

The sibling paper and analysis repository lives at `~/repos/fastperm-paper/`. That repo holds the paper, simulation scripts, and applied analyses built on top of this package.

## Build System

Standard R package tooling via `devtools`:

```
Rscript -e 'devtools::load_all()'
Rscript -e 'devtools::document()'
Rscript -e 'devtools::test()'
Rscript -e 'devtools::check()'
```

GitHub Actions runs `R CMD check` on the standard matrix via `.github/workflows/R-CMD-check.yaml`.

## Code Style

See `CLAUDE_CODING.md` in this repo for the full coding-style guidance. Pure R only for now. C++/Rcpp deferred until profiling identifies hotspots. The pure-R implementation serves as the oracle for future C++/GPU rewrites --- keep its tests.

## Architecture

The internal dispatch routes each call to the cheapest method that is accurate for the user's combination of design, test statistic, and Gamma:

- **Tier 0 --- defer to Rosenbaum.** When the design and statistic match a case Rosenbaum's asymptotics cover, call that software and return its answer. Offer a non-asymptotic cross-check as a toggle.
- **Tier 1 --- saddlepoint approximation (SPA).** For linear or asymptotically-linear test statistics, use SPA. Closed-form, essentially instant, accurate to several digits in the tail for moderate n.
- **Tier 2 --- cross-entropy adaptive importance sampling.** For nonlinear statistics on exchangeable or stratified designs, use the CE method to find a good tilted proposal, then importance-weighted Monte Carlo for the tail estimate.
- **Tier 3 --- brute-force with Besag-Clifford adaptive stopping.** Fallback for exotic statistics or designs where SPA and CE are unstable.

Sensitivity analyses are made efficient by warm-starting the optimal tilt `theta*(Gamma)` across the Gamma grid, rather than recomputing from scratch at each Gamma.

## Reference Literature

Key references for the methodology (see `docs/memo.md` section 3 for the full bibliography):

- Robinson, J. (1982). Saddlepoint Approximations for Permutation Tests and Confidence Intervals. JRSS-B 44(1), 91--101.
- Davison, A. C. and Hinkley, D. V. (1988). Saddlepoint Approximations in Resampling Methods. Biometrika 75, 417--431.
- Booth, J. G. and Butler, R. W. (1990). Randomization Distributions and Saddlepoint Approximations in Generalized Linear Models. Biometrika 77, 787--796.
- Abd-Elfattah, E. F. and Butler, R. W. (2007). The Weighted Log-Rank Class of Permutation Tests: P-Values and Confidence Intervals Using Saddlepoint Methods. Biometrika 94, 543--551.
- Gao et al. (2024, arXiv:2407.08911). The Conditional Saddlepoint Approximation for Fast and Accurate Large-Scale Hypothesis Testing.
- Mehta, C. R., Patel, N. R., and Senchaudhuri, P. (1988). Importance Sampling for Estimating Exact Probabilities in Permutational Inference. JASA 83(404), 999--1005.
- Shi et al. (2016, arXiv:1608.00053) and Shi et al. (2019). Cross-entropy method for small p-value estimation in permutation tests.
- Besag, J. and Clifford, P. (1991). Sequential Monte Carlo p-values.
- Hapfelmeier et al. (2023, CSDA). Sequential permutation tests with Type I error control.
- Yu et al. (2011, Biostatistics). Stochastic approximation Monte Carlo (SAMC) for resampling tests.

## Writing and Revision

The historical session HANDOFF.md lives in the paper repo at `~/repos/fastperm-paper/HANDOFF.md`. This repo's design document is `docs/memo.md`.

## Important Notes

- Do not implement methodology without Jake's go-ahead. This package is a scaffolding; Jake will drive the research agenda deliberately. Ask before adding dependencies to Imports.
- Pure R only for now. No `src/`, no Rcpp.
- Keep `Imports` in DESCRIPTION empty until Jake makes the design decisions. `coin` goes in `Suggests` because it is the benchmark comparison. (`stats` was added to `Imports` 2026-06-25 for the saddlepoint core -- base R, the only addition.)
- DOWNSTREAM DEPENDENCY: the `riposte` package (`~/repos/riposte`) calls `fastperm_spa_linear()`, `fastperm_linear_cgf()`, and `saddlepoint_tail()` by name for its fast unadjusted-Cauchy path. Do NOT change those names or signatures without updating riposte and the shared note in `~/repos/riposte/FASTPERM_INTEGRATION.md`. When the planned `perm_test()` front door lands, keep a stable entry point riposte can call.
