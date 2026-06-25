# fastperm

Fast, general-purpose permutation tests for arbitrary test statistics and as-if-randomized designs, with optional Rosenbaum-style sensitivity analysis over the Gamma bias parameter.

## Status

Early development. No functions implemented yet. The methodology and API are being designed. See `docs/memo.md` for the full design document, and the companion paper repository at `~/repos/fastperm-paper/`.

## Motivation

Rosenbaum's sensitivity-analysis software (`sensitivitymv`, `sensitivitymw`, `sensitivityfull`) is powerful but locked into specific designs and test statistics. This package generalizes that machinery: a user supplies a test statistic T, a design (permutation null), and optionally a Gamma, and the package returns a p-value without being restricted to Rosenbaum's specific templates.

## The four-tier dispatch

Internally, `fastperm` routes each call to the fastest method that is accurate for the user's combination of design, test statistic, and Gamma:

- **Tier 0** --- defer to Rosenbaum's own software when the design and statistic match a case his asymptotics cover.
- **Tier 1** --- saddlepoint approximation, for linear or asymptotically-linear statistics. Closed-form, essentially instant, accurate to several digits in the tail.
- **Tier 2** --- cross-entropy adaptive importance sampling, for nonlinear statistics on standard exchangeable or stratified designs.
- **Tier 3** --- brute-force Monte Carlo with Besag-Clifford adaptive stopping, as a fallback.

Sensitivity analyses are made efficient by warm-starting the importance-sampling tilt across the Gamma grid, rather than recomputing from scratch at each Gamma.

## Installation

Not yet installable. When there is code worth installing:

```r
# remotes::install_github("bowers-illinois-edu/fastperm")
```

## Related

- Companion paper repository (local sibling): `~/repos/fastperm-paper/`
- Design memo: `docs/memo.md`
- Intellectual lineage: Robinson (1982, JRSS-B); Mehta, Patel, and Senchaudhuri (1988, JASA); Besag and Clifford (1991); Rosenbaum's sensitivity-analysis line of work; Gao et al. (2024) on the conditional saddlepoint.

## License

MIT. See `LICENSE` and `LICENSE.md`.
