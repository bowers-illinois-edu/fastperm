# Memo: Fast, General Permutation Testing with Sensitivity Analysis

*Research planning notes from a walking conversation*

**Jake — April 18, 2026**

---

## 1. The core problem

The existing Rosenbaum sensitivity analysis software (`sensitivitymv`, `sensitivitymw`, `sensitivityfull`, and so on) is powerful but locked into specific test statistics and stratified designs. Rosenbaum has proved large-sample properties for particular combinations of design and test statistic, which is excellent work, but it constrains the analyst: if you want a test statistic Rosenbaum has not yet covered, or a design that does not match one of his templates, you are stuck.

The appeal of randomization inference is exactly the opposite: it does not require large samples, and it is agnostic about the test statistic, as long as we know how treatment was assigned (or can act as-if it was randomized conditional on covariates via matching, weighting, stratification, or a design-based argument at a discontinuity).

The goal, then, is a general-purpose permutation engine. The user supplies three things: a test statistic $T$, a design (which specifies the null permutation measure), and a sensitivity parameter $\Gamma$ (optionally). The software returns a $p$-value — for the sharp null of no effect whatsoever, at a specified $\Gamma$ — without being locked into Rosenbaum's specific frameworks.

The computational bottleneck is obvious: for every value of $\Gamma$, a brute-force approach recomputes the full permutation distribution. That is wasteful. Two ideas make this tractable: (a) we only care about the tails of the permutation distribution, not the whole thing; and (b) there is a clean mathematical apparatus — exponential tilting of the permutation measure — for concentrating sampling effort in those tails.

---

## 2. The direct answer: exponential tilting of the permutation measure

The question driving this project is: **which shuffles of the treatment vector $Z$ should we oversample to concentrate $T(Z)$ around a given target region of the distribution?**

The theoretically optimal importance distribution assigns probability to each shuffle $z$ proportional to

$$
q_\theta(z) \;\propto\; \exp\{\theta \cdot T(z)\} \cdot \pi(z)
$$

where $\pi$ is the design's null permutation measure (uniform over the orbit of admissible permutations, or $\Gamma$-biased in Rosenbaum's framework), and $\theta$ is the tilting parameter chosen to solve the saddlepoint equation

$$
\mathbb{E}_\theta[T(Z)] \;=\; t
$$

for whatever target value $t$ we care about — typically the observed test statistic, or a threshold deep in the tail.

Draws from $q_\theta$ are reweighted by the likelihood ratio

$$
\frac{\pi(z)}{q_\theta(z)} \;=\; e^{-\theta \cdot T(z)} \cdot M(\theta)
$$

where $M(\theta) = \mathbb{E}_\pi[e^{\theta T(Z)}]$ is the moment generating function of $T$ under the null. This gives unbiased tail estimates with variance orders of magnitude smaller than crude Monte Carlo.

### 2.1 Why this is beautiful for this project

**Linear test statistics factorize.** For linear $T(Z) = \sum_i c_i Z_i$ — which covers Wilcoxon signed rank, rank-sum, stratified rank statistics, McNemar, and basically every statistic Rosenbaum uses — the tilted measure factorizes. Each $Z_i$ becomes an independently tilted Bernoulli with tilted success probability $p_i(\theta)$. Sampling is direct, no MCMC required. This is Robinson's (1982) construction.

**Nonlinear statistics still tractable.** For nonlinear or complex statistics — robust M-estimators, studentized statistics, custom nonparametric estimands — the tilted measure does not factorize, but we can sample from it via Metropolis–Hastings with pair-swap proposals, or approximate it via the cross-entropy (CE) method, which adaptively searches for a good tilted proposal.

**$\Gamma$-biased null composes cleanly.** Under Rosenbaum's $\Gamma$-biased null, the treatment probability within each stratum is itself exponentially tilted — the odds-ratio bound literally *is* an exponential tilt. So the worst-case $p$-value is a tail probability under a tilted measure, and we want to estimate this efficiently. Composition of tilts is just addition in the canonical parameter: if the $\Gamma$-tilt has canonical parameter $\gamma$ and the importance-sampling tilt has canonical parameter $\theta$, the combined measure is tilted at $\gamma + \theta$. For a linear statistic in a stratified design, the combined measure is still a product of independent tilted Bernoullis per stratum. No MCMC, direct sampling.

---

## 3. The relevant literature, in three strands

### 3.1 Saddlepoint approximation (SPA) for permutation tests

Saddlepoint approximation gives a closed-form tail approximation with relative error of order $O(n^{-3/2})$ — materially better than the central limit theorem's $O(n^{-1/2})$ — and, crucially, the SPA error is *uniformly bounded in the tail* rather than deteriorating as the CLT does. For most linear statistics and most reasonable sample sizes, SPA is accurate to three or four digits in the tail.

Key references:

- Robinson, J. (1982). "Saddlepoint Approximations for Permutation Tests and Confidence Intervals." *JRSS-B* 44(1), 91–101. The canonical reference — exactly this problem.
- Davison, A. C. and Hinkley, D. V. (1988). "Saddlepoint Approximations in Resampling Methods." *Biometrika* 75, 417–431.
- Booth, J. G. and Butler, R. W. (1990). "Randomization Distributions and Saddlepoint Approximations in Generalized Linear Models." *Biometrika* 77, 787–796.
- Abd-Elfattah, E. F. and Butler, R. W. (2007). "The Weighted Log-Rank Class of Permutation Tests: P-Values and Confidence Intervals Using Saddlepoint Methods." *Biometrika* 94, 543–551.
- Gao et al. (2024, arXiv:2407.08911). "The Conditional Saddlepoint Approximation for Fast and Accurate Large-Scale Hypothesis Testing." Rigorously extends SPA to conditional randomization tests with relative-error guarantees in the tail.
- Newer, H. A. (2025, *Scientific Reports*). Double-saddlepoint for rank-based independence tests with clustered data. Empirical demonstration of SPA delivering permutation-accuracy at a tiny fraction of the cost.

### 3.2 Importance sampling for permutation tests directly

This is even more on-the-nose for our project.

- Mehta, C. R., Patel, N. R., and Senchaudhuri, P. (1988). "Importance Sampling for Estimating Exact Probabilities in Permutational Inference." *JASA* 83(404), 999–1005. Precisely our problem, published in our target journal. Worth reading closely to see how the JASA reviewers framed the contribution.
- Shi et al. (2016, arXiv:1608.00053) and Shi et al. (2019) applied the cross-entropy method to small $p$-value estimation in permutation tests. The CE method iteratively minimizes the cross-entropy divergence between the current proposal and the optimal proposal. Reported order-of-magnitude gains for small $p$-values.
- Hu, J. and Su, Z. (2008). Adaptive resampling algorithms for bootstrap tail probabilities — conceptually adjacent, demonstrates the CE-method pattern in a resampling context.

### 3.3 Adaptive and sequential stopping

Complementary to importance sampling: stop permuting once you know the answer.

- Besag, J. and Clifford, P. (1991). The original adaptive permutation scheme. Stop once you have accumulated enough "successes" (permuted statistics at least as extreme as observed) to know you are not in the tail.
- Hapfelmeier et al. (2023, *CSDA*). Modern sequential permutation tests with formal Type I error control, implemented in the `rfvimptest` R package.
- Yu et al. (2011, *Biostatistics*). SAMC (stochastic approximation Monte Carlo) for resampling tests. Reports $100\times$ to $500{,}000\times$ efficiency gains for small tail probabilities — basically a proof-of-concept for the software we want to build.
- Che et al. (2014, PMC4070098). Adaptive permutation for GWAS with explicit precision targeting.

---

## 4. Recommended software architecture: one R package, four tiers

After discussion, the decision is to build a single R package rather than splitting the sensitivity analysis into a separate downstream package. The reasoning: once we have the tilted-sampling machinery in place, iterating over $\Gamma$ is not a separate computational problem — it is the same problem with a loop wrapper. Splitting would create two APIs to learn and maintain, for little benefit. Some users will use the package just for fast permutation testing (a `coin` replacement for workflows that iterate permutation tests); others will use it for $\Gamma$ sensitivity analysis; the core code is the same.

The internal dispatch logic routes each call to the fastest method that gives accurate answers for the user's combination of design, test statistic, and $\Gamma$.

**Tier 0 — defer to Rosenbaum.** If the design, test statistic, and null combination matches a case for which Rosenbaum (or Fogarty, or others in that lineage) already has asymptotic results, we call their software and return that answer. We offer a non-asymptotic cross-check as a toggle so that small-$n$ users can verify the asymptotics. No point reinventing the wheel.

**Tier 1 — saddlepoint approximation.** For linear or asymptotically-linear test statistics, we use SPA directly. This is closed-form, essentially instant, and accurate to three or more digits in the tail for $n \gtrsim 20$. This tier probably covers 70%+ of practical use — Wilcoxon, rank sum, stratified rank statistics, signed-rank, McNemar, Mantel–Haenszel, etc.

**Tier 2 — cross-entropy adaptive importance sampling.** For non-linear statistics on standard (exchangeable or stratified) designs, we use the CE method to find a good tilted proposal, then importance-weighted Monte Carlo for the tail estimate. Order-of-magnitude speedup over crude permutation.

**Tier 3 — brute-force with Besag–Clifford adaptive stopping.** Fallback for anything weird: exotic test statistics, unusual designs, cases where SPA or CE are unstable. Adaptive stopping makes this far cheaper than naive Monte Carlo permutation for the common case where the result is clearly significant or clearly not.

### 4.1 Sensitivity-analysis loop: warm-starting over $\Gamma$

Critical for the sensitivity use case: we do not recompute from scratch for each $\Gamma$. The optimal tilt $\theta^*(\Gamma)$ is a smooth function of $\Gamma$. Warm-starting the tilting parameter (and the CE proposal when Tier 2 is in use) across $\Gamma$ values should give the biggest single practical speedup relative to naive iteration.

Concretely: if we want the sensitivity curve over $\Gamma \in [1.0, 5.0]$ at 41 grid points, we do one full solve at $\Gamma = 1.0$, then propagate the solution forward with small perturbations. The per-step cost becomes marginal.

---

## 5. A stress-test to flag: constrained-orbit designs

The tilting approach works cleanly when the permutation orbit admits a tractable MGF for $T(Z)$. For exchangeable designs, stratified designs, and matched-pair designs, this is standard. For some designs — network designs, cluster-randomized trials with unusual cluster-level constraints, and restricted permutations satisfying multiple linear constraints — the permutation orbit is constrained in ways that make the MGF unavailable in closed form. You then need the Diaconis–Sturmfels machinery (Markov bases from algebraic statistics) to sample from the constrained orbit at all, let alone tilt.

If the paper claims to cover "any as-if-randomized design," we will need to decide how to handle the constrained-orbit case. Two options:

1. **Scope the paper honestly:** cover exchangeable, stratified, and matched designs; note that constrained-orbit designs require additional machinery; leave for future work.
2. **Embrace the harder case:** implement the Diaconis–Sturmfels sampler for specific families of constrained designs (e.g., contingency tables with fixed margins) and combine it with the tilting. This may be the novel methodological contribution that distinguishes us from Mehta et al. 1988 and subsequent work.

Lazzeroni and Lange, Diaconis and Sturmfels, and the algebraic-statistics literature on Markov bases are the relevant references. Worth at minimum a careful scoping conversation before settling the paper's claimed generality.

---

## 6. Unconventional angle worth considering

The obvious computational wins are C++ backend and GPU parallelism for the inner permutation loop. Those matter, but the bigger win is algorithmic: SPA plus tilting plus warm-starting is the real speedup, and may well mean that GPUs are unnecessary for most of the use cases that matter.

A more surprising framing: for the large class of problems where SPA gives a closed-form tail probability (linear statistics, standard designs), the sensitivity curve $\Gamma \mapsto p(\Gamma)$ is *effectively analytic*. That is, we may be writing a paper whose headline is "we don't need permutations at all for the most common cases — the sensitivity curve has a closed-form approximation." That is a cleaner, stronger JASA story than "we permute faster." Worth considering which framing to lead with.

---

## 7. Illustrative example scenarios

The paper needs a mix of scenarios demonstrating the speedup in contexts where it matters. Below are ten candidates; the suggestion is to use five for the paper, picked to span matched-pairs, binary outcomes, clusters/strata, nonlinear statistics, and at least one real applied example.

### Scenario 1 — Stratified matched pairs, Wilcoxon signed-rank, moderate $n$

$n \approx 100$ matched pairs, stratified by region or age band, continuous outcome. Brute-force permutation feasible but slow; SPA instant. Sweep $\Gamma \in [1.0, 5.0]$, plot sensitivity curve. Headline: exact-accuracy $p$-values in milliseconds rather than seconds. This is the workhorse example.

### Scenario 2 — McNemar test for paired binary outcomes

$n \approx 200$ matched pairs, binary response, treatment vs. control. Permutation space is $2^{200}$ — cannot enumerate. Adaptive stopping plus SPA crushes brute-force importance sampling. Highlights the "we do not need that many permutations" angle.

### Scenario 3 — Stratified two-sample rank test, small strata

$\approx 15$ strata, 3–5 per stratum, total $n \approx 50$–$70$. Fixed-margin constrained permutation space per stratum is huge relative to $n$. SPA works beautifully here because even with tiny strata, the CGF-based approximation is accurate. Brute force via `coin` would struggle to iterate $\Gamma$.

### Scenario 4 — Regression discontinuity with local permutation test

Running variable, cutoff, local randomization assumption within a bandwidth. $n$ in-window $\approx 100$–$300$. Vary the bandwidth and re-test; each iteration is a permutation test. SPA plus warm-start makes this interactive. Figure: $x$-axis is bandwidth, $y$-axis is $p$-value, overlay compute time — our method flat, `coin` spiky.

### Scenario 5 — Cluster-randomized trial, Wilcoxon on cluster summaries

30 clusters, 15 treated, 15 control, individual-level outcomes aggregated to cluster medians. Permutation space is $\binom{30}{15}$ — large. SPA on the cluster-level rank statistic is instant. Highly relevant to education and public-health trials, which is exactly the development-economics audience that often mis-labels what they are doing as "randomization inference."

### Scenario 6 — Matched triplets, signed-rank, canonical Rosenbaum sensitivity

Each treated unit matched to two controls. Continuous outcome. Signed-rank test. Iterate $\Gamma \in [1.0, 10.0]$, plot $p$-value trajectory. The canonical Rosenbaum case — we want to show that our package handles it interactively, not in minutes.

### Scenario 7 — Robust M-estimator or trimmed mean (nonlinear statistic)

Skewed outcome (income, medical cost). The ATE test is underpowered; trimmed mean or M-estimator is better. Use CE-adaptive importance sampling (Tier 2). Show we beat `coin`'s approximation on both speed and accuracy relative to asymptotic normality.

### Scenario 8 — Fisher-exact alternative for small $n$ with ties, iterated over $\Gamma$

$n \approx 30$ total, imbalanced, binary outcome, many ties. Our permutation test with SPA recovers Fisher's exact $p$-value (validation), then extends to a sensitivity curve over $\Gamma$. Fisher's exact by itself does not give you a sensitivity analysis.

### Scenario 9 — Simulation study: precision vs. compute budget

Type I error and power under the null, varying the number of permutation samples across methods: our adaptive-stopping-plus-SPA, `coin`, and brute force. Show we hit a target precision (say, $\mathrm{SE} < 0.01 \cdot p$) with dramatically fewer draws. Likely a supplement or a self-contained section.

### Scenario 10 — Real applied example

The interrogation-deception-ban-and-clearance-rates study, or a matched observational study from epidemiology, education, or comparative politics. Run the test, iterate over $\Gamma$, report sensitivity curve, highlight computational savings. Anchors the paper in a real problem — crucial for JASA Applications or a political science venue.

### 7.1 Starting strategy

Start with one or two semi-simulated examples for clean speed benchmarks (reproducible, no confounding narrative). Then pivot quickly to real data — semi-simulation alone will feel thin for a JASA paper unless it is doing something clever like robustness across a parameter grid.

---

## 8. Target venue: to be decided

Several plausible venues, each with trade-offs:

- **JASA (Theory and Methods or Applications).** Strongest methods audience, precedent with Mehta–Patel–Senchaudhuri 1988. High-prestige but slow cycle and demanding reviewers.
- **Political Analysis.** Our core disciplinary audience, friendlier cycle, but the methodological contribution may read as too technical for a political-science venue.
- **Development economics venue (AEJ Applied, *Journal of Development Economics*).** Right audience for practitioners who confuse cluster-randomized trials with randomization inference. But a methods paper may be a poor fit for their editorial taste.
- **Observational Studies.** Rosenbaum-adjacent venue, but smaller reach and impact than JASA.

**TO DO:** revisit venue choice after the first draft is written — the framing of the draft will partly settle this. Worth being open to a two-paper strategy: a methods paper in JASA and a shorter applied piece in *Political Analysis* or a development venue using the same package.

---

## 9. Future directions: sensitivity analysis for null and near-null effects

Flagging this here so it is not lost, even though it likely belongs in a follow-up paper rather than this one.

The standard Rosenbaum sensitivity analysis is directional: it asks, "If there is unobserved confounding of magnitude $\Gamma$, could my *significant* result be overturned?" It answers a question about affirmative findings. But there is a mirror-image question that is equally valid and largely unformalized in the literature: **"If my result is null or small, could unobserved confounding be hiding a real treatment effect?"**

This matters for two reasons. First, analysts routinely over-interpret non-significant $p$-values as evidence for the null — a basic misreading of frequentist logic. Second, even when they do not do that, there is often a substantive claim at stake: "this policy had no effect," or "this intervention did not move the outcome." Substantively null claims deserve sensitivity analysis too. Unobserved confounders can bias point estimates toward zero as easily as away from zero, especially in observational settings with selection on unobservables operating in opposition to the treatment effect.

Rosenbaum's existing framework does not directly address this — credit to his extraordinary analytical contributions, but the framework is built for the affirmative-finding case. The extension would ask: *under what magnitude of unobserved confounding would my null result become significant* (in the direction that substantive theory predicts)? That is the dual sensitivity question.

Argument for deferring to a separate paper:

1. This paper is solving a clean computational problem, and mixing in a new conceptual contribution will muddy it.
2. The null-effects sensitivity question requires new thinking about what $\Gamma$ means and how to interpret the result, not just new machinery.
3. The computational infrastructure from this paper is exactly the enabler we need to write that follow-up paper well. Once the tilted-sampling machinery is in place, asking the dual question is straightforward to implement.

Suggested handling for this paper: one sentence in the introduction or discussion acknowledging this motivation — something like, "One application is sensitivity analysis for null and near-null effects, which we pursue in future work" — and then leave it alone. Keep detailed notes on the follow-up paper's conceptual content separately.

---

## 10. Immediate next steps

1. Pull full PDFs of Robinson (1982), Mehta–Patel–Senchaudhuri (1988), Gao et al. (2024) conditional saddlepoint, and Yu et al. (2011) SAMC. Work through the proofs and algorithm descriptions in detail.
2. Draft a package API sketch: function signatures, argument structure, return types. Decide what the user interface looks like for the four tiers — ideally the user does not have to care which tier fires.
3. Prototype SPA for Wilcoxon signed-rank in a stratified matched-pair design as the minimum viable example. Validate against brute-force permutation on a small problem. Time both.
4. Prototype the warm-start logic across $\Gamma$. This is the biggest practical speedup and the one most novel to our setting.
5. Pick three of the ten example scenarios to implement first. Suggested: Scenario 1 (stratified matched pairs), Scenario 6 (matched triplets, canonical Rosenbaum), and Scenario 10 (real data).
6. Decide on the constrained-orbit question: do we scope it out, or embrace Diaconis–Sturmfels as the novel methodological contribution?
7. Defer venue decision until after a draft exists.
8. Keep separate notes on the null-effects sensitivity paper — do not let it leak into this paper.

---

*— end of memo —*
