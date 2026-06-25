## spa-linear.R --- user-facing saddlepoint p-value for a linear statistic on a
## stratified (within-stratum permutation) design.
##
## This is the convenience that callers (including riposte) use: give the scores,
## the observed treatment, and the strata, and get back the saddlepoint
## permutation p-value for the linear statistic T = sum_i z_i score_i, without
## drawing a single re-randomization. It builds the CGF, evaluates the observed
## statistic, and inverts the requested tail. The two-sided p-value is taken
## about the permutation mean, the continuous-saddlepoint analogue of a two-sided
## permutation p-value.

#' Saddlepoint permutation p-value for a linear statistic
#'
#' Computes the within-stratum permutation p-value of the linear statistic
#' `T = sum_i treatment_i * score_i` by saddlepoint approximation -- closed-form,
#' with no re-randomization. The permutation null fixes each stratum's treated
#' count.
#'
#' @inheritParams fastperm_linear_cgf
#' @param alternative `"two.sided"` (default), `"greater"` (`P(T >= t_obs)`), or
#'   `"less"` (`P(T <= t_obs)`).
#' @return a list with `p.value`, the observed `statistic`, the permutation
#'   `mean` and `sd`, the `alternative`, and the underlying `cgf` object.
#' @examples
#' score <- c(1, 2, 3, 4)
#' fastperm_spa_linear(score, c(0, 1, 0, 1), rep(1, 4), alternative = "greater")$p.value
#' @export
fastperm_spa_linear <- function(score, treatment, strata,
                                alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  cg <- fastperm_linear_cgf(score, treatment, strata)
  t_obs <- sum(as.numeric(treatment) * as.numeric(score))

  p <- switch(
    alternative,
    greater = saddlepoint_tail(cg, t_obs, lower.tail = FALSE),
    less    = saddlepoint_tail(cg, t_obs, lower.tail = TRUE),
    two.sided = {
      dev <- abs(t_obs - cg$mean)                 # two-sided about the mean
      min(1, saddlepoint_tail(cg, cg$mean + dev, lower.tail = FALSE) +
             saddlepoint_tail(cg, cg$mean - dev, lower.tail = TRUE))
    }
  )
  list(p.value = p, statistic = t_obs, mean = cg$mean, sd = sqrt(cg$variance),
       alternative = alternative, cgf = cg)
}
