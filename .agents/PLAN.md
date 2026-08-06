# ertte development plan

This document tracks scoped-out future development for ertte -- work
that's been thought about but not done, or deliberately deferred. It is
not a changelog: once an item here is completed, its write-up should
move to [.agents/HISTORY.md](HISTORY.md) and be removed from this file
rather than marked "done" in place.

The original design/scoping issue,
[#1](https://github.com/djnavarro/ertte/issues/1), remains open --
every other item it scoped, and every issue filed during
stress-testing (#3-#12), has shipped (see `HISTORY.md`). One other
issue is currently open alongside it, tracked below
([#15](https://github.com/djnavarro/ertte/issues/15)).

## Revert CI's temporary erplots branch pin (issue #15)

PR #14 (implementing `er_predict_survival.ertte_model()` for issue
#13) pinned `.github/workflows/R-CMD-check.yaml` and
`.github/workflows/test-coverage.yaml` to install
`github::djnavarro/erplots@feat/er-tte-core-scaffolding` instead of
erplots' default branch, since that branch is currently the only
place `erplots::er_tte()`/`er_predict_survival()` exist --
`tests/testthat/test-er-methods.R`'s `er_tte()`-grammar integration
tests would otherwise `skip_if_not(exists(...))` silently in CI rather
than actually running. `pkgdown.yaml` was left installing the default
branch (unaffected, since it doesn't run tests).

Revert both pinned workflow files back to
`github::djnavarro/erplots` (no branch qualifier) once
`feat/er-tte-core-scaffolding` merges to erplots' default branch.

## `ertte_rmst()`'s CI bounds aren't constrained to `[0, tau]`

`ci_lower`/`ci_upper` on both engines are symmetric Wald intervals on
the RMST scale, unlike `ertte_predict()`'s survival-probability
intervals, which are naturally bounded to `[0, 1]` by their CDF
back-transform. A Wald interval on RMST has no equivalent
transform-induced bound, so it's possible (for a wide interval near
one of RMST's structural limits) for `ci_lower < 0` or `ci_upper >
tau`. Not yet addressed -- would need either a bounded reparameterisation
(analogous to the survival-probability case) or an explicit `pmax`/`pmin`
clamp (cruder, and arguably misleading about the interval's actual
coverage).

## VPC scalar-reduction censoring: IPCW/pseudo-value weighting

`er_simulate.ertte_model()`'s landmark/RMST VPC support (see
`AGENTS.md`) uses a complete-case convention: a simulated replicate
with an outcome that's ambiguous relative to `landmark_time`/`tau`
(censored strictly before it) becomes `NA`, excluded from
`er_vpc_add_simulated()`'s `mean(..., na.rm = TRUE)` aggregation. This
is a reasonable default but not bias-corrected for informative
censoring.

Two corrections were investigated and both rejected *for now*, not
because they're wrong ideas, but because erplots' VPC contract has no
slot for what they need:

- **IPCW** (inverse-probability-of-censoring weighting, Robins-style)
  needs a per-replicate weight column to reach
  `er_vpc_add_simulated()`'s aggregation -- checking erplots'
  `R/er-vpc-layer.R` confirmed every observed/simulated summary is
  currently a plain unweighted `mean(..., na.rm = TRUE)`, with no
  weight column anywhere in the contract.
- **Pseudo-observations** (Andersen-Perme-style) would fit the
  contract mechanically (reduces each censored outcome to a single
  value compatible with ordinary averaging) but needs a leave-one-out
  KM/RMST jackknife recomputation per simulated replicate -- real
  computational cost (`nsim` separate leave-one-out passes over `n`
  observations) for a benefit that's conceptually murky here anyway:
  pseudo-values correct for censoring bias in an *unknown* population,
  but every simulated replicate is already a fully known draw from the
  fitted model.

Revisit only after an erplots-side contract change adds a weighted
aggregation path (mirroring the `predict_args`/`summary_args`/
`simulate_args` splicing pattern erplots#11 already added). In the
meantime, `censor_time` (via `simulate_args`) is the existing lever for
reducing how much ambiguity arises in the first place.

## Cox RMST SE: the independent-increments simplification

`ertte_rmst.ertte_coxph()`'s SE (`.ertte_rmst_pfun_delta()`) treats
`Var[H(t|x)]`'s increments as accumulating independently across jump
times, mirroring the classic Greenwood-formula assumption -- not
exactly right, since the coefficient-uncertainty component of `H(t|x)`
is really a single random direction shared across every `t`, not a sum
of independent per-jump increments. Cross-validated against a
bootstrap across two contrasting covariate profiles and found to track
it substantially more closely than the alternatives considered (see
`HISTORY.md`), but not a fully principled derivation. No concrete
improvement is planned -- noting it here so a future revisit doesn't
have to rediscover the caveat from scratch.
