test_that("ertte_rmst() on ertte_aft returns a tidy tibble with sane RMST values", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  rmst <- ertte_rmst(mod, ertte_data[1:5, ], tau = c(30, 60, 90))
  expect_s3_class(rmst, "tbl_df")
  expect_equal(nrow(rmst), 15L)
  # RMST is an area under a survival curve bounded in [0, 1], so it can
  # never exceed its own horizon
  expect_true(all(rmst$fit_rmst >= 0 & rmst$fit_rmst <= rmst$tau))
  expect_true(all(rmst$ci_lower <= rmst$fit_rmst))
  expect_true(all(rmst$ci_upper >= rmst$fit_rmst))
  # RMST should increase with tau, for a fixed subject
  one_subject <- rmst[rmst$id == rmst$id[1], ]
  expect_true(all(diff(one_subject$fit_rmst) >= 0))
})

test_that("ertte_rmst() on ertte_coxph returns a tidy tibble with sane RMST values", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  rmst <- ertte_rmst(mod, ertte_data[1:5, ], tau = c(30, 60, 90))
  expect_s3_class(rmst, "tbl_df")
  expect_equal(nrow(rmst), 15L)
  expect_true(all(rmst$fit_rmst >= 0 & rmst$fit_rmst <= rmst$tau))
  expect_true(all(rmst$ci_lower <= rmst$fit_rmst))
  expect_true(all(rmst$ci_upper >= rmst$fit_rmst))
  one_subject <- rmst[rmst$id == rmst$id[1], ]
  expect_true(all(diff(one_subject$fit_rmst) >= 0))
})

test_that("ertte_rmst() on ertte_coxph matches a hand-computed step-function integral", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  prof <- ertte_data[1, , drop = FALSE]
  tau <- 90
  rmst <- ertte_rmst(mod, prof, tau = tau)

  sf <- survival::survfit(mod, newdata = prof)
  manual_rmst <- function(time, surv, tau) {
    keep <- time < tau
    t_use <- c(0, time[keep])
    s_use <- c(1, surv[keep])
    sum(s_use * diff(c(t_use, tau)))
  }
  expect_equal(rmst$fit_rmst, manual_rmst(sf$time, sf$surv, tau))
})

test_that("ertte_rmst()'s Cox SE is not just the population-level Greenwood term", {
  # regression test pinning down the fix described in R/ertte-rmst.R:
  # the profile-specific SE should differ meaningfully from `survival`'s
  # own population-level `survmean()` (accessed here only to demonstrate
  # the two disagree, not as a dependency of the package code itself)
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  prof_hi <- ertte_data[which.max(ertte_data$aucss), , drop = FALSE]
  tau <- 90

  rmst <- ertte_rmst(mod, prof_hi, tau = tau)
  sf <- survival::survfit(mod, newdata = prof_hi)
  naive_se <- survival:::survmean(sf, rmean = tau)$matrix["se(rmean)"]

  se_new <- (rmst$ci_upper - rmst$fit_rmst) / stats::qnorm(.975)
  expect_false(isTRUE(all.equal(unname(se_new), unname(naive_se))))
})

test_that("ertte_rmst() output has the expected columns, no leftovers", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  nd <- ertte_data[1:5, ]
  rmst <- ertte_rmst(mod, nd, tau = 90)
  expect_equal(
    setdiff(names(rmst), names(nd)),
    c("tau", "fit_rmst", "ci_lower", "ci_upper")
  )
})

test_that("ertte_rmst() validates tau", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  nd <- ertte_data[1:5, ]
  expect_error(ertte_rmst(mod, nd, tau = -10), "tau")
  expect_error(ertte_rmst(mod, nd, tau = 0), "tau")
  expect_error(ertte_rmst(mod, nd, tau = NA_real_), "tau")
  expect_error(ertte_rmst(mod, nd, tau = "90"), "tau")
  expect_error(ertte_rmst(mod, nd, tau = numeric(0)), "tau")
})

test_that("ertte_rmst() validates conf_level", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  nd <- ertte_data[1:5, ]
  expect_error(ertte_rmst(mod, nd, tau = 90, conf_level = 1.5), "conf_level")
})

test_that("ertte_rmst() defaults newdata to the fitted data", {
  mod_aft <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  mod_cox <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_equal(nrow(ertte_rmst(mod_aft, tau = 90)), nrow(ertte_data))
  expect_equal(nrow(ertte_rmst(mod_cox, tau = 90)), nrow(ertte_data))
})

test_that("ertte_rmst() works across all supported AFT distributions", {
  for (dd in c("exponential", "weibull", "lognormal", "loglogistic")) {
    mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data, dist = dd)
    rmst <- ertte_rmst(mod, ertte_data[1:3, ], tau = 90)
    expect_true(all(is.finite(rmst$fit_rmst)))
    expect_true(all(rmst$fit_rmst >= 0 & rmst$fit_rmst <= 90))
  }
})

test_that("ertte_rmst() on ertte_coxph warns when tau exceeds the observed follow-up range", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  max_time <- max(ertte_data$time)
  expect_warning(
    ertte_rmst(mod, ertte_data[1, , drop = FALSE], tau = max_time + 100),
    "exceeds the last observed follow-up time"
  )
  expect_no_warning(ertte_rmst(mod, ertte_data[1, , drop = FALSE], tau = 30))
})

test_that("ertte_rmst() on ertte_coxph inherits the all-censored-Cox guard", {
  dat_censored <- ertte_data[1:60, ]
  dat_censored$event <- 0
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, dat_censored)
  expect_error(ertte_rmst(mod, tau = 90), "zero observed events")
})

test_that("ertte_rmst() handles a single-row newdata for both engines", {
  mod_aft <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  mod_cox <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  nd <- ertte_data[1, , drop = FALSE]
  expect_equal(nrow(ertte_rmst(mod_aft, nd, tau = c(30, 60, 90))), 3L)
  expect_equal(nrow(ertte_rmst(mod_cox, nd, tau = c(30, 60, 90))), 3L)
})

test_that("ertte_rmst() on ertte_aft doesn't silently return 0 for very large tau (issue #12)", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  prof <- ertte_data[1, , drop = FALSE]

  # before the fix, stats::integrate()'s adaptive quadrature silently
  # failed (returning 0, with no warning) somewhere between tau = 1e5 and
  # tau = 5e5 for this model/profile -- reparameterising the integral onto
  # the log(t) scale fixes this; fit_rmst should stay close to its
  # converged value (well under the RMST horizon) for tau values that used
  # to trigger the bug.
  rmst <- suppressWarnings(ertte_rmst(mod, prof, tau = c(1e5, 5e5, 1e6, 1e9)))
  expect_true(all(rmst$fit_rmst > 1)) # not silently 0
  # RMST is non-decreasing in tau and must have converged well before
  # these horizons (the survival curve is negligible long before tau=1e5)
  expect_equal(rmst$fit_rmst, rep(rmst$fit_rmst[1], nrow(rmst)), tolerance = 1e-4)
})

test_that("ertte_rmst() on ertte_aft warns for tau far beyond the observed follow-up range", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  max_time <- max(ertte_data$time)
  expect_warning(
    ertte_rmst(mod, ertte_data[1, , drop = FALSE], tau = max_time * 1e5),
    "10000x the last observed"
  )
  expect_no_warning(ertte_rmst(mod, ertte_data[1, , drop = FALSE], tau = 90))
})

test_that("ertte_rmst() recomputes fit_rmst/ci_lower/ci_upper even when newdata already has columns of those names", {
  # regression test for erplots#12: `newdata[rep_rows, ] |> mutate(fit_rmst
  # = fit_rmst, ...)` used to resolve the bare `fit_rmst` (and, in turn,
  # `ci_lower`/`ci_upper`) reference against a pre-existing column of the
  # same name in `newdata` (dplyr's data-mask precedence), rather than the
  # freshly computed local vector -- silently passing through stale values
  # whenever a previous `ertte_rmst()` call's own output was reused as
  # `newdata`. This is exactly what happens if a fitted RMST curve is
  # rebuilt from data produced by an earlier call, or -- the scenario that
  # surfaced the bug -- when `erplots::er_plot()`'s model-curve grid copies
  # every non-exposure column (including a stale `fit_rmst`/`ci_lower`/
  # `ci_upper` triple) from the plot's own data into the prediction grid it
  # hands to `er_predict()`.
  grid <- data.frame(aucss = seq(0, 4700, length.out = 10))

  for (mod in list(
    ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data),
    ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  )) {
    fresh <- ertte_rmst(mod, grid, tau = 90)
    reused <- ertte_rmst(mod, fresh, tau = 90)

    expect_equal(reused$fit_rmst, fresh$fit_rmst)
    expect_equal(reused$ci_lower, fresh$ci_lower)
    expect_equal(reused$ci_upper, fresh$ci_upper)
    # a flat/stale `fit_rmst` (the symptom of the bug) would fail this --
    # RMST strictly decreases with `aucss` for this fitted model
    expect_true(all(diff(reused$fit_rmst) < 0))
  }
})
