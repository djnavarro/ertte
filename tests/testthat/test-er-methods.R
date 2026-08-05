test_that("er_predict.ertte_model() returns the same shape as ertte_landmark()", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  pred <- er_predict.ertte_model(mod, ertte_data[1:5, ], landmark_time = 60)
  ref <- ertte_landmark(mod, ertte_data[1:5, ], landmark_time = 60)
  # `expect_equal(ignore_attr = TRUE)` because row-subsetting a tibble
  # doesn't reliably preserve incidental per-column `label` attributes
  # (from `ertte_data`'s documentation metadata) -- irrelevant here.
  expect_equal(pred, ref, ignore_attr = TRUE)
})

test_that("er_predict.ertte_model() requires a landmark_time or tau argument", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(er_predict.ertte_model(mod, ertte_data[1:5, ]), "landmark_time")
})

test_that("er_predict.ertte_model() forwards tau to ertte_rmst(), renaming fit_rmst to fit_resp", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  pred <- er_predict.ertte_model(mod, ertte_data[1:5, ], tau = 90)
  ref <- ertte_rmst(mod, ertte_data[1:5, ], tau = 90)
  expect_equal(pred$fit_resp, ref$fit_rmst)
  expect_equal(pred$ci_lower, ref$ci_lower)
  expect_equal(pred$ci_upper, ref$ci_upper)
  expect_false("fit_rmst" %in% names(pred))
})

test_that("er_predict.ertte_model() rejects landmark_time and tau together", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(
    er_predict.ertte_model(mod, ertte_data[1:5, ], landmark_time = 90, tau = 90),
    "only one of"
  )
})

test_that("er_predict.ertte_model() rejects a vectorised tau", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(
    er_predict.ertte_model(mod, ertte_data[1:5, ], tau = c(60, 90)),
    "single, strictly positive number"
  )
})

test_that("er_simulate.ertte_model() matches .ertte_simulate_draws()", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  sim <- er_simulate.ertte_model(mod, ertte_data, nsim = 5, seed = 99)
  expect_equal(sim, .ertte_simulate_draws(mod, ertte_data, nsim = 5, seed = 99), ignore_attr = TRUE)
})

test_that("er_simulate.ertte_model() with landmark_time adds fit_resp/sim_resp, matching ertte_landmark()", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  prof <- ertte_data[1, , drop = FALSE]
  sim <- er_simulate.ertte_model(mod, prof, nsim = 3000, seed = 5123, landmark_time = 90)

  expect_true(all(c("fit_resp", "sim_resp", "landmark_time") %in% names(sim)))
  expect_true(all(sim$sim_resp %in% c(0, 1, NA)))
  expect_true(all(sim$fit_resp >= 0 & sim$fit_resp <= 1))

  ref <- ertte_landmark(mod, prof, landmark_time = 90)$fit_resp
  expect_equal(mean(sim$sim_resp, na.rm = TRUE), ref, tolerance = 0.02)
  expect_equal(mean(sim$fit_resp), ref, tolerance = 0.01)
})

test_that("er_simulate.ertte_model() with tau adds fit_resp/sim_resp, matching ertte_rmst()", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  prof <- ertte_data[1, , drop = FALSE]
  sim <- er_simulate.ertte_model(mod, prof, nsim = 3000, seed = 77, tau = 90)

  expect_true(all(c("fit_resp", "sim_resp", "tau") %in% names(sim)))
  expect_true(all(sim$sim_resp >= 0 & sim$sim_resp <= 90 | is.na(sim$sim_resp)))
  expect_true(all(sim$fit_resp >= 0 & sim$fit_resp <= 90))

  ref <- ertte_rmst(mod, prof, tau = 90)$fit_rmst
  expect_equal(mean(sim$sim_resp, na.rm = TRUE), ref, tolerance = 1)
  expect_equal(mean(sim$fit_resp), ref, tolerance = 1)
})

test_that("er_simulate.ertte_model() treats simulated censoring before landmark_time/tau as NA", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  nd <- ertte_data[1:5, ]

  sim_lmk <- er_simulate.ertte_model(mod, nd, nsim = 50, seed = 9, landmark_time = 90, censor_time = 40)
  expect_true(any(is.na(sim_lmk$sim_resp)))
  expect_false(any(sim_lmk$sim_resp %in% 0)) # 40 < 90, so "known non-event by t*" can't occur

  sim_rmst <- er_simulate.ertte_model(mod, nd, nsim = 50, seed = 9, tau = 90, censor_time = 40)
  expect_true(any(is.na(sim_rmst$sim_resp)))
  expect_true(all(sim_rmst$sim_resp[!is.na(sim_rmst$sim_resp)] < 40))
})

test_that("er_simulate.ertte_model() rejects landmark_time and tau together", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(
    er_simulate.ertte_model(mod, ertte_data[1, , drop = FALSE], landmark_time = 90, tau = 90),
    "only one of"
  )
})

test_that("er_simulate.ertte_model() rejects a vectorised tau", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(
    er_simulate.ertte_model(mod, ertte_data[1, , drop = FALSE], tau = c(60, 90)),
    "single, strictly positive number"
  )
})

test_that("er_summary.ertte_model() returns p_value/coefficients/glance", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  smm <- er_summary.ertte_model(mod)
  expect_named(smm, c("p_value", "coefficients", "glance"))
  expect_equal(nrow(smm$coefficients), 2L)
  expect_equal(nrow(smm$glance), 1L)
})

test_that("ertte registers er_predict/er_simulate/er_summary with erplots, if installed", {
  skip_if_not_installed("erplots")
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_true(!is.null(getS3method("er_predict", "ertte_model", envir = asNamespace("erplots"))))
})

test_that("er_predict_survival.ertte_model() matches ertte_predict() for strictly positive time_grid", {
  mod_aft <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  mod_cox <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  nd <- ertte_data[1:5, ]
  tg <- c(30, 60, 90)

  expect_equal(
    er_predict_survival.ertte_model(mod_aft, nd, tg),
    ertte_predict(mod_aft, nd, time = tg)
  )
  expect_equal(
    er_predict_survival.ertte_model(mod_cox, nd, tg),
    ertte_predict(mod_cox, nd, time = tg)
  )
})

test_that("er_predict_survival.ertte_model() treats time_grid = 0 as S(0) = 1, for both engines", {
  # erplots::er_tte_add_model()'s default time_grid spans
  # `object$time$limits`, whose lower end is 0 -- but `ertte_predict()`
  # rejects non-positive `time` outright, so this needs special-case
  # handling rather than a pure pass-through (see R/er-methods.R).
  mod_aft <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  mod_cox <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  nd <- ertte_data[1:3, ]
  tg <- c(0, 30, 60)

  for (mod in list(mod_aft, mod_cox)) {
    out <- er_predict_survival.ertte_model(mod, nd, tg)
    expect_equal(nrow(out), nrow(nd) * length(tg))
    zero_rows <- out[out$time == 0, ]
    expect_true(all(zero_rows$fit_survival == 1))
    expect_true(all(zero_rows$ci_lower == 1))
    expect_true(all(zero_rows$ci_upper == 1))

    # positive entries agree with a direct ertte_predict() call, and
    # row order/values are otherwise unaffected by interleaving the
    # zero-time rows back in
    ref_pos <- ertte_predict(mod, nd, time = tg[tg > 0])
    out_pos <- out[out$time > 0, ]
    expect_equal(out_pos$fit_survival, ref_pos$fit_survival)
    expect_equal(out_pos$ci_lower, ref_pos$ci_lower)
    expect_equal(out_pos$ci_upper, ref_pos$ci_upper)
  }
})

test_that("er_predict_survival.ertte_model() handles a time_grid of all zeros", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  nd <- ertte_data[1:2, ]
  out <- er_predict_survival.ertte_model(mod, nd, time_grid = c(0, 0))
  expect_equal(nrow(out), 4L)
  expect_true(all(out$fit_survival == 1))
})

test_that("er_predict_survival.ertte_model() returns a zero-row tibble for zero-row newdata, both engines", {
  mod_aft <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  mod_cox <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  nd0 <- ertte_data[0, ]

  for (mod in list(mod_aft, mod_cox)) {
    out <- er_predict_survival.ertte_model(mod, nd0, time_grid = c(0, 30))
    expect_equal(nrow(out), 0L)
    expect_true(all(c("time", "fit_survival", "ci_lower", "ci_upper") %in% names(out)))
  }
})

test_that("er_predict_survival.ertte_model() rejects a negative or non-numeric time_grid", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  nd <- ertte_data[1:2, ]
  expect_error(er_predict_survival.ertte_model(mod, nd, time_grid = c(-1, 30)), "non-negative")
  expect_error(er_predict_survival.ertte_model(mod, nd, time_grid = "a"), "non-negative")
})

test_that("ertte registers er_predict_survival with erplots, if installed", {
  skip_if_not_installed("erplots")
  skip_if_not(exists("er_predict_survival", where = asNamespace("erplots")))
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_true(!is.null(getS3method("er_predict_survival", "ertte_model", envir = asNamespace("erplots"))))
})

test_that("er_tte_add_model() renders an S(t) curve/ribbon for an ertte_aft model, via erplots' er_tte() grammar", {
  skip_if_not_installed("erplots")
  skip_if_not(exists("er_tte", where = asNamespace("erplots")))

  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  p <- erplots::er_tte(ertte_data, time, event) |>
    erplots::er_tte_add_curve() |>
    erplots::er_tte_add_model(mod)
  built <- erplots::er_tte_build(p)

  curve <- built$layer$model$config$predictions
  expect_true(all(c("time", "fit_survival", "ci_lower", "ci_upper") %in% names(curve)))
  expect_true(all(curve$fit_survival >= 0 & curve$fit_survival <= 1))
  expect_true(is.unsorted(-curve$fit_survival) == FALSE) # non-increasing in time
})

test_that("er_predict_survival.ertte_model() errors informatively for an all-censored ertte_coxph model, at any positive time_grid", {
  # ertte_predict.ertte_coxph()'s existing .ertte_check_coxph_nevent()
  # guard (a zero-event model has no baseline hazard to build survival
  # predictions from) fires unchanged through this pass-through method,
  # for both a direct call and the full erplots::er_tte_add_model()
  # grammar -- no new gap here, but worth pinning down since this
  # method's `time_grid = 0` special-casing (below) means the guard
  # isn't *always* reached.
  dat_censored <- ertte_data
  dat_censored$event <- 0
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, dat_censored)

  expect_error(
    er_predict_survival.ertte_model(mod, dat_censored[1:2, ], time_grid = c(0, 30)),
    "zero observed events"
  )

  skip_if_not_installed("erplots")
  skip_if_not(exists("er_tte", where = asNamespace("erplots")))
  # `er_tte_add_model()` evaluates the model layer eagerly (unlike some
  # other erplots layers), so the error surfaces here rather than at
  # `er_tte_build()`.
  expect_error(
    {
      p <- erplots::er_tte(dat_censored, time, event) |>
        erplots::er_tte_add_curve() |>
        erplots::er_tte_add_model(mod)
    },
    "zero observed events"
  )
})

test_that("er_predict_survival.ertte_model()'s time_grid = 0 special case bypasses ertte_predict()'s guards entirely", {
  # S(0) = 1 holds by definition regardless of model validity, so a
  # time_grid consisting only of 0 never reaches ertte_predict() (and
  # therefore never triggers .ertte_check_coxph_nevent()) even for an
  # otherwise-unusable all-censored ertte_coxph fit. Documented here as
  # intentional -- not a gap to close -- since it's a trivially true
  # fact, not a model-dependent computation.
  dat_censored <- ertte_data
  dat_censored$event <- 0
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, dat_censored)

  out <- er_predict_survival.ertte_model(mod, dat_censored[1:2, ], time_grid = 0)
  expect_equal(out$fit_survival, c(1, 1))
  expect_equal(out$ci_lower, c(1, 1))
  expect_equal(out$ci_upper, c(1, 1))
})

test_that("er_predict_survival.ertte_model() degrades like ertte_predict() for an all-censored ertte_aft model (no error, but unstable/non-monotonic values)", {
  # Unlike ertte_coxph(), ertte_aft() has no zero-event guard --
  # survreg() itself just warns "did not converge" and returns a fit
  # ertte_predict() still evaluates (see AGENTS.md's "Stress-test
  # findings"). Confirming that carries through this pass-through
  # method unchanged, rather than silently producing something worse.
  dat_censored <- ertte_data
  dat_censored$event <- 0
  mod <- suppressWarnings(ertte_aft(survival::Surv(time, event) ~ aucss, dat_censored))

  out <- suppressWarnings(
    er_predict_survival.ertte_model(mod, dat_censored[1:2, ], time_grid = c(0, 30, 60))
  )
  expect_equal(nrow(out), 6L)
  expect_false(anyNA(out$fit_survival))
  expect_equal(out$fit_survival[out$time == 0], c(1, 1))
})

test_that("er_predict_survival.ertte_model() propagates NA for a single-level-factor ertte_aft model, but returns finite values for ertte_coxph", {
  # Matches ertte_predict()'s own documented per-engine difference
  # (AGENTS.md's "Stress-test findings"): predict.survreg() propagates
  # an aliased factor level's NA coefficient straight through, while
  # survival::survfit() silently drops the aliased column. `time = 0`
  # rows are unaffected either way (S(0) = 1 needs no coefficient).
  dat_single <- ertte_data[ertte_data$sex == levels(ertte_data$sex)[1], ]
  mod_aft <- ertte_aft(survival::Surv(time, event) ~ sex + aucss, dat_single)
  mod_cox <- ertte_coxph(survival::Surv(time, event) ~ sex + aucss, dat_single)
  expect_true(anyNA(coef(mod_aft)))
  expect_true(anyNA(coef(mod_cox)))

  nd <- dat_single[1:2, ]
  out_aft <- er_predict_survival.ertte_model(mod_aft, nd, time_grid = c(0, 30))
  out_cox <- er_predict_survival.ertte_model(mod_cox, nd, time_grid = c(0, 30))

  expect_equal(out_aft$fit_survival[out_aft$time == 0], c(1, 1))
  expect_true(all(is.na(out_aft$fit_survival[out_aft$time > 0])))

  expect_equal(out_cox$fit_survival[out_cox$time == 0], c(1, 1))
  expect_false(anyNA(out_cox$fit_survival[out_cox$time > 0]))
})

test_that("er_tte_add_model() renders one curve per stratum for an ertte_coxph model", {
  skip_if_not_installed("erplots")
  skip_if_not(exists("er_tte", where = asNamespace("erplots")))

  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss + sex, ertte_data)
  p <- erplots::er_tte(ertte_data, time, event, stratify_by = sex) |>
    erplots::er_tte_add_curve() |>
    erplots::er_tte_add_model(mod)
  built <- erplots::er_tte_build(p)

  curve <- built$layer$model$config$predictions
  expect_true("sex" %in% names(curve))
  expect_equal(sort(unique(as.character(curve$sex))), sort(levels(ertte_data$sex)))
})

test_that("landmark_time reaches er_predict.ertte_model() through erplots' er_plot_add_model(predict_args = ...)", {
  # erplots#10/#11: er_plot_add_model()'s `...` used to reach only its
  # style builder, never `er_predict()` -- so `landmark_time` (which
  # has no other slot in `er_predict()`'s fixed signature) couldn't
  # reach `ertte_landmark()` this way. Fixed upstream via a dedicated
  # `predict_args` argument; skip gracefully on an erplots install that
  # predates it, rather than erroring on an unrecognised argument.
  skip_if_not_installed("erplots")
  skip_if_not("predict_args" %in% names(formals(erplots::er_plot_add_model)))

  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  dat_landmark <- ertte_data |>
    dplyr::mutate(
      event_by_90 = dplyr::case_when(
        event == 1 & time <= 90 ~ 1,
        time > 90 ~ 0,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::filter(!is.na(event_by_90))

  p <- erplots::er_plot(dat_landmark, aucss, event_by_90) |>
    erplots::er_plot_add_model(mod, predict_args = list(landmark_time = 90))
  built <- erplots::er_plot_build(p)
  expect_s3_class(built$output, "patchwork")
})

test_that("er_plot()'s landmark model curve is correct when built from a reused ertte_landmark() output (erplots#12)", {
  # analogous to the erplots#12 dplyr data-mask shadowing bug fixed in
  # ertte_rmst() -- see the matching RMST test below. ertte_landmark()'s
  # mutate() never bare-references its own output column names on the
  # RHS, so this isn't expected to reproduce the bug -- kept for test
  # parity with the RMST coverage and to check curve *values*, not just
  # `built$output`'s class (the check that let the RMST bug slip through
  # review originally).
  skip_if_not_installed("erplots")
  skip_if_not("predict_args" %in% names(formals(erplots::er_plot_add_model)))

  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  lmk_ref <- ertte_landmark(mod, ertte_data, landmark_time = 90)

  p <- erplots::er_plot(lmk_ref, aucss, fit_resp) |>
    erplots::er_plot_add_model(mod, predict_args = list(landmark_time = 90))
  built <- erplots::er_plot_build(p)

  curve <- built$layer$model$config$predictions
  # non-decreasing, not strictly increasing: fit_resp saturates to (a
  # floating-point) 1 at the high end of the exposure grid, producing
  # ties there rather than a decrease
  expect_false(is.unsorted(curve$fit_resp))
  expect_gt(diff(range(curve$fit_resp)), 0.3)
})

test_that("tau reaches er_predict.ertte_model() through erplots' er_plot_add_model(predict_args = ...)", {
  skip_if_not_installed("erplots")
  skip_if_not("predict_args" %in% names(formals(erplots::er_plot_add_model)))

  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  rmst_ref <- ertte_rmst(mod, ertte_data, tau = 90)

  p <- erplots::er_plot(rmst_ref, aucss, fit_rmst) |>
    erplots::er_plot_add_model(mod, predict_args = list(tau = 90))
  built <- erplots::er_plot_build(p)
  expect_s3_class(built$output, "patchwork")
})

test_that("er_plot()'s RMST model curve is correct when built from a reused ertte_rmst() output (erplots#12)", {
  # `built$output`'s class alone (checked above) doesn't catch a curve
  # that renders but is numerically wrong. erplots' `.get_model_predictions()`
  # fills its exposure grid with reference values for every other column of
  # the plot's own data -- including, here, `rmst_ref`'s own `fit_rmst`/
  # `ci_lower`/`ci_upper` columns -- which used to trigger a dplyr
  # data-mask shadowing bug in `ertte_rmst()` (fixed in R/ertte-rmst.R)
  # that made the model curve flat instead of varying with exposure.
  skip_if_not_installed("erplots")
  skip_if_not("predict_args" %in% names(formals(erplots::er_plot_add_model)))

  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  rmst_ref <- ertte_rmst(mod, ertte_data, tau = 90)

  p <- erplots::er_plot(rmst_ref, aucss, fit_rmst) |>
    erplots::er_plot_add_model(mod, predict_args = list(tau = 90))
  built <- erplots::er_plot_build(p)

  curve <- built$layer$model$config$predictions
  expect_true(all(diff(curve$fit_resp) < 0))
  expect_gt(diff(range(curve$fit_resp)), 30)
})

test_that("landmark_time reaches a landmark-binary VPC through erplots' er_vpc_add_simulated(simulate_args = ...)", {
  skip_if_not_installed("erplots")
  skip_if_not("simulate_args" %in% names(formals(erplots::er_vpc_add_simulated)))

  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  dat_landmark <- ertte_data |>
    dplyr::mutate(
      event_by_90 = dplyr::case_when(
        event == 1 & time <= 90 ~ 1,
        time > 90 ~ 0,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::filter(!is.na(event_by_90))

  p <- erplots::er_vpc(dat_landmark, aucss, event_by_90, plot_by = aucss) |>
    erplots::er_vpc_add_observed() |>
    erplots::er_vpc_add_simulated(
      model = mod, seed = 42, nsim = 100,
      simulate_args = list(landmark_time = 90)
    )
  built <- erplots::er_vpc_build(p)
  expect_s3_class(built$output, "ggplot")
})

test_that("tau reaches an RMST VPC through erplots' er_vpc_add_simulated(simulate_args = ...)", {
  skip_if_not_installed("erplots")
  skip_if_not("simulate_args" %in% names(formals(erplots::er_vpc_add_simulated)))

  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  dat_rmst <- ertte_data |>
    dplyr::mutate(
      rmst_obs = dplyr::case_when(
        event == 1 ~ pmin(time, 90),
        time >= 90 ~ 90,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::filter(!is.na(rmst_obs))

  p <- erplots::er_vpc(dat_rmst, aucss, rmst_obs, plot_by = aucss) |>
    erplots::er_vpc_add_observed() |>
    erplots::er_vpc_add_simulated(
      model = mod, seed = 42, nsim = 100,
      simulate_args = list(tau = 90)
    )
  built <- erplots::er_vpc_build(p)
  expect_s3_class(built$output, "ggplot")
})
