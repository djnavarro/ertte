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
