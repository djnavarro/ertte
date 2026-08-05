test_that("ertte_landmark() is 1 - ertte_predict() with CI bounds swapped (AFT)", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  nd <- ertte_data[1:5, ]
  lmk <- ertte_landmark(mod, nd, landmark_time = 90)
  ref <- ertte_predict(mod, nd, time = 90)

  expect_equal(lmk$fit_resp, 1 - ref$fit_survival)
  expect_equal(lmk$ci_lower, 1 - ref$ci_upper)
  expect_equal(lmk$ci_upper, 1 - ref$ci_lower)
  expect_equal(lmk$landmark_time, ref$time)
})

test_that("ertte_landmark() is 1 - ertte_predict() with CI bounds swapped (Cox PH)", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  nd <- ertte_data[1:5, ]
  lmk <- ertte_landmark(mod, nd, landmark_time = 90)
  ref <- ertte_predict(mod, nd, time = 90)

  expect_equal(lmk$fit_resp, 1 - ref$fit_survival)
  expect_equal(lmk$ci_lower, 1 - ref$ci_upper)
  expect_equal(lmk$ci_upper, 1 - ref$ci_lower)
  expect_equal(lmk$landmark_time, ref$time)
})

test_that("ertte_landmark()'s interval brackets its point estimate", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  lmk <- ertte_landmark(mod, ertte_data[1:20, ], landmark_time = 120)
  expect_true(all(lmk$ci_lower <= lmk$fit_resp))
  expect_true(all(lmk$fit_resp <= lmk$ci_upper))
})

test_that("ertte_landmark()'s output has the expected columns, no leftovers", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  nd <- ertte_data[1:5, ]
  lmk <- ertte_landmark(mod, nd, landmark_time = 60)
  expect_equal(
    setdiff(names(lmk), names(nd)),
    c("landmark_time", "fit_resp", "ci_lower", "ci_upper")
  )
  expect_false(any(c("time", "fit_survival") %in% names(lmk)))
})

test_that("ertte_landmark() validates landmark_time", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  nd <- ertte_data[1:5, ]
  expect_error(ertte_landmark(mod, nd, landmark_time = c(30, 60)), "landmark_time")
  expect_error(ertte_landmark(mod, nd, landmark_time = -10), "landmark_time")
  expect_error(ertte_landmark(mod, nd, landmark_time = 0), "landmark_time")
  expect_error(ertte_landmark(mod, nd, landmark_time = NA_real_), "landmark_time")
  expect_error(ertte_landmark(mod, nd, landmark_time = "90"), "landmark_time")
})

test_that("ertte_landmark() defaults newdata to the fitted data", {
  mod <- ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data)
  lmk <- ertte_landmark(mod, landmark_time = 90)
  expect_equal(nrow(lmk), nrow(ertte_data))
})

test_that("ertte_landmark() inherits ertte_predict()'s all-censored-Cox guard", {
  dat_censored <- ertte_data[1:60, ]
  dat_censored$event <- 0
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, dat_censored)
  expect_error(ertte_landmark(mod, landmark_time = 90), "zero observed events")
})

test_that("ertte_landmark() recomputes landmark_time/fit_resp/ci_lower/ci_upper even when newdata already has columns of those names", {
  # analogous to the erplots#12 dplyr data-mask shadowing bug fixed in
  # ertte_rmst() (`fit_rmst = fit_rmst` resolved against a pre-existing
  # `newdata` column instead of the freshly computed local vector -- see
  # test-ertte-rmst.R). ertte_landmark()'s own mutate() never bare-
  # references its own output column names on the RHS (`landmark_time =
  # time`, `fit_resp = 1 - fit_survival`, etc.), so this isn't expected to
  # reproduce the bug -- kept as a regression test for parity with the
  # RMST coverage, and as a guard against a future refactor reintroducing
  # a self-referential assignment.
  grid <- data.frame(aucss = seq(0, 4700, length.out = 10))

  for (mod in list(
    ertte_aft(survival::Surv(time, event) ~ aucss, ertte_data),
    ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  )) {
    fresh <- ertte_landmark(mod, grid, landmark_time = 90)
    reused <- ertte_landmark(mod, fresh, landmark_time = 90)

    expect_equal(reused$fit_resp, fresh$fit_resp)
    expect_equal(reused$ci_lower, fresh$ci_lower)
    expect_equal(reused$ci_upper, fresh$ci_upper)
    expect_equal(reused$landmark_time, fresh$landmark_time)
    # a flat/stale `fit_resp` (the symptom the RMST bug produced) would
    # fail this -- P(event by t*) strictly increases with `aucss` here
    expect_true(all(diff(reused$fit_resp) > 0))
  }
})
