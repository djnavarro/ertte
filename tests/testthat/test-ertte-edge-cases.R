# Edge-case coverage from the design issue's test-suite wishlist (see
# AGENTS.md): all-censored data, single-level ("single stratum")
# categorical covariates, and heavy ties in event times. The
# zero-exposure/placebo group is already exercised throughout the rest
# of the test suite via `ertte_data`'s `dose == 0`/`aucss == 0` rows, so
# isn't repeated here.

# --- all-censored data ------------------------------------------------

test_that("ertte_aft() fits all-censored data (with a non-convergence warning) and predicts", {
  dat <- ertte_data[1:60, ]
  dat$event <- 0
  expect_warning(
    mod <- ertte_aft(Surv(time, event) ~ aucss, dat),
    "did not converge"
  )
  pred <- ertte_predict(mod, dat[1:3, ], time = c(30, 60))
  expect_equal(nrow(pred), 6L)
  expect_true(all(pred$fit_survival >= 0 & pred$fit_survival <= 1))
})

test_that("ertte_coxph() fits all-censored data with NA coefficients (0 events)", {
  dat <- ertte_data[1:60, ]
  dat$event <- 0
  mod <- ertte_coxph(Surv(time, event) ~ aucss, dat)
  expect_equal(mod$nevent, 0)
  expect_true(is.na(stats::coef(mod)["aucss"]))
})

test_that("ertte_predict()/ertte_fun()/simulate() error informatively for a zero-event ertte_coxph model", {
  dat <- ertte_data[1:60, ]
  dat$event <- 0
  mod <- ertte_coxph(Surv(time, event) ~ aucss, dat)

  expect_error(
    ertte_predict(mod, dat[1:3, ], time = c(30, 60)),
    "zero observed events"
  )
  expect_error(ertte_fun(mod), "zero observed events")
  expect_error(
    simulate(mod, newdata = dat[1:3, ], nsim = 2, seed = 5),
    "zero observed events"
  )
})

# --- single-level ("single stratum") categorical covariate ------------

test_that("ertte_aft()/ertte_coxph() fit a single-level factor with an NA (aliased) coefficient", {
  dat <- ertte_data[ertte_data$sex == "Female", ]
  expect_equal(length(unique(dat$sex)), 1L)

  mod_aft <- ertte_aft(Surv(time, event) ~ aucss + sex, dat)
  expect_true(is.na(stats::coef(mod_aft)["sexMale"]))

  mod_cox <- ertte_coxph(Surv(time, event) ~ aucss + sex, dat)
  expect_true(is.na(stats::coef(mod_cox)["sexMale"]))
})

test_that("ertte_predict.ertte_aft() propagates NA for an aliased single-level factor term", {
  dat <- ertte_data[ertte_data$sex == "Female", ]
  mod <- ertte_aft(Surv(time, event) ~ aucss + sex, dat)
  pred <- ertte_predict(mod, dat[1:3, ], time = c(30, 60))
  # `predict.survreg()` propagates the aliased coefficient's NA straight
  # through to the linear predictor, and therefore to every downstream
  # survival probability/CI -- a real (if degenerate) consequence of the
  # singularity, not a bug ertte introduces.
  expect_true(all(is.na(pred$fit_survival)))
})

test_that("ertte_predict.ertte_coxph() still returns finite predictions for an aliased single-level factor term", {
  dat <- ertte_data[ertte_data$sex == "Female", ]
  mod <- ertte_coxph(Surv(time, event) ~ aucss + sex, dat)
  # `survival::survfit()` silently drops the aliased column rather than
  # propagating its NA -- genuinely different from the AFT engine's
  # behaviour above (and from `ertte_fun.ertte_coxph()`'s manual
  # `mm %*% param`, which *does* propagate the NA), but not a crash.
  pred <- ertte_predict(mod, dat[1:3, ], time = c(30, 60))
  expect_true(all(!is.na(pred$fit_survival)))
  expect_true(all(pred$fit_survival >= 0 & pred$fit_survival <= 1))

  mod_fun <- ertte_fun(mod)
  expect_true(all(is.na(mod_fun(time = 60))))
})

test_that("ertte_add_term()/SCM don't select an aliased single-level candidate", {
  dat <- ertte_data[ertte_data$sex == "Female", ]
  mod0 <- ertte_aft(Surv(time, event) ~ aucss, dat)
  # Adding `sex` here doesn't error (`anova()` reports a trivial p = 1,
  # 0 deviance, rather than an NA p-value) since there's no variation
  # left in a Female-only subset -- it just never looks significant, so
  # SCM correctly never selects it.
  mod1 <- ertte_scm_forward(mod0, candidates = c("sex", "dose"), seed = 77)
  expect_false("sex" %in% attr(stats::terms(mod1), "term.labels"))
})

# --- heavy ties in event times ----------------------------------------

test_that("ertte_coxph() fits normally under heavy ties in event times", {
  dat <- ertte_data[1:100, ]
  dat$time[dat$event == 1] <- 50 # tie every event at the same time
  mod <- ertte_coxph(Surv(time, event) ~ aucss, dat)
  expect_equal(mod$method, "efron")
  expect_false(anyNA(stats::coef(mod)))

  pred <- ertte_predict(mod, dat[1:5, ], time = c(30, 60))
  expect_true(all(pred$fit_survival >= 0 & pred$fit_survival <= 1))
  # survival at the tied event time should differ from just before it
  # for at least one subject (the drop shouldn't be smoothed away)
  pred2 <- ertte_predict(mod, dat[1:5, ], time = c(49, 51))
  expect_false(isTRUE(all.equal(pred2$fit_survival[pred2$time == 49], pred2$fit_survival[pred2$time == 51])))
})

test_that("simulate() reproduces sensible event times under heavy ties", {
  dat <- ertte_data[1:100, ]
  dat$time[dat$event == 1] <- 50
  mod <- ertte_coxph(Surv(time, event) ~ aucss, dat)
  sim <- simulate(mod, newdata = dat[1:10, ], nsim = 5, seed = 33)
  expect_equal(nrow(sim), 50L)
  # only censored rows are capped at their observed exit time by default
  ref <- dat[1:10, ]
  is_censored <- ref$event[sim$dat_id] == 0
  expect_true(all(sim$sim_time[is_censored] <= ref$time[sim$dat_id][is_censored] + 1e-8))
})

test_that("ertte_aft() fits normally under heavy ties in event times", {
  dat <- ertte_data[1:100, ]
  dat$time[dat$event == 1] <- 50
  mod <- ertte_aft(Surv(time, event) ~ aucss, dat)
  expect_false(anyNA(stats::coef(mod)))
  pred <- ertte_predict(mod, dat[1:5, ], time = c(30, 60))
  expect_true(all(pred$fit_survival >= 0 & pred$fit_survival <= 1))
})

# --- zero-row newdata (issue #10) --------------------------------------

test_that("ertte_predict()/ertte_rmst()/ertte_landmark() return a zero-row tibble for empty newdata, both engines", {
  mod_aft <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
  mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
  nd0 <- ertte_data[0, ]

  pred_aft <- ertte_predict(mod_aft, nd0, time = c(30, 60))
  pred_cox <- ertte_predict(mod_cox, nd0, time = c(30, 60))
  expect_equal(nrow(pred_aft), 0L)
  expect_equal(nrow(pred_cox), 0L)
  expect_true(all(c("time", "fit_survival", "ci_lower", "ci_upper") %in% names(pred_aft)))
  expect_true(all(c("time", "fit_survival", "ci_lower", "ci_upper") %in% names(pred_cox)))

  rmst_aft <- ertte_rmst(mod_aft, nd0, tau = c(30, 60))
  rmst_cox <- ertte_rmst(mod_cox, nd0, tau = c(30, 60))
  expect_equal(nrow(rmst_aft), 0L)
  expect_equal(nrow(rmst_cox), 0L)
  expect_true(all(c("tau", "fit_rmst", "ci_lower", "ci_upper") %in% names(rmst_aft)))
  expect_true(all(c("tau", "fit_rmst", "ci_lower", "ci_upper") %in% names(rmst_cox)))

  lmk_aft <- ertte_landmark(mod_aft, nd0, landmark_time = 30)
  lmk_cox <- ertte_landmark(mod_cox, nd0, landmark_time = 30)
  expect_equal(nrow(lmk_aft), 0L)
  expect_equal(nrow(lmk_cox), 0L)
})

# --- conf_level = 0/1 boundaries, ertte_coxph (issue #11) ---------------

test_that("ertte_predict.ertte_coxph() handles conf_level = 0/1 like ertte_predict.ertte_aft() does", {
  mod_aft <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
  mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
  nd3 <- ertte_data[1:3, ]

  pred_cox0 <- ertte_predict(mod_cox, nd3, time = 30, conf_level = 0)
  expect_equal(pred_cox0$ci_lower, pred_cox0$fit_survival)
  expect_equal(pred_cox0$ci_upper, pred_cox0$fit_survival)

  pred_cox1 <- ertte_predict(mod_cox, nd3, time = 30, conf_level = 1)
  expect_true(all(pred_cox1$ci_lower == 0))
  expect_true(all(pred_cox1$ci_upper == 1))

  # matches ertte_predict.ertte_aft()'s own boundary behaviour (its
  # CDF back-transform naturally collapses/bounds the interval there)
  pred_aft0 <- ertte_predict(mod_aft, nd3, time = 30, conf_level = 0)
  pred_aft1 <- ertte_predict(mod_aft, nd3, time = 30, conf_level = 1)
  expect_equal(pred_aft0$ci_lower, pred_aft0$fit_survival)
  expect_true(all(pred_aft1$ci_lower == 0))
  expect_true(all(pred_aft1$ci_upper == 1))
})

test_that("ertte_rmst.ertte_coxph() handles conf_level = 0/1 without erroring, matching the AFT engine's pattern", {
  mod_aft <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
  mod_cox <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
  nd3 <- ertte_data[1:3, ]

  rmst_cox0 <- ertte_rmst(mod_cox, nd3, tau = 60, conf_level = 0)
  expect_equal(rmst_cox0$ci_lower, rmst_cox0$fit_rmst)
  expect_equal(rmst_cox0$ci_upper, rmst_cox0$fit_rmst)

  rmst_cox1 <- ertte_rmst(mod_cox, nd3, tau = 60, conf_level = 1)
  expect_true(all(is.infinite(rmst_cox1$ci_lower)))
  expect_true(all(is.infinite(rmst_cox1$ci_upper)))

  rmst_aft0 <- ertte_rmst(mod_aft, nd3, tau = 60, conf_level = 0)
  rmst_aft1 <- ertte_rmst(mod_aft, nd3, tau = 60, conf_level = 1)
  expect_equal(rmst_aft0$ci_lower, rmst_aft0$fit_rmst)
  expect_true(all(is.infinite(rmst_aft1$ci_lower)))
  expect_true(all(is.infinite(rmst_aft1$ci_upper)))

  # ordinary confidence levels are unaffected by the conf.int-placeholder
  # refactor in ertte_rmst.ertte_coxph()
  rmst_cox_90 <- ertte_rmst(mod_cox, nd3, tau = 60, conf_level = 0.9)
  expect_true(all(rmst_cox_90$ci_lower < rmst_cox_90$fit_rmst))
  expect_true(all(rmst_cox_90$ci_upper > rmst_cox_90$fit_rmst))
})
