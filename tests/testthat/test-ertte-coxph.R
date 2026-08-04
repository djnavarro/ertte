test_that("ertte_coxph() reproduces direct coxph() results", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  ref <- survival::coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_equal(unname(coef(mod)), unname(coef(ref)))
  expect_equal(as.numeric(stats::logLik(mod)), as.numeric(stats::logLik(ref)))
  expect_s3_class(mod, "ertte_coxph")
  expect_s3_class(mod, "ertte_model")
  expect_s3_class(mod, "coxph")
  expect_identical(mod$ertte$type, "coxph")
})

test_that("ertte_coxph() stores the fitting data on the object", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_identical(mod$data, ertte_data)
})

test_that("ertte_fun() on ertte_coxph reproduces the fitted model's own predictions", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss + sex, ertte_data)
  mod_fun <- ertte_fun(mod)
  s1 <- unname(mod_fun(time = 60))
  s2 <- ertte_predict(mod, ertte_data, time = 60)$fit_survival
  expect_equal(s1, s2)
})

test_that("ertte_fun() on ertte_coxph works with custom newdata and vectorised time", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss + sex, ertte_data)
  mod_fun <- ertte_fun(mod)
  nd <- ertte_data[1:5, ]
  times <- c(30, 60, 90, 30, 60)
  s1 <- unname(mod_fun(data = nd, time = times))
  s2 <- vapply(seq_along(times), function(i) {
    ertte_predict(mod, nd[i, ], time = times[i])$fit_survival
  }, numeric(1))
  expect_equal(s1, s2)
})

test_that("ertte_fun() on ertte_coxph responds to a custom param", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  mod_fun <- ertte_fun(mod)
  s1 <- mod_fun(time = 60)
  par2 <- coef(mod)
  par2["aucss"] <- par2["aucss"] * 2
  s2 <- mod_fun(param = par2, time = 60)
  expect_false(isTRUE(all.equal(unname(s1), unname(s2))))
})

test_that("ertte_fun() on ertte_coxph checks param length", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  mod_fun <- ertte_fun(mod)
  expect_error(mod_fun(param = c(1, 2, 3), time = 60), "length")
})

test_that("ertte_coxph() fits with model = TRUE (required by survfit())", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_false(is.null(mod$model))
})

test_that("ertte_predict() on ertte_coxph returns a tidy tibble with sane survival probabilities", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  pred <- ertte_predict(mod, ertte_data[1:5, ], time = c(30, 60, 90))
  expect_s3_class(pred, "tbl_df")
  expect_equal(nrow(pred), 15L)
  expect_true(all(pred$fit_survival >= 0 & pred$fit_survival <= 1))
  expect_true(all(pred$ci_lower <= pred$fit_survival))
  expect_true(all(pred$ci_upper >= pred$fit_survival))
  # survival probability should decrease with time, for a fixed subject
  one_subject <- pred[pred$id == pred$id[1], ]
  expect_true(all(diff(one_subject$fit_survival) <= 0))
})

test_that("ertte_predict() on ertte_coxph handles a single-row newdata", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  pred <- ertte_predict(mod, ertte_data[1, , drop = FALSE], time = c(30, 60, 90))
  expect_equal(nrow(pred), 3L)
  expect_true(all(diff(pred$fit_survival) <= 0))
})

test_that("ertte_predict() on ertte_coxph extrapolates times beyond the observed range", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  pred <- ertte_predict(mod, ertte_data[1, , drop = FALSE], time = c(180, 10000))
  # extrapolated survival should be held constant (flat step-function
  # extension), not zero/NA
  expect_equal(pred$fit_survival[1], pred$fit_survival[2])
})

test_that("ertte_predict() on ertte_coxph validates conf_level and time", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(ertte_predict(mod, ertte_data, time = 60, conf_level = 1.5), "`conf_level`")
  expect_error(ertte_predict(mod, ertte_data, time = -1), "`time`")
  expect_error(ertte_predict(mod, ertte_data, time = "a"), "`time`")
})

test_that("ertte_add_term()/ertte_remove_term() work for ertte_coxph models", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  mod2 <- ertte_add_term(mod, ~sex)
  expect_s3_class(mod2, "ertte_coxph")
  expect_true("sex" %in% attr(stats::terms(mod2), "term.labels"))

  mod3 <- ertte_remove_term(mod2, ~sex)
  expect_s3_class(mod3, "ertte_coxph")
  expect_false("sex" %in% attr(stats::terms(mod3), "term.labels"))
})

test_that("simulate() on ertte_coxph has the expected shape", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  sim <- simulate(mod, nsim = 5, seed = 321)
  expect_s3_class(sim, "tbl_df")
  expect_equal(nrow(sim), nrow(ertte_data) * 5L)
  expect_true(all(c("dat_id", "sim_id", "sim_time", "sim_event", "coef_aucss") %in% names(sim)))
  expect_true(all(sim$sim_event %in% c(0, 1)))
  # by default (no `censor_time`), only *censored* rows' simulated draws
  # are capped at their observed exit time (see test-ertte-simulate.R's
  # equivalent comment for the AFT engine)
  is_censored <- ertte_data$event[sim$dat_id] == 0
  expect_true(all(sim$sim_time[is_censored] <= ertte_data$time[sim$dat_id][is_censored] + 1e-8))
})

test_that("simulate() on ertte_coxph is reproducible given a seed", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  sim_a <- simulate(mod, nsim = 5, seed = 42)
  sim_b <- simulate(mod, nsim = 5, seed = 42)
  expect_equal(sim_a, sim_b)
})

test_that("simulate() on ertte_coxph reports an auto-picked seed", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_message(simulate(mod, nsim = 5), "Using seed")
})

test_that(".ertte_simulate_draws() on ertte_coxph requires response columns in newdata", {
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  expect_error(
    .ertte_simulate_draws(mod, newdata = ertte_data[c("aucss")], nsim = 5, seed = 1),
    "must contain"
  )
})

test_that("simulate() on ertte_coxph reproduces the fitted baseline survival curve", {
  # simulating many draws at the reference (mean) covariate profile
  # should reproduce the fitted S0(t) closely -- a sanity check on the
  # baseline-hazard inversion, not a tight statistical test. Uses a
  # row with long observed follow-up (administratively censored at 180
  # days) so the sim_time <- pmin(sim_time_raw, obs_time) censoring cap
  # doesn't interfere with the comparison at time = 90.
  mod <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  ref_row <- which(ertte_data$time >= 170)[1]
  nd <- ertte_data[rep(ref_row, 10000), ]
  nd$aucss <- mean(ertte_data$aucss)
  sim <- simulate(mod, nsim = 1, newdata = nd, seed = 777)
  emp_surv_90 <- mean(sim$sim_time > 90)
  theo_surv_90 <- ertte_predict(mod, nd[1, ], time = 90)$fit_survival
  expect_equal(emp_surv_90, theo_surv_90, tolerance = 0.05)
})

test_that("ertte_scm_forward()/ertte_scm_backward() work for ertte_coxph models", {
  mod0 <- ertte_coxph(survival::Surv(time, event) ~ aucss, ertte_data)
  mod1 <- ertte_scm_forward(mod0, candidates = c("sex", "dose", "age", "weight"), seed = 4821)
  expect_s3_class(mod1, "ertte_coxph")
  # ertte_data's ground truth includes a sex effect
  expect_true("sex" %in% attr(stats::terms(mod1), "term.labels"))

  mod2 <- ertte_coxph(survival::Surv(time, event) ~ aucss + sex + dose, ertte_data)
  mod3 <- ertte_scm_backward(mod2, candidates = c("sex", "dose"), seed = 912)
  labs <- attr(stats::terms(mod3), "term.labels")
  expect_true("aucss" %in% labs)
  expect_false("dose" %in% labs) # not part of the ground truth model
})
