# Tests for the refined administrative-censoring simulation behaviour
# (see AGENTS.md's "Administrative-censoring simulation" bullet) and the
# `censor_time` argument to `simulate()`.

test_that(".ertte_check_censor_time() validates length and positivity", {
  expect_null(.ertte_check_censor_time(NULL, 10))
  expect_equal(.ertte_check_censor_time(5, 3), c(5, 5, 5))
  expect_equal(.ertte_check_censor_time(c(1, 2, 3), 3), c(1, 2, 3))
  expect_error(.ertte_check_censor_time(c(1, 2), 3), "length 1")
  expect_error(.ertte_check_censor_time(-1, 3), "strictly positive")
  expect_error(.ertte_check_censor_time(NA, 3), "strictly positive")
  expect_error(.ertte_check_censor_time("a", 3), "strictly positive")
})

test_that("default simulation leaves event rows uncensored (ertte_aft)", {
  mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
  sim <- simulate(mod, nsim = 30, seed = 501)
  is_event <- ertte_data$event[sim$dat_id] == 1
  is_censored <- !is_event
  # censored rows: unchanged, exact cap
  expect_true(all(sim$sim_time[is_censored] <= ertte_data$time[sim$dat_id][is_censored] + 1e-8))
  # event rows: at least some simulated draws should exceed the row's
  # own observed event day, since they're no longer capped there
  expect_true(any(sim$sim_time[is_event] > ertte_data$time[sim$dat_id][is_event] + 1e-8))
  # every event-row draw is marked as an event (uncensored)
  expect_true(all(sim$sim_event[is_event] == 1))
})

test_that("default simulation leaves event rows uncensored (ertte_coxph)", {
  mod <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
  sim <- simulate(mod, nsim = 30, seed = 502)
  is_event <- ertte_data$event[sim$dat_id] == 1
  is_censored <- !is_event
  expect_true(all(sim$sim_time[is_censored] <= ertte_data$time[sim$dat_id][is_censored] + 1e-8))
  # some (not necessarily all, since the baseline hazard's support is
  # finite -- see the Inf-handling note below) event-row draws should
  # exceed the observed event day
  expect_true(any(sim$sim_time[is_event] > ertte_data$time[sim$dat_id][is_event] + 1e-8))
})

test_that("censor_time overrides the event/censored default split uniformly", {
  mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
  sim <- simulate(mod, nsim = 20, seed = 503, censor_time = 90)
  # every row (event or censored) is now capped at 90, regardless of its
  # own observed event status
  expect_true(all(sim$sim_time <= 90 + 1e-8))
  expect_true(all(sim$sim_event == as.numeric(sim$sim_time < 90 - 1e-8) | sim$sim_time == 90))
})

test_that("censor_time accepts a per-row vector matching newdata's admin_censor column", {
  mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
  sim <- simulate(mod, nsim = 10, seed = 504, censor_time = ertte_data$admin_censor)
  expect_true(all(sim$sim_time <= ertte_data$admin_censor[sim$dat_id] + 1e-8))
  # this is (approximately) how `ertte_data` was actually generated --
  # simulated censoring at exactly the 180-day cutoff should happen for a
  # non-trivial share of draws, consistent with the fitted survival
  # curve at t = 180 being well below 1
  prop_censored_at_cutoff <- mean(sim$sim_time >= 180 - 1e-8 & sim$sim_event == 0)
  expect_true(prop_censored_at_cutoff > 0)
})

test_that("censor_time errors informatively on the wrong length", {
  mod <- ertte_aft(Surv(time, event) ~ aucss, ertte_data)
  expect_error(simulate(mod, nsim = 5, seed = 505, censor_time = c(90, 100)), "length 1")
})

test_that("coxph simulation treats draws beyond the baseline hazard's support as censored, not events", {
  mod <- ertte_coxph(Surv(time, event) ~ aucss, ertte_data)
  # a large nsim makes it near-certain some draws exceed the baseline
  # hazard's support (`sim_time_raw == Inf` before capping)
  sim <- simulate(mod, nsim = 200, seed = 506)
  bh_max <- max(survival::basehaz(mod, centered = TRUE)$time)
  extrapolated <- sim$sim_time >= bh_max - 1e-8 & ertte_data$event[sim$dat_id] == 1
  # any row whose simulated time reaches the hazard's support boundary
  # while its observed row was an event must be marked non-event -- an
  # extrapolation artifact, not a real simulated event
  expect_true(all(sim$sim_event[extrapolated] == 0))
})
