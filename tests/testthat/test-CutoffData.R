make_toy <- function() {
  data.frame(sim = c(1, 1, 2), group = "A", accrual_time = c(0, 5, 1),
             pfs_time = c(3, 2, 10), os_time = c(8, 2, 20),
             progression = c(1, 0, 1), response = c(1, 0, 1), ttr = c(2, NA, 4),
             dropout_time = c(Inf, Inf, 6))
}

test_that("censoring at a cutoff per simulated trial", {
  out <- CutoffData(make_toy(), cutoff = c(6, 4))
  expect_equal(nrow(out), 3)
  expect_equal(out$pfs_tte, c(3, 1, 3))
  expect_equal(out$pfs_event, c(1L, 0L, 0L))
  expect_equal(out$os_tte, c(6, 1, 3))
  expect_equal(out$os_event, c(0L, 0L, 0L))
  expect_equal(out$response_obs, c(1L, 0L, 0L))
  expect_equal(out$cutoff, c(6, 6, 4))
  expect_equal(out$os_calendar_time, c(6, 6, 4))
})

test_that("patients enrolled after the cutoff are dropped", {
  out <- CutoffData(make_toy(), cutoff = 0.5)
  expect_equal(nrow(out), 1)
  expect_equal(out$pfs_tte, 0.5)
  expect_equal(out$response_obs, 0L)
})

test_that("no cutoff keeps all follow-up and dropout still applies", {
  out <- CutoffData(make_toy(), cutoff = NA)
  expect_equal(out$pfs_tte, c(3, 2, 6))
  expect_equal(out$pfs_event, c(1L, 1L, 0L))
  expect_equal(out$response_obs, c(1L, 0L, 1L))
})

test_that("a subset of simulated trials takes one cutoff per trial present", {
  d <- make_toy()
  d$sim <- d$sim + 1
  out <- CutoffData(d, cutoff = c(6, 4))
  ref <- CutoffData(make_toy(), cutoff = c(6, 4))
  expect_equal(out$sim, c(2, 2, 3))
  expect_equal(out[names(out) != "sim"], ref[names(ref) != "sim"])
})

test_that("a landmark response is observed only after the first assessment", {
  d <- data.frame(sim = 1, group = "A", accrual_time = c(0, 0, 0, 9),
                  pfs_time = 5, os_time = 8, progression = 1,
                  response = c(1, 1, 0, 1), ttr = NA_real_,
                  dropout_time = c(Inf, 1, Inf, Inf), pfs_cens_time = Inf,
                  resp_tau = 1.5)
  out <- CutoffData(d, cutoff = 10)
  # patient 2 dropped out before the assessment at month 1.5, and patient 4
  # entered 1 month before the cutoff
  expect_equal(out$response_obs, c(1L, 0L, 0L, 0L))
  expect_equal(out$response_eval, c(1L, 0L, 1L, 0L))
})

test_that("censoring of PFS censors PFS and the response but not OS", {
  d <- data.frame(sim = 1, group = "A", accrual_time = 0, pfs_time = 6,
                  os_time = 12, progression = 1, response = 1, ttr = c(2, 4),
                  dropout_time = Inf, pfs_cens_time = 3, resp_tau = 0)
  out <- CutoffData(d, cutoff = NA)
  expect_equal(out$pfs_tte, c(3, 3))
  expect_equal(out$pfs_event, c(0L, 0L))
  expect_equal(out$os_tte, c(12, 12))
  expect_equal(out$os_event, c(1L, 1L))
  # the second response occurs after PFS is censored (another therapy)
  expect_equal(out$response_obs, c(1L, 0L))
  expect_true(all(is.na(out$response_eval)))
})

test_that("observed responses with generated data", {
  a <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
               os.median = 15, resp.timing = "landmark", resp.tau = 1.5)
  d <- rOncoEndpoints(nsim = 1, n = 20000, arms = a, a.time = c(0, 12),
                      d.hazard = 0.1, seed = 21)
  out <- CutoffData(d, cutoff = 12)
  seen <- out$dropout_time >= 1.5 & out$accrual_time + 1.5 <= 12
  expect_equal(out$response_eval, as.integer(seen))
  expect_equal(out$response_obs, as.integer(out$response == 1 & seen))
})

test_that("outputs of version 0.2.0 are reproduced", {
  skip_on_cran()
  ref <- reg_reference()$cutoff
  reg_compare(reg_cutoff(), ref[!grepl("^event_time", names(ref))], label = "cutoff")
})

test_that("invalid inputs give errors", {
  expect_error(CutoffData(make_toy(), cutoff = c(1, 2, 3)), "one value per")
  expect_error(CutoffData(data.frame(a = 1), cutoff = 1), "rOncoEndpoints")
})
