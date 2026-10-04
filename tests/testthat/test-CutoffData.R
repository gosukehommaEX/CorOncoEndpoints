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

test_that("invalid inputs give errors", {
  expect_error(CutoffData(make_toy(), cutoff = c(1, 2, 3)), "one value per")
  expect_error(CutoffData(data.frame(a = 1), cutoff = 1), "rOncoEndpoints")
})
