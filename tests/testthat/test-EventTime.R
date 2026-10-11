test_that("calendar time of the k-th event per simulated trial", {
  d <- data.frame(sim = c(1, 1, 1, 2, 2), os_event = c(1, 1, 0, 1, 0),
                  os_calendar_time = c(5, 3, 1, 7, 2), pfs_event = c(1, 1, 1, 1, 1),
                  pfs_calendar_time = c(4, 2, 1, 3, 2))
  expect_equal(unname(EventTime(d, 1, "os")), c(3, 7))
  expect_equal(unname(EventTime(d, 2, "os")), c(5, NA))
  expect_equal(unname(EventTime(d, 2, "pfs")), c(2, 3))
  expect_named(EventTime(d, 1, "os"), c("1", "2"))
})

test_that("event-driven cutoff with generated data", {
  a <- arm_A()
  d <- rOncoEndpoints(nsim = 2, n = 300, arms = a, a.time = c(0, 12), seed = 5)
  et <- EventTime(d, 150, "os")
  cut <- CutoffData(d, et)
  expect_equal(as.vector(table(cut$sim[cut$os_event == 1])), c(150L, 150L))
})

test_that("events censored by the censoring of PFS are not counted", {
  a <- arm_A()
  d <- rOncoEndpoints(nsim = 1, n = 500, arms = a, pfs.cens.hazard = 0.05, seed = 22)
  et <- EventTime(d, 100, "pfs")
  expect_equal(unname(et), sort(d$pfs_calendar_time[d$pfs_event == 1])[100])
})

test_that("outputs of version 0.2.0 are reproduced", {
  skip_on_cran()
  ref <- reg_reference()$cutoff
  cur <- reg_cutoff()
  reg_compare(cur[c("event_time_pfs", "event_time_os")],
              ref[c("event_time_pfs", "event_time_os")], label = "EventTime")
})

test_that("invalid inputs give errors", {
  expect_error(EventTime(data.frame(sim = 1), 1), "rOncoEndpoints")
  expect_error(EventTime(data.frame(sim = 1, os_event = 1, os_calendar_time = 1), 0),
               "positive integer")
})

test_that("an endpoint that the groups do not have gives an error", {
  a <- OncoArm(os.median = 15, orr = 0.3, resp.cor = 0.35, endpoints = c("os", "response"))
  d <- rOncoEndpoints(nsim = 1, n = 50, arms = a, seed = 35)
  expect_error(EventTime(d, 10, "pfs"), "no PFS")
  expect_equal(unname(EventTime(d, 10, "os")), sort(d$os_calendar_time)[10])
})
