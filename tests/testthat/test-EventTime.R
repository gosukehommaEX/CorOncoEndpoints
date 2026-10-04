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

test_that("invalid inputs give errors", {
  expect_error(EventTime(data.frame(sim = 1), 1), "rOncoEndpoints")
  expect_error(EventTime(data.frame(sim = 1, os_event = 1, os_calendar_time = 1), 0),
               "positive integer")
})
