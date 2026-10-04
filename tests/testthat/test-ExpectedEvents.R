test_that("expected events match independent reference values", {
  arms <- design_arms()
  ee <- ExpectedEvents(arms, n = c(350, 350), a.time = c(0, 24), endpoint = "os",
                       time = c(12, 24, 36))
  expect_named(ee, c("time", "events", "Control", "Treatment"))
  expect_equal(ee$events, c(58.143235, 242.259173, 449.658220), tolerance = 1e-5)
  expect_equal(ee$Control + ee$Treatment, ee$events, tolerance = 1e-12)
  ep <- ExpectedEvents(arms, n = c(350, 350), a.time = c(0, 24), endpoint = "pfs",
                       time = c(12, 24, 36))
  expect_equal(ep$events, c(143.321157, 427.211127, 611.890605), tolerance = 1e-5)
})

test_that("expected events with dropout", {
  arms <- design_arms()
  ee <- ExpectedEvents(arms, n = c(350, 350), a.time = c(0, 24), d.hazard = 0.01,
                       endpoint = "os", time = c(12, 24, 36))
  expect_equal(ee$events, c(55.594701, 223.470780, 400.653395), tolerance = 1e-5)
})

test_that("all patients entering at time 0 gives n F(t) for exponential PFS", {
  a <- arm_A()
  ee <- ExpectedEvents(a, n = 100, endpoint = "pfs", time = c(3, 6, 10))
  expect_equal(ee$events, 100 * (1 - exp(-log(2) / 6 * c(3, 6, 10))), tolerance = 1e-5)
})

test_that("invalid inputs give errors", {
  expect_error(ExpectedEvents(arm_A(), n = 100, time = -1), "positive")
  expect_error(ExpectedEvents(arm_A(), n = c(100, 100), time = 5), "one per group")
  expect_error(ExpectedEvents(arm_A(), n = 100, a.time = c(1, 5), time = 5), "a.time")
})
