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

test_that("expected events with accrual rates, piecewise dropout and censoring of PFS", {
  # reference values: inst/validation/python/design_ref.py (version 0.3.0 section)
  arms <- design_arms()
  # 10 per month on [0, 6), then 40 per month until 700 patients (month 22);
  # dropout hazard 0.02 before month 12 and 0.005 after
  args <- list(arms = arms, n = c(350, 350), a.time = c(0, 6), a.rate = c(10, 40),
               d.hazard = c(0.02, 0.005), d.time = c(0, 12, Inf), time = c(12, 24, 36))
  eo <- do.call(ExpectedEvents, c(args, list(endpoint = "os")))
  expect_equal(eo$events, c(30.503955, 194.554769, 366.650184), tolerance = 1e-5)
  ep <- do.call(ExpectedEvents, c(args, list(endpoint = "pfs",
                                             pfs.cens.hazard = c(0.01, 0.03))))
  expect_equal(ep$events, c(84.368845, 355.003186, 472.098251), tolerance = 1e-5)
  # the censoring of PFS does not change the expected number of OS events
  eo2 <- do.call(ExpectedEvents, c(args, list(endpoint = "os", pfs.cens.hazard = 0.05)))
  expect_equal(eo2$events, eo$events)
})

test_that("accrual that starts after time 0 shifts the expected events", {
  arms <- design_arms()
  e3 <- ExpectedEvents(arms, n = c(350, 350), a.time = c(3, 27), endpoint = "os",
                       time = c(12, 24, 36))
  expect_equal(e3$events, c(30.917577, 186.629330, 405.346941), tolerance = 1e-5)
})

test_that("outputs of version 0.2.0 are reproduced", {
  skip_on_cran()
  reg_compare(reg_expected(), reg_reference()$expected, label = "expected")
})

test_that("invalid inputs give errors", {
  expect_error(ExpectedEvents(arm_A(), n = 100, time = -1), "positive")
  expect_error(ExpectedEvents(arm_A(), n = c(100, 100), time = 5), "one per group")
  expect_error(ExpectedEvents(arm_A(), n = 100, a.time = c(5, 1), time = 5), "a.time")
  expect_error(ExpectedEvents(arm_A(), n = 100, a.time = c(0, 10), a.rate = 5, time = 5),
               "accrues 50")
})

test_that("expected OS events of a group with OS and response", {
  a <- OncoArm(os.median = 15, orr = 0.3, resp.cor = 0.35, endpoints = c("os", "response"))
  ee <- ExpectedEvents(a, n = 100, endpoint = "os", time = c(6, 15))
  expect_equal(ee$events, 100 * (1 - exp(-log(2) / 15 * c(6, 15))), tolerance = 1e-5)
  expect_error(ExpectedEvents(a, n = 100, endpoint = "pfs", time = 6), "no PFS")
})
