test_that("required OS events match independent reference values", {
  arms <- design_arms()
  r <- RequiredEvents(arms, n = c(350, 350), a.time = c(0, 24), endpoint = "os")
  expect_identical(r$events, 366L)
  expect_equal(r$time, 30.615, tolerance = 1e-4)
  expect_equal(r$ahr, 0.746035, tolerance = 1e-5)
  # Schoenfeld formula with the returned average hazard ratio
  z2 <- (qnorm(0.975) + qnorm(0.8))^2
  expect_equal(r$events.exact, z2 / (0.25 * log(r$ahr)^2), tolerance = 1e-12)
})

test_that("required OS events with dropout", {
  arms <- design_arms()
  r <- RequiredEvents(arms, n = c(350, 350), a.time = c(0, 24), d.hazard = 0.01,
                      endpoint = "os")
  expect_identical(r$events, 370L)
  expect_equal(r$time, 33.435, tolerance = 1e-4)
  expect_equal(r$ahr, 0.747109, tolerance = 1e-5)
})

test_that("PFS with a constant hazard ratio gives the Schoenfeld number", {
  arms <- design_arms()
  r <- RequiredEvents(arms, n = c(350, 350), a.time = c(0, 24), endpoint = "pfs")
  expect_identical(r$events, 247L)
  expect_equal(r$ahr, 0.7, tolerance = 1e-10)
  expect_equal(r$time, 16.761, tolerance = 1e-4)
})

test_that("required events with accrual rates, piecewise dropout and censoring of PFS", {
  # reference values: inst/validation/python/design_ref.py (version 0.3.0 section)
  arms <- design_arms()
  args <- list(arms = arms, n = c(350, 350), a.time = c(0, 6), a.rate = c(10, 40),
               d.hazard = c(0.02, 0.005), d.time = c(0, 12, Inf))
  r <- do.call(RequiredEvents, c(args, list(endpoint = "os", power = 0.85)))
  expect_identical(r$events, 457L)
  expect_equal(r$time, 47.36872, tolerance = 1e-4)
  expect_equal(r$ahr, 0.755476, tolerance = 1e-5)
  p <- do.call(RequiredEvents, c(args, list(endpoint = "pfs",
                                            pfs.cens.hazard = c(0.01, 0.03))))
  expect_identical(p$events, 247L)
  expect_equal(p$time, 19.49442, tolerance = 1e-4)
})

test_that("outputs of version 0.2.0 are reproduced", {
  skip_on_cran()
  reg_compare(reg_required(), reg_reference()$required, label = "RequiredEvents")
})

test_that("an unattainable number of events gives an error", {
  arms <- design_arms()
  expect_error(RequiredEvents(arms, n = c(50, 50), a.time = c(0, 24), endpoint = "os"),
               "exceeds")
})
