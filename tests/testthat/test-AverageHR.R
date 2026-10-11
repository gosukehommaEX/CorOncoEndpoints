test_that("average hazard ratio of OS matches independent reference values", {
  arms <- design_arms()
  ah <- AverageHR(arms, n = c(350, 350), a.time = c(0, 24), endpoint = "os",
                  time = c(12, 24, 36))
  expect_equal(ah$ahr, c(0.728589, 0.740479, 0.750512), tolerance = 1e-5)
  ahd <- AverageHR(arms, n = c(350, 350), a.time = c(0, 24), d.hazard = 0.01,
                   endpoint = "os", time = c(12, 24, 36))
  expect_equal(ahd$ahr, c(0.728259, 0.739612, 0.749039), tolerance = 1e-5)
})

test_that("a constant hazard ratio is returned unchanged", {
  arms <- design_arms()
  ah <- AverageHR(arms, n = c(350, 350), a.time = c(0, 24), endpoint = "pfs",
                  time = c(10, 30))
  expect_equal(ah$ahr, c(0.7, 0.7), tolerance = 1e-10)
})

test_that("average hazard ratio with accrual rates and piecewise dropout", {
  # reference values: inst/validation/python/design_ref.py (version 0.3.0 section)
  arms <- design_arms()
  args <- list(arms = arms, n = c(350, 350), a.time = c(0, 6), a.rate = c(10, 40),
               d.hazard = c(0.02, 0.005), d.time = c(0, 12, Inf))
  ah <- do.call(AverageHR, c(args, list(endpoint = "os", time = c(12, 24, 36))))
  expect_equal(ah$ahr, c(0.724166, 0.736230, 0.747681), tolerance = 1e-5)
  ap <- do.call(AverageHR, c(args, list(endpoint = "pfs", time = 24,
                                        pfs.cens.hazard = c(0.01, 0.03))))
  expect_equal(ap$ahr, 0.7, tolerance = 1e-10)
})

test_that("outputs of version 0.2.0 are reproduced", {
  skip_on_cran()
  reg_compare(reg_ahr(), reg_reference()$ahr, label = "AverageHR")
})

test_that("two groups are required", {
  expect_error(AverageHR(arm_A(), n = 100, time = 10), "two OncoArm")
})
