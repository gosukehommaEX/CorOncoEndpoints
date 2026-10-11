test_that("cumulative hazard of one or several pieces", {
  two <- list(time = c(0, 6), hazard = c(0.1, 0.02), cum = c(0, 0.6))
  expect_equal(.piecewise_cumhaz(c(0, 3, 6, 10), two), c(0, 0.3, 0.6, 0.68))
  one <- list(time = 0, hazard = 0.05, cum = 0)
  expect_identical(.piecewise_cumhaz(c(2, 4), one), 0.05 * c(2, 4))
})
