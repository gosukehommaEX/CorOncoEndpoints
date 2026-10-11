test_that("rates are recycled to the number of groups and checked", {
  expect_equal(.group_rates(0.1, 3, "x"), c(0.1, 0.1, 0.1))
  expect_equal(.group_rates(c(0, 0.2), 2, "x"), c(0, 0.2))
  expect_error(.group_rates(c(0.1, 0.2), 3, "x"), "'x' must be")
  expect_error(.group_rates(-1, 1, "x"), "non-negative")
})
