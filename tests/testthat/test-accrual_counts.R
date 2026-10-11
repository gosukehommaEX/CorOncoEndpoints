test_that("largest-remainder counts add up to the group size", {
  expect_identical(.accrual_counts(7, rep(1 / 3, 3)), c(3L, 2L, 2L))
  expect_identical(.accrual_counts(100, c(0.29, 0.71)), c(29L, 71L))
  expect_identical(.accrual_counts(10, c(0.06, 0.94)), c(1L, 9L))
  expect_identical(.accrual_counts(3, c(0, 1)), c(0L, 3L))
  expect_identical(.accrual_counts(5, numeric(0)), integer(0))
})
