test_that("print method reports parameters and implied quantities", {
  a <- arm_A()
  expect_output(print(a), "OncoArm")
  expect_output(print(a), "Implied correlations")
  # capture the printed text so that it does not appear in the test log
  out <- capture.output(vis <- withVisible(print(a)))
  expect_false(vis$visible)
  expect_identical(vis$value, a)
  expect_true(any(grepl("Partial correlation", out)))
})
