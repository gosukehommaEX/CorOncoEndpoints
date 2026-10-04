test_that("print method reports parameters and implied quantities", {
  a <- arm_A()
  expect_output(print(a), "OncoArm")
  expect_output(print(a), "Implied correlations")
  expect_invisible(print(a))
})
