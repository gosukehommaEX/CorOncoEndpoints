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

test_that("print method for groups with two endpoints", {
  expect_output(print(OncoArm(pfs.median = 6, death.prop = 0.15, os.median = 15,
                              endpoints = c("pfs", "os"))), "Corr\\(PFS, OS\\)")
  expect_output(print(OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4,
                              endpoints = c("pfs", "response"))), "Corr\\(PFS, R\\)")
  expect_output(print(OncoArm(os.median = 15, orr = 0.3, resp.cor = 0.35,
                              endpoints = c("os", "response"))), "OS hazard")
})
