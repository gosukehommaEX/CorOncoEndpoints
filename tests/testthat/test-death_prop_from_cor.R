test_that("the proportion of deaths for a target correlation of PFS and OS", {
  # reference value: Python root finding with reference.py
  build <- function(dp) {
    OncoArm(pfs.median = 6, death.prop = dp, os.median = 15, endpoints = c("pfs", "os"))
  }
  upper <- 0.5 / (1 - exp(-log(2) / 6 * 15))
  dp <- .death_prop_from_cor(0.55, build, c(1e-4, 1 - 1e-3) * upper)
  expect_equal(dp, 0.2405415229, tolerance = 1e-6)
  expect_error(.death_prop_from_cor(0.9, build, c(1e-4, 1 - 1e-3) * upper),
               "attainable range")
})
