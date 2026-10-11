test_that("groups without an endpoints element have all three endpoints", {
  a <- unclass(arm_A())
  a$endpoints <- NULL
  class(a) <- "OncoArm"
  expect_identical(.arm_endpoints(a), c("pfs", "os", "response"))
  po <- OncoArm(pfs.median = 6, death.prop = 0.15, os.median = 15,
                endpoints = c("os", "pfs"))
  expect_identical(.arm_endpoints(po), c("pfs", "os"))
})
