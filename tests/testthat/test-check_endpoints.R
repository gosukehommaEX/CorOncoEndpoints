test_that("endpoints are put in a fixed order and checked", {
  expect_identical(.check_endpoints(c("response", "pfs")), c("pfs", "response"))
  expect_identical(.check_endpoints(c("os", "pfs", "os")), c("pfs", "os"))
  expect_identical(.check_endpoints(c("pfs", "os", "response")), c("pfs", "os", "response"))
  expect_error(.check_endpoints("os"), "at least two")
  expect_error(.check_endpoints(c("pfs", "ttp")), "must contain")
  expect_error(.check_endpoints(c("pfs", NA)), "must contain")
})
