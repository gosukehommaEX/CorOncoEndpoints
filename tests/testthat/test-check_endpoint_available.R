test_that("requests for a missing endpoint or response give errors", {
  po <- OncoArm(pfs.median = 6, death.prop = 0.15, os.median = 15,
                endpoints = c("pfs", "os"))
  expect_true(.check_endpoint_available(po, "os"))
  expect_error(.check_endpoint_available(po, "os", "responders"), "no response")
  osr <- OncoArm(os.median = 15, orr = 0.3, resp.cor = 0.35, endpoints = c("os", "response"))
  expect_error(.check_endpoint_available(osr, "pfs"), "no PFS")
})
