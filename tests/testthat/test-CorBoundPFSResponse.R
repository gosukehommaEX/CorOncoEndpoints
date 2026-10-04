test_that("Frechet bounds without a landmark", {
  b <- CorBoundPFSResponse(orr = 0.3)
  expect_equal(unname(b["upper"]), -sqrt(0.3 / 0.7) * log(0.3), tolerance = 1e-12)
  expect_equal(unname(b["lower"]), sqrt(0.7 / 0.3) * log(0.7), tolerance = 1e-12)
  expect_equal(unname(b), c(-0.5448299764, 0.7881852158), tolerance = 1e-9)
})

test_that("lower bound with a landmark and agreement with the Gaussian copula limit", {
  b <- CorBoundPFSResponse(orr = 0.3, pfs.median = 6, resp.tau = 1.5)
  expect_equal(unname(b["lower"]), -0.4073680014, tolerance = 1e-9)
  expect_equal(unname(b["upper"]), 0.7881852158, tolerance = 1e-9)
  # Gaussian copula at theta = -0.9999 (independent value -0.4073044647)
  z_tau <- qnorm(1 - exp(-log(2) / 6 * 1.5))
  v <- .cor_pfs_resp(-0.9999, 0.3, z_tau)$cor
  expect_equal(v, -0.4073044647, tolerance = 1e-7)
  expect_lt(abs(v - b[["lower"]]), 1e-3)
})

test_that("invalid inputs give errors", {
  expect_error(CorBoundPFSResponse(orr = 0), "orr")
  expect_error(CorBoundPFSResponse(orr = 0.3, resp.tau = 1), "exactly one")
  expect_error(CorBoundPFSResponse(orr = 0.9, pfs.median = 1, resp.tau = 2), "smaller")
})
