test_that("orthant probabilities agree with closed forms", {
  # P(Z1 > 0, Z2 > 0) = 1/4 + asin(theta) / (2 pi)
  for (th in c(-0.9, -0.3, 0.5, 0.95)) {
    expect_equal(.bvn_upper(0, 0, th), 0.25 + asin(th) / (2 * pi), tolerance = 1e-9)
  }
  expect_equal(.bvn_upper(-Inf, 0.4, 0.5), pnorm(0.4, lower.tail = FALSE))
  expect_equal(.bvn_upper(0.3, -0.2, 0), pnorm(0.3, lower.tail = FALSE) *
                 pnorm(-0.2, lower.tail = FALSE))
  expect_equal(.bvn_upper(Inf, 0, 0.3), 0)
})

test_that("response threshold gives the target response rate with a landmark", {
  z_tau <- qnorm(1 - exp(-log(2) / 6 * 1.5))
  cc <- .resp_threshold(0.47, 0.3, z_tau)
  expect_equal(.bvn_upper(z_tau, cc, 0.47), 0.3, tolerance = 1e-9)
})
