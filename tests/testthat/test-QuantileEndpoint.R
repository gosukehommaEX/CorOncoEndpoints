test_that("medians match independent reference values", {
  a <- arm_A()
  expect_equal(QuantileEndpoint(a, endpoint = "os"), 15, tolerance = 1e-7)
  expect_equal(QuantileEndpoint(a, endpoint = "os", response = "responders"), 24.06564768,
               tolerance = 1e-6)
  expect_equal(QuantileEndpoint(a, endpoint = "os", response = "nonresponders"), 12.22740831,
               tolerance = 1e-6)
  expect_equal(QuantileEndpoint(a, endpoint = "pfs", response = "responders"), 11.38721409,
               tolerance = 1e-6)
  expect_equal(QuantileEndpoint(a, endpoint = "pfs", response = "nonresponders"), 4.39304265,
               tolerance = 1e-6)
  b <- arm_B()
  expect_equal(QuantileEndpoint(b, endpoint = "pfs", response = "responders"), 11,
               tolerance = 1e-7)
  expect_equal(QuantileEndpoint(b, endpoint = "os"), 15.08537599, tolerance = 1e-6)
  cc <- arm_C()
  expect_equal(QuantileEndpoint(cc, endpoint = "os"), 15, tolerance = 1e-10)
  expect_equal(QuantileEndpoint(cc, endpoint = "os", response = "responders"), 23.68403189,
               tolerance = 1e-5)
  expect_equal(QuantileEndpoint(cc, endpoint = "pfs", response = "responders"), 11.28857374,
               tolerance = 1e-6)
})

test_that("several probabilities and PFS quantiles", {
  a <- arm_A()
  expect_equal(QuantileEndpoint(a, c(0.25, 0.5), endpoint = "pfs"),
               -log(c(0.75, 0.5)) / (log(2) / 6), tolerance = 1e-12)
})

test_that("Fleischer et al. (2009) Example 3 medians are reproduced", {
  # control (0.342, 0.0654, 0.0797): observed median OS 9.2
  m_c <- QuantileEndpoint(arm_fleischer(0.342, 0.0654, 0.0797), endpoint = "os")
  expect_equal(round(m_c, 1), 9.2)
  # lambda_3 of the control replaced by 0.0876: published median OS 8.6
  m_s <- QuantileEndpoint(arm_fleischer(0.342, 0.0654, 0.0876), endpoint = "os")
  expect_equal(round(m_s, 1), 8.6)
  expect_equal(m_s, 8.631816411, tolerance = 1e-7)
  # treatment (0.139, 0.0648, 0.0876): observed median OS 9.3
  m_t <- QuantileEndpoint(arm_fleischer(0.139, 0.0648, 0.0876), endpoint = "os")
  expect_equal(round(m_t, 1), 9.3)
})

test_that("TrialSimulator documented medians (4.62 and 9.61) are reproduced", {
  ts <- arm_fleischer(0.1, 0.05, 0.12)
  expect_equal(round(QuantileEndpoint(ts, endpoint = "pfs"), 2), 4.62)
  expect_equal(round(QuantileEndpoint(ts, endpoint = "os"), 2), 9.61)
})

test_that("invalid probabilities give errors", {
  expect_error(QuantileEndpoint(arm_A(), prob = 1), "prob")
})
