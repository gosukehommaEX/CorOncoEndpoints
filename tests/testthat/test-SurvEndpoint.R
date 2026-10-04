test_that("OS survival, density and hazard match independent reference values", {
  a <- arm_A()
  expect_equal(SurvEndpoint(a, 12, "os"), 0.5966965852, tolerance = 1e-7)
  expect_equal(SurvEndpoint(a, 12, "os", "responders"), 0.8022802053, tolerance = 1e-7)
  expect_equal(SurvEndpoint(a, 12, "os", "nonresponders"), 0.5085893195, tolerance = 1e-7)
  expect_equal(SurvEndpoint(a, 12, "os", type = "density"), 0.03415425313, tolerance = 1e-7)
  expect_equal(SurvEndpoint(a, 12, "os", type = "hazard"), 0.0572388949, tolerance = 1e-7)
  b <- arm_B()
  expect_equal(SurvEndpoint(b, 12, "os"), 0.602081528, tolerance = 1e-7)
  expect_equal(SurvEndpoint(b, 12, "os", "responders"), 0.7547041168, tolerance = 1e-7)
  expect_equal(SurvEndpoint(b, 12, "os", "nonresponders"), 0.5366718471, tolerance = 1e-7)
  expect_equal(SurvEndpoint(b, 12, "os", type = "density"), 0.03483770969, tolerance = 1e-7)
  cc <- arm_C()
  # the expexp values depend on the tabulated post-progression hazard (finer
  # grid in the Python reference), hence the looser tolerance
  expect_equal(SurvEndpoint(cc, 12, "os", "responders"), 0.7672331998, tolerance = 1e-5)
  expect_equal(SurvEndpoint(cc, 12, "os", "nonresponders"), 0.4916845913, tolerance = 1e-5)
})

test_that("OS is exactly exponential in the exponential-exponential model", {
  cc <- arm_C()
  tt <- c(1, 6, 12, 30)
  expect_equal(SurvEndpoint(cc, tt, "os"), exp(-log(2) / 15 * tt), tolerance = 1e-12)
  expect_equal(SurvEndpoint(cc, tt, "os", type = "hazard"), rep(log(2) / 15, 4),
               tolerance = 1e-12)
  # the mixture over response groups also gives the exponential survival
  s_mix <- 0.3 * SurvEndpoint(cc, tt, "os", "responders") +
    0.7 * SurvEndpoint(cc, tt, "os", "nonresponders")
  expect_equal(s_mix, exp(-log(2) / 15 * tt), tolerance = 1e-5)
})

test_that("PFS survival and density overall and by response", {
  a <- arm_A()
  tt <- c(0, 3, 12)
  expect_equal(SurvEndpoint(a, tt, "pfs"), exp(-log(2) / 6 * tt), tolerance = 1e-12)
  expect_equal(SurvEndpoint(a, 12, "pfs", "responders"), 0.4745194414, tolerance = 1e-7)
  expect_equal(SurvEndpoint(a, 12, "pfs", "responders", "density"), 0.04081042162,
               tolerance = 1e-7)
  expect_equal(SurvEndpoint(a, 0, "pfs", "responders"), 1, tolerance = 1e-10)
  s_mix <- 0.3 * SurvEndpoint(a, tt, "pfs", "responders") +
    0.7 * SurvEndpoint(a, tt, "pfs", "nonresponders")
  expect_equal(s_mix, exp(-log(2) / 6 * tt), tolerance = 1e-10)
  # with a landmark, responders do not have PFS events before resp.tau
  b <- arm_B()
  expect_equal(SurvEndpoint(b, c(0.5, 1.4), "pfs", "responders"), c(1, 1), tolerance = 1e-10)
})

test_that("mixture of response groups equals the overall OS survival", {
  for (a in list(arm_A(), arm_B())) {
    tt <- c(2, 10, 25)
    s_mix <- a$orr * SurvEndpoint(a, tt, "os", "responders") +
      (1 - a$orr) * SurvEndpoint(a, tt, "os", "nonresponders")
    expect_equal(s_mix, SurvEndpoint(a, tt, "os"), tolerance = 1e-8)
  }
})

test_that("Fleischer et al. (2009) Theorem 5 closed form is used without response effects", {
  f <- arm_fleischer(0.342, 0.0654, 0.0797)
  a <- 0.342 + 0.0654
  tt <- c(3, 9, 20)
  s5 <- 0.342 / (a - 0.0797) * exp(-0.0797 * tt) -
    (0.0797 - 0.0654) / (a - 0.0797) * exp(-a * tt)
  expect_equal(SurvEndpoint(f, tt, "os"), s5, tolerance = 1e-12)
})

test_that("invalid times give errors", {
  expect_error(SurvEndpoint(arm_A(), -1), "non-negative")
})
