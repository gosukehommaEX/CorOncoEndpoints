test_that("calibration reproduces independent reference values (idm, correlation input)", {
  a <- arm_A()
  expect_s3_class(a, "OncoArm")
  expect_equal(a$lam_p, log(2) / 6, tolerance = 1e-12)
  expect_equal(a$theta, 0.5375882765, tolerance = 1e-7)
  expect_equal(a$c_resp, 0.5244005127, tolerance = 1e-7)
  expect_equal(a$gam0, 0.09702488823, tolerance = 1e-7)
  expect_equal(a$gam1, 0.6 * a$gam0, tolerance = 1e-12)
  expect_equal(a$z_tau, -Inf)
})

test_that("calibration with landmark and responders' median PFS", {
  b <- arm_B()
  expect_equal(b$theta, 0.4702658637, tolerance = 1e-7)
  expect_equal(b$c_resp, 0.4824116033, tolerance = 1e-7)
  expect_equal(b$z_tau, -0.9981488825, tolerance = 1e-8)
  expect_equal(b$gam0, log(2) / 8, tolerance = 1e-12)
  expect_equal(b$gam1, b$gam0, tolerance = 1e-12)
})

test_that("exponential-exponential model parameters", {
  cc <- arm_C()
  expect_equal(cc$theta, 0.5001446401, tolerance = 1e-7)
  expect_equal(cc$c_resp, 0.4878836073, tolerance = 1e-7)
  expect_equal(cc$lam_o, log(2) / 15, tolerance = 1e-12)
  expect_equal(cc$c_dec, 0.1925408835, tolerance = 1e-8)
  expect_equal(cc$ttr_rate, log(2) / (2.5 - 1.5), tolerance = 1e-12)
  # death.prop equal to mPFS / mOS gives Fleischer's maximal independence model
  d <- arm_D()
  expect_equal(d$c_dec, 0, tolerance = 1e-12)
  # the tabulated cumulative hazard is then lam_o * t
  expect_equal(d$grid$G[1001], d$lam_o * d$grid$t[1001], tolerance = 1e-10)
})

test_that("the post-progression hazard is calibrated to the OS median", {
  a <- arm_A()
  expect_equal(SurvEndpoint(a, 15, endpoint = "os"), 0.5, tolerance = 1e-8)
})

test_that("a zero correlation gives a zero copula correlation", {
  a <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0, death.prop = 0.15,
               os.median = 15)
  expect_equal(a$theta, 0, tolerance = 1e-8)
})

test_that("invalid inputs give errors", {
  expect_error(OncoArm(pfs.median = 6, pfs.hazard = 0.1, orr = 0.3, resp.cor = 0.4,
                       death.prop = 0.1, os.median = 15), "exactly one")
  expect_error(OncoArm(pfs.median = 6, orr = 1.2, resp.cor = 0.4,
                       death.prop = 0.1, os.median = 15), "orr")
  expect_error(OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, resp.pfs.median = 10,
                       death.prop = 0.1, os.median = 15), "exactly one")
  expect_error(OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4,
                       death.prop = 0.1, os.median = 5), "larger than the PFS median")
  expect_error(OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.1,
                       os.median = 15, pps.median = 8), "exactly one")
  expect_error(OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.1,
                       os.model = "expexp", os.median = 15, pps.hr.resp = 0.5),
               "must be 1")
  expect_error(OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.5,
                       os.model = "expexp", os.median = 15), "death.prop")
  expect_error(OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.1,
                       os.median = 15, resp.tau = 1), "must be 0")
  expect_error(OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.1,
                       os.median = 15, resp.timing = "landmark"), "positive")
  expect_error(OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.1,
                       os.median = 15, resp.timing = "ttr"), "ttr.median")
  expect_error(OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.79, death.prop = 0.1,
                       os.median = 15), "attainable range")
  # without deaths after progression S_OS(20) = 1 - 0.6 * (1 - 2^(-20 / 6)) < 0.5
  expect_error(OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.6,
                       os.median = 20), "cannot be attained")
  expect_error(OncoArm(pfs.median = 6, orr = NA, resp.cor = 0.4, death.prop = 0.1,
                       os.median = 15), "orr")
  expect_error(OncoArm(pfs.median = 6, orr = 0.3, resp.cor = NA, death.prop = 0.1,
                       os.median = 15), "resp.cor")
  expect_error(OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = NA,
                       os.median = 15), "death.prop")
})
