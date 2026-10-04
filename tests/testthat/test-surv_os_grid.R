test_that("grid values agree with adaptive quadrature", {
  arms <- list(arm_A(), design_arms()$Treatment,
               OncoArm(pfs.median = 6, orr = 0.3, resp.pfs.median = 11,
                       death.prop = 0.15, pps.median = 8, pps.hr.resp = 0.5,
                       resp.timing = "landmark", resp.tau = 1.5))
  tt <- seq(0, 120, length.out = 801)
  pick <- c(2, 3, 11, 100, 400, 801)
  for (a in arms) {
    g <- .surv_os_grid(tt, a)
    expect_equal(g$S[pick], .surv_os(tt[pick], a, type = "survival"), tolerance = 1e-9)
    expect_equal(g$f[pick], .surv_os(tt[pick], a, type = "density"), tolerance = 1e-9)
    expect_equal(g$S[1], 1)
  }
})

test_that("closed forms are used without response-dependent post-progression survival", {
  for (a in list(arm_B(), arm_C())) {
    tt <- c(0, 2, 12, 30)
    g <- .surv_os_grid(tt, a)
    expect_identical(g$S, .surv_os(tt, a, type = "survival"))
    expect_identical(g$f, .surv_os(tt, a, type = "density"))
  }
})

test_that("a grid that does not start at 0 gives an error", {
  expect_error(.surv_os_grid(c(1, 2), arm_A()), "start at 0")
})
