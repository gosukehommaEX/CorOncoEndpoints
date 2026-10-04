test_that("correlations match independent reference values", {
  ref <- list(
    A = c(0.4, 0.3808328878, 0.6034092891, 0.1908293248),
    B = c(0.380762792, 0.2301205512, 0.604367223, 0),
    C = c(0.4, 0.2435830595, 0.5756098965, 0.01779835806),
    D = c(0.4, 0.16, 0.4, 0)
  )
  arms <- list(A = arm_A(), B = arm_B(), C = arm_C(), D = arm_D())
  for (k in names(arms)) {
    v <- CorEndpoints(arms[[k]])
    expect_named(v, c("cor.pfs.resp", "cor.os.resp", "cor.pfs.os", "pcor.os.resp"))
    expect_equal(unname(v), ref[[k]], tolerance = 1e-6, info = k)
  }
})

test_that("without response-dependent post-progression survival Corr(OS,R) = Corr(PFS,R) Corr(PFS,OS)", {
  v <- CorEndpoints(arm_B())
  expect_equal(v[["cor.os.resp"]], v[["cor.pfs.resp"]] * v[["cor.pfs.os"]], tolerance = 1e-10)
  expect_equal(v[["pcor.os.resp"]], 0, tolerance = 1e-10)
})

test_that("a longer post-progression survival of responders gives a positive partial correlation", {
  a <- arm_A()
  expect_gt(CorEndpoints(a)[["pcor.os.resp"]], 0)
  h <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
               os.median = 15, pps.hr.resp = 1.5)
  expect_lt(CorEndpoints(h)[["pcor.os.resp"]], 0)
})

test_that("Fleischer et al. (2009) published correlations are reproduced", {
  # Section 3: lambda = (2, 0.5, 1.5) gives a correlation of 0.52
  v1 <- CorEndpoints(arm_fleischer(2, 0.5, 1.5))[["cor.pfs.os"]]
  expect_equal(round(v1, 2), 0.52)
  expect_equal(v1, 1.5 / sqrt(2^2 + 2 * 2 * 0.5 + 1.5^2), tolerance = 1e-8)
  # Example 1: lambda = (0.284, 0.075, 0.128) gives a correlation of 0.34
  v2 <- CorEndpoints(arm_fleischer(0.284, 0.075, 0.128))[["cor.pfs.os"]]
  expect_equal(round(v2, 2), 0.34)
  # maximal independence model: Corr(PFS, OS) = median PFS / median OS (Theorem 1)
  expect_equal(CorEndpoints(arm_D())[["cor.pfs.os"]], 6 / 15, tolerance = 1e-6)
})

test_that("TrialSimulator documented example (h01 = 0.1, h02 = 0.05, h12 = 0.12) is reproduced", {
  v <- CorEndpoints(arm_fleischer(0.1, 0.05, 0.12))[["cor.pfs.os"]]
  expect_equal(round(v, 2), 0.65)
  expect_equal(v, 0.6469966392, tolerance = 1e-8)
})

test_that("non-OncoArm input gives an error", {
  expect_error(CorEndpoints(list()), "OncoArm")
})
