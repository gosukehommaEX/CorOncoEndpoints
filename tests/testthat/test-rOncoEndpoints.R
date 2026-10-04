test_that("output structure and ordering", {
  arms <- design_arms()
  d <- rOncoEndpoints(nsim = 3, n = c(5, 4), arms = arms, a.time = c(0, 12), seed = 1)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 3 * 9)
  expect_named(d, c("sim", "group", "accrual_time", "pfs_time", "os_time",
                    "progression", "response", "ttr", "dropout_time", "pfs_tte",
                    "pfs_event", "os_tte", "os_event", "pfs_calendar_time",
                    "os_calendar_time"))
  expect_equal(d$sim, rep(1:3, each = 9))
  expect_equal(d$group, rep(rep(c("Control", "Treatment"), c(5, 4)), 3))
  expect_true(all(d$accrual_time >= 0 & d$accrual_time <= 12))
})

test_that("the same seed reproduces the data and another seed does not", {
  a <- arm_A()
  d1 <- rOncoEndpoints(nsim = 2, n = 50, arms = a, seed = 123)
  d2 <- rOncoEndpoints(nsim = 2, n = 50, arms = a, seed = 123)
  d3 <- rOncoEndpoints(nsim = 2, n = 50, arms = a, seed = 124)
  expect_identical(d1, d2)
  expect_false(identical(d1$pfs_time, d3$pfs_time))
})

test_that("structural properties hold for every patient", {
  for (a in list(arm_A(), arm_B(), arm_C())) {
    d <- rOncoEndpoints(nsim = 1, n = 20000, arms = a, seed = 7)
    expect_true(all(d$pfs_time <= d$os_time))
    expect_true(all((d$os_time == d$pfs_time) == (d$progression == 0)))
    if (a$tau > 0) expect_true(all(d$pfs_time[d$response == 1] > a$tau))
    if (a$resp.timing == "ttr") {
      r <- d$response == 1
      expect_true(all(d$ttr[r] >= a$tau & d$ttr[r] < d$pfs_time[r]))
      expect_true(all(is.na(d$ttr[!r])))
    } else {
      expect_true(all(is.na(d$ttr)))
    }
  }
})

test_that("Monte Carlo moments agree with the theory (illness-death model)", {
  a <- arm_A()
  d <- rOncoEndpoints(nsim = 1, n = 200000, arms = a, seed = 2026)
  # Monte Carlo standard errors are about 0.001 for the proportions and the
  # correlations and about 0.03 to 0.05 for the medians
  expect_lt(abs(mean(d$response) - 0.30), 0.005)
  expect_lt(abs(mean(d$progression == 0) - 0.15), 0.005)
  expect_lt(abs(median(d$pfs_time) - 6), 0.15)
  expect_lt(abs(median(d$os_time) - 15), 0.25)
  cr <- CorEndpoints(a)
  expect_lt(abs(cor(d$pfs_time, d$response) - cr[["cor.pfs.resp"]]), 0.01)
  expect_lt(abs(cor(d$os_time, d$response) - cr[["cor.os.resp"]]), 0.01)
  expect_lt(abs(cor(d$pfs_time, d$os_time) - cr[["cor.pfs.os"]]), 0.01)
})

test_that("Monte Carlo moments agree with the theory (exponential-exponential model)", {
  cc <- arm_C()
  d <- rOncoEndpoints(nsim = 1, n = 200000, arms = cc, seed = 2027)
  expect_lt(abs(mean(d$response) - 0.30), 0.005)
  expect_lt(abs(mean(d$progression == 0) - 0.15), 0.005)
  # OS is exponential with median 15 (mean 15 / log(2))
  expect_lt(abs(mean(d$os_time) - 15 / log(2)), 0.3)
  expect_lt(abs(mean(d$os_time > 15) - 0.5), 0.005)
  expect_lt(abs(mean(d$os_time > 30) - 0.25), 0.005)
  cr <- CorEndpoints(cc)
  expect_lt(abs(cor(d$os_time, d$response) - cr[["cor.os.resp"]]), 0.01)
  expect_lt(abs(cor(d$pfs_time, d$os_time) - cr[["cor.pfs.os"]]), 0.01)
  # median of the untruncated time to response is 2.5; truncation shortens it
  expect_lt(median(d$ttr, na.rm = TRUE), 2.5)
})

test_that("accrual and dropout follow their distributions", {
  a <- arm_A()
  d <- rOncoEndpoints(nsim = 1, n = 100000, arms = a, a.time = c(0, 6, 24),
                      a.rate = c(1, 3), d.hazard = 0.05, seed = 3)
  # P(accrual < 6) = 6 / (6 + 3 * 18) = 0.1
  expect_lt(abs(mean(d$accrual_time < 6) - 0.1), 0.005)
  # P(dropout before PFS) = 0.05 / (0.05 + log(2) / 6)
  expect_lt(abs(mean(d$pfs_event == 0) - 0.05 / (0.05 + log(2) / 6)), 0.006)
  expect_equal(d$pfs_tte, pmin(d$pfs_time, d$dropout_time))
  expect_equal(d$os_calendar_time, d$accrual_time + d$os_tte)
})

test_that("invalid inputs give errors", {
  a <- arm_A()
  expect_error(rOncoEndpoints(nsim = 0, n = 10, arms = a), "nsim")
  expect_error(rOncoEndpoints(n = c(10, 10), arms = a), "one element per group")
  expect_error(rOncoEndpoints(n = 10, arms = list(1)), "OncoArm")
  expect_error(rOncoEndpoints(n = 10, arms = a, d.hazard = -1), "d.hazard")
  expect_error(rOncoEndpoints(n = 10, arms = a, a.rate = 1), "a.time")
})
