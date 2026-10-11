test_that("output structure and ordering", {
  arms <- design_arms()
  d <- rOncoEndpoints(nsim = 3, n = c(5, 4), arms = arms, a.time = c(0, 12), seed = 1)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 3 * 9)
  expect_named(d, c("sim", "group", "accrual_time", "pfs_time", "os_time",
                    "progression", "response", "ttr", "dropout_time", "pfs_tte",
                    "pfs_event", "os_tte", "os_event", "pfs_calendar_time",
                    "os_calendar_time", "pfs_cens_time", "resp_tau"))
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

test_that("accrual has a fixed number of patients in each interval", {
  a <- arm_A()
  # proportions: 10% of the patients in [0, 6) and 90% in [6, 24)
  d <- rOncoEndpoints(nsim = 2, n = c(1000, 500), arms = list(A = a, B = a),
                      a.time = c(0, 6, 24), a.prop = c(0.1, 0.9), seed = 3)
  tab <- table(d$sim, d$group, d$accrual_time < 6)
  expect_equal(as.vector(tab[, , "TRUE"]), c(100, 100, 50, 50))
  expect_true(all(d$accrual_time >= 0 & d$accrual_time < 24))
  # absolute rates: 5 per month on [0, 6) and 15 per month on [6, 24) accrue
  # 30 + 270 = 300 patients
  r <- rOncoEndpoints(nsim = 1, n = c(150, 150), arms = list(A = a, B = a),
                      a.time = c(0, 6, 24), a.rate = c(5, 15), seed = 4)
  expect_equal(as.vector(table(r$group, r$accrual_time < 6)[, "TRUE"]), c(15, 15))
  # open final interval: 10 per month on [0, 6), then 20 per month until 200
  # patients, that is, until month 6 + 140 / 20 = 13
  o <- rOncoEndpoints(nsim = 1, n = c(100, 100), arms = list(A = a, B = a),
                      a.time = c(0, 6), a.rate = c(10, 20), seed = 5)
  expect_equal(sum(o$accrual_time < 6), 60)
  expect_true(max(o$accrual_time) <= 13 && max(o$accrual_time) > 12.5)
  # one rate from time 0: 25 per month until 100 patients (month 4)
  s <- rOncoEndpoints(nsim = 1, n = 100, arms = a, a.time = 0, a.rate = 25, seed = 6)
  expect_true(all(s$accrual_time < 4))
  # accrual that starts after time 0
  e <- rOncoEndpoints(nsim = 1, n = 50, arms = a, a.time = c(12, 24), seed = 7)
  expect_true(all(e$accrual_time >= 12 & e$accrual_time < 24))
})

test_that("invalid accrual specifications give errors", {
  a <- arm_A()
  expect_error(rOncoEndpoints(n = 100, arms = a, a.time = c(0, 6, 24), a.rate = c(1, 3)),
               "accrues 60")
  expect_error(rOncoEndpoints(n = 100, arms = a, a.time = c(0, 12), a.rate = 5, a.prop = 1),
               "at most one")
  expect_error(rOncoEndpoints(n = 100, arms = a, a.time = c(0, 12), a.rate = c(10, 1)),
               "already accrue")
  expect_error(rOncoEndpoints(n = 100, arms = a, a.time = c(0, 6, 24), a.prop = 1), "a.prop")
  expect_error(rOncoEndpoints(n = 100, arms = a, a.time = 0), "at least two")
  expect_error(rOncoEndpoints(n = 100, arms = a, a.time = c(6, 0)), "increasing")
})

test_that("exponential dropout follows its distribution", {
  a <- arm_A()
  d <- rOncoEndpoints(nsim = 1, n = 100000, arms = a, a.time = c(0, 24),
                      d.hazard = 0.05, seed = 3)
  # P(dropout before PFS) = 0.05 / (0.05 + log(2) / 6)
  expect_lt(abs(mean(d$pfs_event == 0) - 0.05 / (0.05 + log(2) / 6)), 0.006)
  expect_equal(d$pfs_tte, pmin(d$pfs_time, d$dropout_time))
  expect_equal(d$os_calendar_time, d$accrual_time + d$os_tte)
})

test_that("dropout given as a median or as one piece reproduces the exponential case", {
  a <- arm_A()
  d1 <- rOncoEndpoints(nsim = 2, n = 50, arms = a, a.time = c(0, 12), d.hazard = 0.05,
                       seed = 8)
  d2 <- rOncoEndpoints(nsim = 2, n = 50, arms = a, a.time = c(0, 12),
                       d.median = log(2) / 0.05, seed = 8)
  d3 <- rOncoEndpoints(nsim = 2, n = 50, arms = a, a.time = c(0, 12), d.hazard = 0.05,
                       d.time = c(0, Inf), seed = 8)
  expect_equal(d2, d1, tolerance = 1e-12)
  expect_identical(d3, d1)
})

test_that("piecewise exponential dropout follows its distribution", {
  a <- arm_A()
  d <- rOncoEndpoints(nsim = 1, n = 200000, arms = a, d.hazard = c(0.1, 0.02),
                      d.time = c(0, 6, Inf), seed = 9)
  # P(D <= 6) = 1 - exp(-0.6) and P(D <= 12) = 1 - exp(-0.6 - 0.12)
  expect_lt(abs(mean(d$dropout_time <= 6) - (1 - exp(-0.6))), 0.004)
  expect_lt(abs(mean(d$dropout_time <= 12) - (1 - exp(-0.72))), 0.004)
  # pieces per group through lists
  g <- rOncoEndpoints(nsim = 1, n = c(100000, 100000), arms = list(A = a, B = a),
                      d.hazard = list(c(0.1, 0.02), 0),
                      d.time = list(c(0, 6, Inf), c(0, Inf)), seed = 10)
  expect_lt(abs(mean(g$dropout_time[g$group == "A"] <= 6) - (1 - exp(-0.6))), 0.006)
  expect_true(all(is.infinite(g$dropout_time[g$group == "B"])))
})

test_that("invalid dropout specifications give errors", {
  a <- arm_A()
  expect_error(rOncoEndpoints(n = 10, arms = a, d.hazard = 0.1, d.median = 5), "at most one")
  expect_error(rOncoEndpoints(n = 10, arms = a, d.hazard = c(0.1, 0.02)), "one per group")
  expect_error(rOncoEndpoints(n = 10, arms = a, d.hazard = c(0.1, 0.02), d.time = c(0, 6)),
               "Inf")
  expect_error(rOncoEndpoints(n = 10, arms = a, d.hazard = 0.1, d.time = c(0, 6, Inf)),
               "length")
  expect_error(rOncoEndpoints(n = 10, arms = a, d.median = 0), "positive medians")
  expect_error(rOncoEndpoints(n = 10, arms = a, pfs.cens.hazard = -1), "pfs.cens.hazard")
})

test_that("censoring of PFS applies to PFS only", {
  a <- arm_A()
  d0 <- rOncoEndpoints(nsim = 1, n = 100000, arms = a, d.hazard = 0.02, seed = 11)
  d1 <- rOncoEndpoints(nsim = 1, n = 100000, arms = a, d.hazard = 0.02,
                       pfs.cens.hazard = 0.05, seed = 11)
  # the censoring times are drawn after the other random numbers of the group
  same <- c("accrual_time", "pfs_time", "os_time", "response", "dropout_time",
            "os_tte", "os_event")
  expect_identical(d1[same], d0[same])
  expect_true(all(is.infinite(d0$pfs_cens_time)))
  expect_equal(d1$pfs_tte, pmin(d1$pfs_time, d1$dropout_time, d1$pfs_cens_time))
  # P(PFS censored before dropout and PFS) = c / (c + d + lambda)
  lam <- log(2) / 6
  expect_lt(abs(mean(d1$pfs_cens_time < pmin(d1$pfs_time, d1$dropout_time)) -
                  0.05 / (0.05 + 0.02 + lam)), 0.005)
})

test_that("streams give reproducible and different random numbers", {
  a <- arm_A()
  s0 <- rOncoEndpoints(nsim = 2, n = 20, arms = a, seed = 12)
  s00 <- rOncoEndpoints(nsim = 2, n = 20, arms = a, seed = 12, stream = 0)
  s1 <- rOncoEndpoints(nsim = 2, n = 20, arms = a, seed = 12, stream = 1)
  s1b <- rOncoEndpoints(nsim = 2, n = 20, arms = a, seed = 12, stream = 1)
  expect_identical(s1, s1b)
  expect_identical(s00, s0)
  expect_false(identical(s1$pfs_time, s0$pfs_time))
  expect_error(rOncoEndpoints(n = 20, arms = a, stream = 1), "requires 'seed'")
  expect_error(rOncoEndpoints(n = 20, arms = a, seed = 1, stream = -1), "whole number")
})

test_that("outputs of version 0.2.0 are reproduced", {
  skip_on_cran()
  ref <- reg_reference()$generate
  reg_compare(reg_generate(), ref, label = "generate")
  # the comparison detects a change: another seed gives other data
  other <- rOncoEndpoints(nsim = 2, n = 20, arms = reg_arms()$expexp, a.time = c(0, 24),
                          seed = 103)
  expect_false(isTRUE(all.equal(other[names(ref$expexp)], ref$expexp, tolerance = 1e-8)))
})

test_that("invalid inputs give errors", {
  a <- arm_A()
  expect_error(rOncoEndpoints(nsim = 0, n = 10, arms = a), "nsim")
  expect_error(rOncoEndpoints(n = c(10, 10), arms = a), "one element per group")
  expect_error(rOncoEndpoints(n = 10, arms = list(1)), "OncoArm")
  expect_error(rOncoEndpoints(n = 10, arms = a, d.hazard = -1), "d.hazard")
  expect_error(rOncoEndpoints(n = 10, arms = a, a.rate = 1), "a.time")
})

test_that("groups with two endpoints share the random numbers of three endpoints", {
  full <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
                  os.median = 15)
  d3 <- rOncoEndpoints(nsim = 2, n = 50, arms = full, a.time = c(0, 12), seed = 31)
  po <- OncoArm(pfs.median = 6, death.prop = 0.15, os.median = 15,
                endpoints = c("pfs", "os"))
  d_po <- rOncoEndpoints(nsim = 2, n = 50, arms = po, a.time = c(0, 12), seed = 31)
  cols <- c("pfs_time", "os_time", "progression")
  expect_identical(d_po[cols], d3[cols])
  expect_true(all(is.na(d_po$response)) && all(is.na(d_po$ttr)))
  pr <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, endpoints = c("pfs", "response"))
  d_pr <- rOncoEndpoints(nsim = 2, n = 50, arms = pr, a.time = c(0, 12), seed = 31)
  cols <- c("pfs_time", "response", "pfs_tte")
  expect_identical(d_pr[cols], d3[cols])
  expect_true(all(is.na(d_pr$os_time)) && all(is.na(d_pr$os_event)) &&
                all(is.na(d_pr$progression)))
})

test_that("OS and response: OS is exponential and linked to response", {
  a <- OncoArm(os.median = 15, orr = 0.3, resp.cor = 0.35, resp.timing = "ttr",
               resp.tau = 1.5, ttr.median = 3, endpoints = c("os", "response"))
  d <- rOncoEndpoints(nsim = 1, n = 200000, arms = a, seed = 32)
  expect_true(all(is.na(d$pfs_time)) && all(is.na(d$progression)))
  r <- d$response == 1
  expect_true(all(d$os_time[r] > 1.5) && all(d$ttr[r] < d$os_time[r]))
  expect_lt(abs(mean(d$response) - 0.3), 0.005)
  expect_lt(abs(mean(d$os_time > 15) - 0.5), 0.005)
  expect_lt(abs(cor(d$os_time, d$response) - 0.35), 0.01)
  expect_error(rOncoEndpoints(n = 10, arms = a, pfs.cens.hazard = 0.1), "requires PFS")
  other <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4,
                   endpoints = c("pfs", "response"))
  expect_error(rOncoEndpoints(n = c(10, 10), arms = list(a, other)), "same 'endpoints'")
})
