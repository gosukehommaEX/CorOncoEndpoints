# Regression settings: calls whose outputs must not change between versions
#
# The calls use only the interface of version 0.2.0 (accrual over one interval
# or none, scalar or per-group dropout hazards, all three endpoints). Their
# outputs from version 0.2.0 are stored in fixtures/reference_v0.2.0.rds,
# created once with
#   devtools::load_all()
#   source("tests/testthat/helper-regression.R")
#   saveRDS(reg_snapshot(), "tests/testthat/fixtures/reference_v0.2.0.rds",
#           version = 3)
# and the tests of each function compare the current outputs with them. The
# comparison uses a relative tolerance of 1e-8 because the last bits of
# floating-point results can differ between platforms.

reg_arms <- function() {
  ctl <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, death.prop = 0.15,
                 os.median = 15, pps.hr.resp = 0.6, label = "Control")
  list(
    idm_ctl = ctl,
    idm_trt = OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.40,
                      death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6,
                      label = "Treatment"),
    idm_k1 = OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.20, death.prop = 0.30,
                     os.median = 15),
    expexp = OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, death.prop = 0.15,
                     os.model = "expexp", os.median = 15),
    landmark = OncoArm(pfs.median = 6, orr = 0.30, resp.pfs.median = 11,
                       death.prop = 0.15, pps.median = 8, resp.timing = "landmark",
                       resp.tau = 1.5),
    ttr = OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, death.prop = 0.15,
                  os.median = 15, pps.hr.resp = 0.6, resp.timing = "ttr",
                  resp.tau = 1.5, ttr.median = 2.5),
    negative = OncoArm(pfs.hazard = log(2) / 4, orr = 0.45, resp.cor = -0.2,
                       death.prop = 0.05, os.median = 10, pps.hr.resp = 2),
    expexp_ttr = OncoArm(pfs.median = 5, orr = 0.25, resp.cor = 0.3, death.prop = 0.2,
                         os.model = "expexp", os.median = 14, resp.timing = "ttr",
                         resp.tau = 0, ttr.median = 2)
  )
}

# Calibrated parameters of a group (the tabulated post-progression hazard of the
# exp-exp model is kept at every 500th grid point)
reg_arm_fields <- function(arm) {
  keep <- c("lam_p", "orr", "theta", "c_resp", "z_tau", "tau", "death.prop",
            "pps.hr.resp", "ttr_rate", "gam0", "gam1", "lam_o", "c_dec")
  out <- unclass(arm)[intersect(keep, names(arm))]
  if (!is.null(arm$grid)) {
    idx <- seq(1L, length(arm$grid$G), by = 500L)
    out$grid <- list(h = arm$grid$h, G = arm$grid$G[idx], m = arm$grid$m[idx])
  }
  out
}

reg_properties <- function(arms = reg_arms()) {
  tt <- c(0.5, 3, 12, 30)
  lapply(arms, function(a) {
    list(fields = reg_arm_fields(a),
         cor = CorEndpoints(a),
         surv = c(SurvEndpoint(a, tt, endpoint = "pfs"),
                  SurvEndpoint(a, tt, endpoint = "os"),
                  SurvEndpoint(a, tt, endpoint = "os", response = "responders"),
                  SurvEndpoint(a, tt, endpoint = "pfs", response = "nonresponders"),
                  SurvEndpoint(a, tt, endpoint = "os", type = "hazard")),
         quantile = c(QuantileEndpoint(a, endpoint = "os"),
                      QuantileEndpoint(a, endpoint = "pfs", response = "responders")))
  })
}

reg_generate <- function(arms = reg_arms()) {
  list(
    two_groups = rOncoEndpoints(nsim = 2, n = c(15, 20),
                                arms = list(Control = arms$idm_ctl,
                                            Treatment = arms$idm_trt),
                                a.time = c(0, 24), seed = 101),
    expexp = rOncoEndpoints(nsim = 2, n = 20, arms = arms$expexp,
                            a.time = c(0, 24), seed = 102),
    landmark = rOncoEndpoints(nsim = 2, n = 20, arms = arms$landmark, seed = 103),
    ttr_dropout = rOncoEndpoints(nsim = 3, n = c(20, 20),
                                 arms = list(A = arms$ttr, B = arms$expexp_ttr),
                                 a.time = c(0, 12), d.hazard = 0.01, seed = 104),
    per_group_dropout = rOncoEndpoints(nsim = 2, n = c(10, 15),
                                       arms = list(arms$negative, arms$idm_k1),
                                       a.time = c(0, 380 / 30),
                                       d.hazard = c(0.02, 0), seed = 105),
    three_groups = rOncoEndpoints(nsim = 2, n = c(10, 10, 10),
                                  arms = list(arms$idm_ctl, arms$idm_trt,
                                              arms$expexp),
                                  a.time = c(0, 24), seed = 106))
}

# Data at analysis cutoffs (groups without a landmark only)
reg_cutoff <- function(gen = reg_generate()) {
  d1 <- gen$ttr_dropout
  d2 <- gen$two_groups
  list(
    event_time_pfs = EventTime(d1, events = 10, endpoint = "pfs"),
    event_time_os = EventTime(d2, events = 12, endpoint = "os"),
    ttr_scalar = CutoffData(d1, cutoff = 9),
    ttr_events = CutoffData(d1, EventTime(d1, events = 10, endpoint = "pfs")),
    two_groups_os = CutoffData(d2, EventTime(d2, events = 12, endpoint = "os")),
    no_cutoff = CutoffData(gen$per_group_dropout, cutoff = NA))
}

reg_two_groups <- function(arms = reg_arms()) {
  list(Control = arms$idm_ctl, Treatment = arms$idm_trt)
}

reg_expected <- function(arms = reg_arms()) {
  two <- reg_two_groups(arms)
  list(
    os = ExpectedEvents(two, n = c(350, 350), a.time = c(0, 24), endpoint = "os",
                        time = c(12, 24, 36)),
    pfs_dropout = ExpectedEvents(two, n = c(350, 350), a.time = c(0, 24),
                                 d.hazard = 0.01, endpoint = "pfs",
                                 time = c(12, 24, 36)),
    no_accrual = ExpectedEvents(arms$landmark, n = 100, endpoint = "os",
                                time = c(3, 6, 10)))
}

reg_ahr <- function(arms = reg_arms()) {
  list(os = AverageHR(reg_two_groups(arms), n = c(350, 350), a.time = c(0, 24),
                      d.hazard = c(0.01, 0.02), endpoint = "os", time = c(24, 36)))
}

reg_required <- function(arms = reg_arms()) {
  two <- reg_two_groups(arms)
  ee <- list(Control = arms$expexp,
             Treatment = OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.40,
                                 death.prop = 0.15, os.model = "expexp",
                                 os.median = 20))
  list(
    os = RequiredEvents(two, n = c(350, 350), a.time = c(0, 24), endpoint = "os"),
    pfs_unequal = RequiredEvents(two, n = c(233, 467), a.time = c(0, 24),
                                 endpoint = "pfs"),
    expexp_dropout = RequiredEvents(ee, n = c(300, 300), a.time = c(0, 18),
                                    d.hazard = 0.01, endpoint = "os"))
}

reg_bounds <- function() {
  c(CorBoundPFSResponse(orr = 0.3),
    CorBoundPFSResponse(orr = 0.3, pfs.median = 6, resp.tau = 1.5))
}

# All outputs, as stored in fixtures/reference_v0.2.0.rds
reg_snapshot <- function() {
  arms <- reg_arms()
  gen <- reg_generate(arms)
  list(properties = reg_properties(arms), generate = gen,
       cutoff = reg_cutoff(gen), expected = reg_expected(arms),
       ahr = reg_ahr(arms), required = reg_required(arms),
       bounds = reg_bounds(), version = "0.2.0")
}

reg_reference <- function() {
  readRDS(test_path("fixtures", "reference_v0.2.0.rds"))
}

# Compare current outputs with the stored ones: data frames on the stored
# columns, lists element by element
reg_compare <- function(current, reference, label = "output") {
  if (is.data.frame(reference)) {
    expect_true(all(names(reference) %in% names(current)), label = label)
    expect_equal(current[names(reference)], reference, tolerance = 1e-8,
                 label = label)
  } else if (is.list(reference)) {
    for (nm in names(reference)) {
      reg_compare(current[[nm]], reference[[nm]], label = paste0(label, "$", nm))
    }
  } else {
    expect_equal(current, reference, tolerance = 1e-8, label = label)
  }
}
