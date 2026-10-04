# Scenarios shared by the tests. Expected values in the tests were computed
# independently in Python (inst/validation/python/reference.py and
# make_expected.py) with scipy (Owen's T function for bivariate normal
# probabilities and adaptive quadrature).

arm_A <- function() {
  OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, death.prop = 0.15,
          os.median = 15, pps.hr.resp = 0.6)
}
arm_B <- function() {
  OncoArm(pfs.median = 6, orr = 0.30, resp.pfs.median = 11, death.prop = 0.15,
          pps.median = 8, resp.timing = "landmark", resp.tau = 1.5)
}
arm_C <- function() {
  OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, death.prop = 0.15,
          os.model = "expexp", os.median = 15, resp.timing = "ttr",
          resp.tau = 1.5, ttr.median = 2.5)
}
arm_D <- function() {
  OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, death.prop = 0.40,
          os.model = "expexp", os.median = 15)
}
# Fleischer et al. (2009) parameterization: TTP hazard l1, death hazard l2,
# post-progression death hazard l3
arm_fleischer <- function(l1, l2, l3) {
  OncoArm(pfs.hazard = l1 + l2, orr = 0.30, resp.cor = 0, death.prop = l2 / (l1 + l2),
          pps.hazard = l3)
}
# Two-group design example: treatment improves PFS (HR 0.7) and response
design_arms <- function() {
  ctl <- arm_A()
  trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.40,
                 death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6)
  list(Control = ctl, Treatment = trt)
}
