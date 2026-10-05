# Control and experimental groups used in the article
#
# The control group has PFS median pfs.median, response rate orr, Corr(PFS, R)
# = resp.cor, a proportion death.prop of PFS events that are deaths, OS median
# os.median and post-progression hazard ratio pps.hr.resp of responders versus
# non-responders (illness-death model). The experimental group has PFS hazard
# ratio hr.pfs and response rate orr.trt, and shares the post-progression
# hazard of non-responders (gam0), pps.hr.resp, death.prop and resp.cor with
# the control group. With null = TRUE the two groups are identical.
#
# Value
#   list(Control = , Treatment = ) of OncoArm objects
two_group_arms <- function(resp.cor = 0.40, pps.hr.resp = 0.6, death.prop = 0.15,
                           hr.pfs = 0.7, orr.trt = 0.45, null = FALSE,
                           pfs.median = 6, orr = 0.30, os.median = 15,
                           resp.timing = "none", resp.tau = 0, ttr.median = NULL) {
  ctl <- OncoArm(pfs.median = pfs.median, orr = orr, resp.cor = resp.cor,
                 death.prop = death.prop, os.median = os.median,
                 pps.hr.resp = pps.hr.resp, resp.timing = resp.timing,
                 resp.tau = resp.tau, ttr.median = ttr.median, label = "Control")
  if (null) {
    trt <- ctl
    trt$label <- "Treatment"
  } else {
    trt <- OncoArm(pfs.median = pfs.median / hr.pfs, orr = orr.trt,
                   resp.cor = resp.cor, death.prop = death.prop,
                   pps.hazard = ctl$gam0, pps.hr.resp = pps.hr.resp,
                   resp.timing = resp.timing, resp.tau = resp.tau,
                   ttr.median = ttr.median, label = "Treatment")
  }
  list(Control = ctl, Treatment = trt)
}
