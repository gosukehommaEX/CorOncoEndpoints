# Success of 2-in-1 trials under two decision rules
#
# The trial is expanded when X >= cutpoint.
#   rule "chen": success if (X < c and Y_PFS > z_{1 - alpha}) or
#                (X >= c and Z_OS > z_{1 - alpha})          (Chen et al. 2018)
#   rule "g1":   PFS and OS tested at alpha / 2 each, the level of a rejected
#                hypothesis passed to the other (graphical procedure G1 of Jin
#                and Zhang 2021); at least one rejection if
#                max(Y_PFS, Y_OS) > z_{1 - alpha / 2} (small trial) or
#                max(Z_PFS, Z_OS) > z_{1 - alpha / 2} (large trial)
#
# Arguments
#   st       data frame returned by two_in_one_statistics()
#   cutpoint expansion cutpoint for X
#   rule     "chen" or "g1"
#   alpha    one-sided significance level
# Value
#   logical vector, TRUE for successful trials
two_in_one_rules <- function(st, cutpoint, rule = c("chen", "g1"), alpha = 0.025) {
  rule <- match.arg(rule)
  go <- st$x >= cutpoint
  if (rule == "chen") {
    z <- stats::qnorm(1 - alpha)
    ifelse(go, st$z_os > z, st$y_pfs > z)
  } else {
    z <- stats::qnorm(1 - alpha / 2)
    ifelse(go, pmax(st$z_pfs, st$z_os) > z, pmax(st$y_pfs, st$y_os) > z)
  }
}
