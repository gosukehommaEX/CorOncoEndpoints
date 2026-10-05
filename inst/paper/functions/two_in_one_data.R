# Simulated patients of a batch of 2-in-1 trials
#
# Patients of the small (phase 2) part are enrolled uniformly over
# [0, design$accrual_small]; the additional patients of the expansion are
# enrolled uniformly over [design$accrual_small, design$accrual_small +
# design$accrual_expand]. Both parts are generated for every trial; whether
# the expansion takes place is decided later from the interim data, so the
# same trials serve every expansion cutpoint. The expansion part uses the seed
# plus 500000.
#
# Arguments
#   nsim   number of trials in the batch
#   arms   list(Control = , Treatment = ) of OncoArm objects
#   design list with n_small (two sample sizes), accrual_small, n_expand (two
#          sample sizes) and accrual_expand
#   seed   integer seed of the small part
# Value
#   data frame as returned by rOncoEndpoints() with the added column part
#   (1 = small part, 2 = expansion)
two_in_one_data <- function(nsim, arms, design, seed) {
  p1 <- rOncoEndpoints(nsim = nsim, n = design$n_small, arms = arms,
                       a.time = c(0, design$accrual_small), seed = seed)
  p2 <- rOncoEndpoints(nsim = nsim, n = design$n_expand, arms = arms,
                       a.time = c(0, design$accrual_expand), seed = seed + 500000L)
  p2$accrual_time <- p2$accrual_time + design$accrual_small
  p2$pfs_calendar_time <- p2$accrual_time + p2$pfs_tte
  p2$os_calendar_time <- p2$accrual_time + p2$os_tte
  p1$part <- 1L
  p2$part <- 2L
  d <- rbind(p1, p2)
  d <- d[order(d$sim, d$accrual_time), , drop = FALSE]
  rownames(d) <- NULL
  d
}
