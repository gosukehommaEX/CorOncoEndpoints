# Test statistics of a batch of 2-in-1 trials
#
# Interim analysis: the first design$m_interim patients of the small part, at
# design$fu_interim months after the enrolment of the last of them; X is the
# pooled Z statistic for the difference in observed response proportions
# (responses observed by the interim cutoff, see CutoffData()).
# End of the small trial: all patients of the small part, when
# design$pfs_events_small PFS events have been observed; Y_PFS and Y_OS are
# log-rank statistics for PFS and OS at that cutoff.
# End of the large trial: all patients, when design$os_events_large deaths have
# been observed; Z_PFS and Z_OS are log-rank statistics at that cutoff.
# All statistics are positive when they favor the experimental group.
#
# Arguments
#   d      data frame returned by two_in_one_data()
#   design list with m_interim, fu_interim, pfs_events_small, os_events_large
#   trt_label label of the experimental group (default "Treatment")
# Value
#   data frame with one row per trial: sim, x, y_pfs, y_os, z_pfs, z_os,
#   t_x, t_y, t_z (calendar times of the three analyses), resp_trt and
#   resp_ctl (responders observed at the interim analysis)
two_in_one_statistics <- function(d, design, trt_label = "Treatment") {
  sims <- sort(unique(d$sim))
  # interim analysis
  small <- d[d$part == 1L, , drop = FALSE]
  di <- small[first_enrolled(small$sim, small$accrual_time, design$m_interim), ,
              drop = FALSE]
  t_x <- as.vector(tapply(di$accrual_time, factor(di$sim, levels = sims), max)) +
    design$fu_interim
  ci <- CutoffData(di, t_x)
  x <- orr_z(ci$sim, ci$group == trt_label, ci$response_obs)
  # end of the small trial
  t_y <- as.vector(EventTime(small, design$pfs_events_small, "pfs"))
  cy <- CutoffData(small, t_y)
  y_pfs <- logrank_z(cy$sim, cy$group == trt_label, cy$pfs_tte, cy$pfs_event)
  y_os <- logrank_z(cy$sim, cy$group == trt_label, cy$os_tte, cy$os_event)
  # end of the large trial
  t_z <- as.vector(EventTime(d, design$os_events_large, "os"))
  cz <- CutoffData(d, t_z)
  z_pfs <- logrank_z(cz$sim, cz$group == trt_label, cz$pfs_tte, cz$pfs_event)
  z_os <- logrank_z(cz$sim, cz$group == trt_label, cz$os_tte, cz$os_event)
  for (r in list(x, y_pfs, y_os, z_pfs, z_os)) {
    if (!identical(as.vector(r$sim), as.vector(sims))) {
      stop("Trial identifiers do not match.", call. = FALSE)
    }
  }
  data.frame(sim = sims, x = x$z, y_pfs = y_pfs$z, y_os = y_os$z,
             z_pfs = z_pfs$z, z_os = z_os$z, t_x = t_x, t_y = t_y, t_z = t_z,
             resp_trt = x$x1, resp_ctl = x$x0, stringsAsFactors = FALSE)
}
