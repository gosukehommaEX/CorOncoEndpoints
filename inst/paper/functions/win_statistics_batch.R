# Win statistics for every trial of a batch of data at an analysis cutoff
#
# Arguments
#   d         data frame returned by CutoffData()
#   trt_label label of the experimental group
# Value
#   data frame with one row per trial (ordered by sim): sim and the elements
#   returned by win_statistics()
win_statistics_batch <- function(d, trt_label = "Treatment") {
  idx <- split(seq_len(nrow(d)), d$sim)
  sims <- names(idx)
  res <- do.call(rbind, lapply(idx, function(i) {
    win_statistics(trt = d$group[i] == trt_label,
                   os_tte = d$os_tte[i], os_event = d$os_event[i],
                   pfs_tte = d$pfs_tte[i], pfs_event = d$pfs_event[i],
                   resp = d$response_obs[i])
  }))
  out <- data.frame(sim = as.numeric(sims), res, row.names = NULL)
  out[order(out$sim), , drop = FALSE]
}
