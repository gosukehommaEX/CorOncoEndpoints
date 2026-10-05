# Two-group log-rank statistics for many simulated trials at once
#
# The log-rank statistic is computed separately for each trial (value of
# `sim`) with vectorized operations. Tied times are handled as in
# survival::survdiff(): all subjects with time >= t are at risk at t, and the
# hypergeometric variance includes the factor (n - d) / (n - 1).
#
# Arguments
#   sim   trial identifier (numeric or character)
#   trt   logical or 0/1 vector, TRUE (1) for the experimental group
#   tte   observed time
#   event 1 for an event, 0 for censoring
# Value
#   data frame with one row per trial, ordered by sim, and columns sim,
#   o_minus_e (observed minus expected events in the experimental group),
#   var (variance of o_minus_e), z = -o_minus_e / sqrt(var) (positive values
#   favor the experimental group) and events (total number of events)
logrank_z <- function(sim, trt, tte, event) {
  stopifnot(length(sim) == length(trt), length(sim) == length(tte),
            length(sim) == length(event))
  trt <- as.numeric(trt)
  event <- as.numeric(event)
  ord <- order(sim, -tte)
  s <- sim[ord]
  t <- tte[ord]
  d <- event[ord]
  g <- trt[ord]
  n <- length(t)
  first <- c(TRUE, s[-1] != s[-n])
  start <- which(first)
  k <- length(start)
  trial <- cumsum(first)
  # numbers at risk: counts of subjects with time >= t within the trial
  at_risk <- seq_len(n) - (start - 1L)[trial]
  cs <- cumsum(g)
  at_risk_trt <- cs - c(0, cs)[start][trial]
  # runs of tied times within a trial; the run end holds the risk set
  run_end <- c(s[-1] != s[-n] | t[-1] != t[-n], TRUE)
  run_id <- cumsum(c(TRUE, run_end[-n]))
  d_run <- as.vector(rowsum(d, run_id, reorder = FALSE))
  d1_run <- as.vector(rowsum(d * g, run_id, reorder = FALSE))
  nn <- at_risk[run_end]
  n1 <- at_risk_trt[run_end]
  tr <- trial[run_end]
  keep <- d_run > 0
  e1 <- d_run * n1 / nn
  v <- ifelse(nn > 1, d_run * (n1 / nn) * (1 - n1 / nn) * (nn - d_run) / (nn - 1), 0)
  o_e <- numeric(k)
  vv <- numeric(k)
  ev <- numeric(k)
  if (any(keep)) {
    a <- rowsum(cbind(d1_run - e1, v, d_run)[keep, , drop = FALSE], tr[keep])
    idx <- as.integer(rownames(a))
    o_e[idx] <- a[, 1]
    vv[idx] <- a[, 2]
    ev[idx] <- a[, 3]
  }
  z <- ifelse(vv > 0, -o_e / sqrt(vv), NA_real_)
  data.frame(sim = s[start], o_minus_e = o_e, var = vv, z = z, events = ev,
             stringsAsFactors = FALSE)
}
