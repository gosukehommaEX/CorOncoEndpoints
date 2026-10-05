# Pooled Z statistic for the difference in response proportions, per trial
#
# Z = (p1 - p0) / sqrt(p (1 - p) (1 / n1 + 1 / n0)), where p1 and p0 are the
# proportions of responders in the experimental and control groups and p is
# the pooled proportion. Z is 0 when p is 0 or 1, and NA when a group is empty.
#
# Arguments
#   sim  trial identifier
#   trt  logical or 0/1 vector, TRUE (1) for the experimental group
#   resp 0/1 response indicator
# Value
#   data frame with one row per trial, ordered by sim, and columns sim, n1, n0,
#   x1, x0 (numbers of patients and responders) and z
orr_z <- function(sim, trt, resp) {
  stopifnot(length(sim) == length(trt), length(sim) == length(resp))
  trt <- as.logical(trt)
  resp <- as.numeric(resp)
  sims <- sort(unique(sim))
  f <- factor(sim, levels = sims)
  sum_by <- function(x) {
    v <- tapply(x, f, sum)
    v[is.na(v)] <- 0
    as.vector(v)
  }
  n1 <- sum_by(as.numeric(trt))
  n0 <- sum_by(as.numeric(!trt))
  x1 <- sum_by(resp * trt)
  x0 <- sum_by(resp * !trt)
  p <- (x1 + x0) / (n1 + n0)
  se <- sqrt(p * (1 - p) * (1 / n1 + 1 / n0))
  z <- ifelse(n1 == 0 | n0 == 0, NA_real_,
              ifelse(se > 0, (x1 / n1 - x0 / n0) / se, 0))
  data.frame(sim = sims, n1 = n1, n0 = n0, x1 = x1, x0 = x0, z = z,
             stringsAsFactors = FALSE)
}
