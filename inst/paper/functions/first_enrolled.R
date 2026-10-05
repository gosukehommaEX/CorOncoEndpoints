# Indicator of the first m patients enrolled in each simulated trial
#
# Arguments
#   sim          trial identifier
#   accrual_time enrolment time
#   m            number of patients
# Value
#   logical vector in the original row order, TRUE for the first m patients
#   (by enrolment time) of the patient's trial
first_enrolled <- function(sim, accrual_time, m) {
  stopifnot(length(sim) == length(accrual_time), length(m) == 1L, m >= 1)
  ord <- order(sim, accrual_time)
  s <- sim[ord]
  n <- length(s)
  first <- c(TRUE, s[-1] != s[-n])
  rank_in_trial <- seq_len(n) - (which(first) - 1L)[cumsum(first)]
  out <- logical(n)
  out[ord] <- rank_in_trial <= m
  out
}
