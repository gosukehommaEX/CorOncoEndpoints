# Pairwise Gehan scores for a right-censored time-to-event endpoint
#
# For every pair of an experimental patient i and a control patient j, the
# experimental patient wins when the control patient's event is observed and
# the experimental patient's observed time is longer, and loses when the
# experimental patient's event is observed and the control patient's observed
# time is longer. Other pairs (both censored, the censored time shorter than
# the other's event time, or equal times) are undecided.
#
# Arguments
#   time1, event1 observed times and event indicators of experimental patients
#   time0, event0 the same for control patients
# Value
#   list(win = , loss = ) of logical matrices with length(time1) rows and
#   length(time0) columns
gehan_scores <- function(time1, event1, time0, event0) {
  n1 <- length(time1)
  n0 <- length(time0)
  ev0 <- matrix(event0 == 1, n1, n0, byrow = TRUE)
  ev1 <- matrix(event1 == 1, n1, n0)
  list(win = outer(time1, time0, ">") & ev0,
       loss = outer(time1, time0, "<") & ev1)
}
