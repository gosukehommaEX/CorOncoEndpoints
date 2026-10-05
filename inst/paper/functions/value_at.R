# Value of a tabulated curve at the grid point nearest to x
#
# Arguments
#   t grid, y values on the grid, x requested points
value_at <- function(t, y, x) {
  vapply(x, function(xx) y[which.min(abs(t - xx))], numeric(1))
}
