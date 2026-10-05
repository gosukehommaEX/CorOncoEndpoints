# Overall win, loss and tie probabilities under independent endpoints
#
# Barnhart et al. (2025): with marginal probabilities w[k] and l[k] that the
# experimental patient wins or loses on endpoint k alone (k in order of
# priority) and t[k] = 1 - w[k] - l[k], independence of the endpoints gives
#   p_win = sum_k w[k] prod_{m < k} t[m],  p_loss = sum_k l[k] prod_{m < k} t[m],
#   p_tie = prod_k t[k].
#
# Arguments
#   w, l numeric vectors of marginal win and loss probabilities
# Value
#   c(p_win = , p_loss = , p_tie = )
win_independence <- function(w, l) {
  stopifnot(length(w) == length(l), all(w + l <= 1 + 1e-12))
  t <- 1 - w - l
  carry <- cumprod(c(1, t))[seq_along(t)]
  c(p_win = sum(w * carry), p_loss = sum(l * carry), p_tie = prod(t))
}
