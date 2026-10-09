# Power of the win ratio or win odds test from win, loss and tie probabilities
#
# Yu and Ganju (2022) and Barnhart et al. (2025):
#   Var(log WR) = 4 (1 + p_tie) / {3 k (1 - k) (1 - p_tie) N},
#   Var(log WO) = 4 (1 + p_tie) (1 - p_tie) / {3 k (1 - k) N},
# where N is the total sample size and k the proportion of patients in one
# group. Power = Phi(log(measure) / sqrt(Var) - z_{1 - alpha}) for a
# one-sided level alpha, where measure > 1 favors the experimental group.
#
# Arguments
#   p_win, p_loss, p_tie probabilities (p_win + p_loss + p_tie = 1)
#   n       total sample size
#   k       allocation proportion (0.5 for 1:1)
#   alpha   one-sided significance level
#   measure "wr" or "wo"
win_power <- function(p_win, p_loss, p_tie, n, k = 0.5, alpha = 0.025,
                      measure = c("wr", "wo")) {
  measure <- match.arg(measure)
  if (measure == "wr") {
    est <- log(p_win / p_loss)
    s2 <- 4 * (1 + p_tie) / (3 * k * (1 - k) * (1 - p_tie))
  } else {
    est <- log((p_win + p_tie / 2) / (p_loss + p_tie / 2))
    s2 <- 4 * (1 + p_tie) * (1 - p_tie) / (3 * k * (1 - k))
  }
  stats::pnorm(sqrt(n / s2) * est - stats::qnorm(1 - alpha))
}
