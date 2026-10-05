# Total sample size for the win ratio test (Yu and Ganju, 2022, formula 2)
#
# N = 4 (1 + p_tie) (z_{1 - alpha} + z_{1 - beta})^2 /
#     {3 k (1 - k) (1 - p_tie) log(WR)^2}, rounded up.
#
# Arguments
#   wr     assumed win ratio
#   p_tie  assumed proportion of ties
#   k      allocation proportion (0.5 for 1:1)
#   alpha  one-sided significance level
#   power  target power
win_sample_size <- function(wr, p_tie, k = 0.5, alpha = 0.025, power = 0.9) {
  s2 <- 4 * (1 + p_tie) / (3 * k * (1 - k) * (1 - p_tie))
  ceiling(s2 * (stats::qnorm(1 - alpha) + stats::qnorm(power)) ^ 2 / log(wr) ^ 2)
}
