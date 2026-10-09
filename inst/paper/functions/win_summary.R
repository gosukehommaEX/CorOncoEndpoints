# Summary of the simulated win statistics of one scenario
#
# Rejection rates (one-sided level alpha) of the win ratio, win odds, log-rank
# (OS, PFS) and response tests, mean win/loss/tie proportions overall and by
# endpoint, and the power given by the formulas of Yu and Ganju (2022) and
# Barnhart et al. (2025) with (i) the mean overall probabilities and (ii) the
# overall probabilities implied by the mean marginal probabilities under
# independence of the endpoints. It also returns the standard deviation of
# log(WR) over the trials, the mean of its estimated standard errors and the
# standard deviation sqrt(Var / N) given by the variance formula of Yu and
# Ganju (2022) with the mean probabilities and 1:1 allocation.
#
# Arguments
#   st    data frame of per-trial statistics (04_win_statistics_simulation.R)
#   n     total sample size
#   alpha one-sided significance level
# Value
#   named numeric vector
win_summary <- function(st, n, alpha = 0.025) {
  z <- stats::qnorm(1 - alpha)
  rate <- function(x) mean(x > z, na.rm = TRUE)
  m <- colMeans(st[, c("p_win", "p_loss", "p_tie", "win_os", "loss_os", "win_pfs",
                       "loss_pfs", "win_resp", "loss_resp", "mwin_os", "mloss_os",
                       "mwin_pfs", "mloss_pfs", "mwin_resp", "mloss_resp")])
  ind <- win_independence(m[c("mwin_os", "mwin_pfs", "mwin_resp")],
                          m[c("mloss_os", "mloss_pfs", "mloss_resp")])
  c(rej_wr = rate(st$z_wr), rej_wo = rate(st$z_wo),
    rej_logrank_os = rate(st$z_logrank_os), rej_logrank_pfs = rate(st$z_logrank_pfs),
    rej_orr = rate(st$z_orr), m,
    wr_mean_prob = unname(m["p_win"] / m["p_loss"]),
    wo_mean_prob = unname((m["p_win"] + m["p_tie"] / 2) / (m["p_loss"] + m["p_tie"] / 2)),
    power_wr_formula = win_power(m[["p_win"]], m[["p_loss"]], m[["p_tie"]], n, alpha = alpha,
                                 measure = "wr"),
    power_wo_formula = win_power(m[["p_win"]], m[["p_loss"]], m[["p_tie"]], n, alpha = alpha,
                                 measure = "wo"),
    p_win_indep = ind[["p_win"]], p_loss_indep = ind[["p_loss"]], p_tie_indep = ind[["p_tie"]],
    power_wr_indep = win_power(ind[["p_win"]], ind[["p_loss"]], ind[["p_tie"]], n,
                               alpha = alpha, measure = "wr"),
    power_wo_indep = win_power(ind[["p_win"]], ind[["p_loss"]], ind[["p_tie"]], n,
                               alpha = alpha, measure = "wo"),
    sd_log_wr = stats::sd(log(st$wr)), mean_se_log_wr = mean(st$se_log_wr),
    sigma_wr_formula = sqrt(4 * (1 + m[["p_tie"]]) / (3 * 0.25 * (1 - m[["p_tie"]]) * n)),
    nsim = nrow(st), mean_cutoff = mean(st$cutoff),
    mean_deaths = if (is.null(st$deaths)) NA_real_ else mean(st$deaths))
}
