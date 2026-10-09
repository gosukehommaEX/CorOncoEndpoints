# Figure 4: type I error rate of the 2-in-1 design (rule of Chen et al. 2018)
# against the expansion cutpoint, simulated and approximated by the bivariate
# normal formula with the simulated correlations of the test statistics
#
# Reads inst/paper/data/two_in_one_simulation.rds. Writes
# inst/paper/output/figures/figure4_two_in_one_type1.eps and .pdf and
# inst/paper/output/numbers/figure4_two_in_one_type1.csv.
#
# Approximation: with X ~ N(mu_X, 1), P(success) = p_Y - P(X >= c, Y > w_Y)
# + P(X >= c, Z > w_Z), where (X, Y) and (X, Z) are bivariate normal with the
# simulated correlations, Y = Y_PFS, Z = Z_OS, mu_X is the simulated mean of X,
# p_Y and p_Z are the simulated rejection rates of Y and Z alone at z_{0.975},
# and w_Y = qnorm(1 - p_Y) and w_Z = qnorm(1 - p_Z), so that the approximation
# equals p_Y as c tends to infinity and p_Z as c tends to minus infinity.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
src <- "two_in_one_simulation.rds"
sim <- readRDS(file.path(data_dir, src))
cuts <- seq(-1, 4, by = 0.05)
bvn <- CorOncoEndpoints:::.bvn_upper
null_ids <- which(sim$scenarios$type %in% c("N1", "N2"))
curves <- do.call(rbind, lapply(null_ids, function(s) {
  st <- sim$results[[s]]$stats
  sc <- sim$scenarios[s, ]
  r_xy <- stats::cor(st$x, st$y_pfs)
  r_xz <- stats::cor(st$x, st$z_os)
  mu <- mean(st$x)
  p_y <- mean(st$y_pfs > z_975)
  p_z <- mean(st$z_os > z_975)
  sim_rate <- vapply(cuts, function(cc) mean(two_in_one_rules(st, cc, "chen")), numeric(1))
  approx <- vapply(cuts, function(cc) {
    p_y - bvn(cc - mu, stats::qnorm(1 - p_y), r_xy) +
      bvn(cc - mu, stats::qnorm(1 - p_z), r_xz)
  }, numeric(1))
  data.frame(scenario = s, type = sc$type, resp_cor = sc$resp_cor, kappa = sc$kappa,
             cut = cuts, simulated = sim_rate, approx = approx, r_xy = r_xy, r_xz = r_xz,
             mu_x = mu, p_y = p_y, p_z = p_z, nsim = nrow(st))
}))
curves$panel <- factor(ifelse(curves$type == "N2", "paste('Response effect only, ', kappa == 1)",
                              paste0("paste('Global null, ', kappa == ", curves$kappa, ")")),
                       levels = c("paste('Global null, ', kappa == 1)",
                                  "paste('Global null, ', kappa == 0.6)",
                                  "paste('Response effect only, ', kappa == 1)"))
long <- rbind(data.frame(curves[, c("panel", "resp_cor", "cut")], rate = curves$simulated,
                         method = "Simulated"),
              data.frame(curves[, c("panel", "resp_cor", "cut")], rate = curves$approx,
                         method = "Bivariate normal approximation"))
long$resp_cor <- factor(long$resp_cor)
long$method <- factor(long$method, levels = c("Simulated", "Bivariate normal approximation"))
lab_r <- setNames(paste0("Corr(PFS, R) = ", levels(long$resp_cor)), levels(long$resp_cor))
p <- ggplot(long, aes(cut, rate, colour = resp_cor, linetype = method)) +
  geom_hline(yintercept = 0.025, colour = "grey40", linetype = "dotted") +
  geom_vline(xintercept = stats::qnorm(0.95), colour = "grey40", linetype = "dotted") +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ panel, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = pal[c(2, 3, 4)], labels = lab_r, name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
  labs(x = "Expansion cutpoint c for X", y = "Type I error rate") +
  theme_paper + theme(legend.position = "bottom", legend.box = "vertical",
                      legend.margin = margin(0, 0, 0, 0))
save_figure(p, "figure4_two_in_one_type1", width = fig_width, height = 3.4)

# numbers: correlations, mean of X and rates at c = 1.645 and the maximum
# over the cutpoints; rejection rates of Y and Z alone; largest excess of the
# simulated rate over the largest of 0.025, p_Y and p_Z; largest absolute
# difference between the simulated and approximated rates
c0 <- stats::qnorm(0.95)
one <- do.call(rbind, lapply(null_ids, function(s) {
  cs <- curves[curves$scenario == s, ]
  st <- sim$results[[s]]$stats
  data.frame(scenario = s, r_xy = cs$r_xy[1], r_xz = cs$r_xz[1], mu_x = cs$mu_x[1],
             rate_c0 = mean(two_in_one_rules(st, c0, "chen")),
             max_sim = max(cs$simulated), max_sim_cut = cs$cut[which.max(cs$simulated)],
             max_approx = max(cs$approx), rate_y = cs$p_y[1], rate_z = cs$p_z[1],
             max_excess = max(cs$simulated - max(0.025, cs$p_y[1], cs$p_z[1])),
             max_abs_diff = max(abs(cs$simulated - cs$approx)), nsim = nrow(st))
}))
lab_s <- with(sim$scenarios[null_ids, ], paste0(type, "_cor", resp_cor, "_kappa", kappa))
vars <- c("r_xy", "r_xz", "mu_x", "rate_c0", "max_sim", "max_sim_cut", "max_approx", "rate_y",
          "rate_z", "max_excess", "max_abs_diff", "nsim")
desc <- c("correlation of X and Y_PFS", "correlation of X and Z_OS", "mean of X",
          "simulated type I error rate at c = 1.645", "largest simulated rate over c",
          "cutpoint of the largest simulated rate", "largest approximated rate over c",
          "rejection rate of Y_PFS alone", "rejection rate of Z_OS alone",
          "largest excess of the simulated rate over max(0.025, Y_PFS alone, Z_OS alone)",
          "largest absolute difference between the simulated and approximated rates",
          "number of simulated trials")
write_numbers("figure4_two_in_one_type1",
              key = as.vector(outer(vars, lab_s, paste, sep = "_")),
              value = as.vector(t(as.matrix(one[, vars]))),
              description = as.vector(outer(desc, lab_s, paste, sep = ", ")),
              source = src)
