# Figure 4: type I error rate of the 2-in-1 design (rule of Chen et al. 2018)
# against the expansion cutpoint, simulated and approximated by the bivariate
# normal formula with the simulated correlations of the test statistics
#
# Reads inst/paper/data/two_in_one_simulation.rds. Writes
# inst/paper/output/figures/figure4_two_in_one_type1.eps and .pdf and
# inst/paper/output/numbers/figure4_two_in_one_type1.csv.
#
# Approximation: with X ~ N(mu_X, 1), P(success) = P(Y > w) - P(X >= c, Y > w)
# + P(X >= c, Z > w), where (X, Y) and (X, Z) are bivariate normal with the
# simulated correlations, w = z_{0.975}, Y = Y_PFS, Z = Z_OS and mu_X is the
# simulated mean of X.

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
  sim_rate <- vapply(cuts, function(cc) mean(two_in_one_rules(st, cc, "chen")), numeric(1))
  approx <- vapply(cuts, function(cc) {
    0.025 - bvn(cc - mu, z_975, r_xy) + bvn(cc - mu, z_975, r_xz)
  }, numeric(1))
  data.frame(scenario = s, type = sc$type, resp_cor = sc$resp_cor, kappa = sc$kappa,
             cut = cuts, simulated = sim_rate, approx = approx, r_xy = r_xy, r_xz = r_xz,
             mu_x = mu, nsim = nrow(st))
}))
curves$panel <- factor(ifelse(curves$type == "N2", "Response effect only, \u03ba = 1",
                              paste0("Global null, \u03ba = ", curves$kappa)),
                       levels = c("Global null, \u03ba = 1", "Global null, \u03ba = 0.6",
                                  "Response effect only, \u03ba = 1"))
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
  facet_wrap(~ panel, nrow = 1) +
  scale_colour_manual(values = pal[c(2, 3, 4)], labels = lab_r, name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
  labs(x = "Expansion cutpoint c for X", y = "Type I error rate") +
  theme_paper + theme(legend.position = "bottom", legend.box = "vertical",
                      legend.margin = margin(0, 0, 0, 0))
save_figure(p, "figure4_two_in_one_type1", width = fig_width, height = 3.4)

# numbers: correlations, mean of X and rates at c = 1.645 and the maximum
# over the cutpoints
c0 <- stats::qnorm(0.95)
one <- do.call(rbind, lapply(null_ids, function(s) {
  cs <- curves[curves$scenario == s, ]
  st <- sim$results[[s]]$stats
  data.frame(scenario = s, r_xy = cs$r_xy[1], r_xz = cs$r_xz[1], mu_x = cs$mu_x[1],
             rate_c0 = mean(two_in_one_rules(st, c0, "chen")),
             max_sim = max(cs$simulated), max_sim_cut = cs$cut[which.max(cs$simulated)],
             max_approx = max(cs$approx), nsim = nrow(st))
}))
lab_s <- with(sim$scenarios[null_ids, ], paste0(type, "_cor", resp_cor, "_kappa", kappa))
vars <- c("r_xy", "r_xz", "mu_x", "rate_c0", "max_sim", "max_sim_cut", "max_approx", "nsim")
desc <- c("correlation of X and Y_PFS", "correlation of X and Z_OS", "mean of X",
          "simulated type I error rate at c = 1.645", "largest simulated rate over c",
          "cutpoint of the largest simulated rate", "largest approximated rate over c",
          "number of simulated trials")
write_numbers("figure4_two_in_one_type1",
              key = as.vector(outer(vars, lab_s, paste, sep = "_")),
              value = as.vector(t(as.matrix(one[, vars]))),
              description = as.vector(outer(desc, lab_s, paste, sep = ", ")),
              source = src)
