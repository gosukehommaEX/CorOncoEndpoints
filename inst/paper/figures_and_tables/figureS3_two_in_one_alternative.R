# Figure S3: probability of expansion, probability of success (rule of Chen et
# al. 2018) and expected sample size of the 2-in-1 design under the
# alternative, against the expansion cutpoint
#
# Reads inst/paper/data/two_in_one_simulation.rds. Writes
# inst/paper/output/figures/figureS3_two_in_one_alternative.eps and .pdf.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
sim <- readRDS(file.path(data_dir, "two_in_one_simulation.rds"))
cuts <- seq(-1, 4, by = 0.05)
n_small <- sum(sim$design$n_small)
n_add <- sum(sim$design$n_expand)
ids <- which(sim$scenarios$type == "A")
d <- do.call(rbind, lapply(ids, function(s) {
  st <- sim$results[[s]]$stats
  sc <- sim$scenarios[s, ]
  go <- vapply(cuts, function(cc) mean(st$x >= cc), numeric(1))
  succ <- vapply(cuts, function(cc) mean(two_in_one_rules(st, cc, "chen")), numeric(1))
  out <- rbind(data.frame(cut = cuts, v = go, measure = "Probability of expansion"),
               data.frame(cut = cuts, v = succ, measure = "Probability of success"),
               data.frame(cut = cuts, v = n_small + go * n_add, measure = "Expected sample size"))
  out$resp_cor <- sc$resp_cor
  out$kappa <- sc$kappa
  out
}))
d$measure <- factor(d$measure, levels = c("Probability of expansion", "Probability of success",
                                          "Expected sample size"))
p <- ggplot(d, aes(cut, v, colour = factor(resp_cor), linetype = factor(kappa))) +
  geom_vline(xintercept = stats::qnorm(0.95), colour = "grey40", linetype = "dotted") +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ measure, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = pal[c(2, 3, 4)], name = "Corr(PFS, R)") +
  scale_linetype_manual(values = c("dashed", "solid"), name = "\u03ba") +
  labs(x = "Expansion cutpoint c for X", y = NULL) +
  theme_paper + theme(legend.position = "bottom")
save_figure(p, "figureS3_two_in_one_alternative", width = fig_width, height = 3.0)
