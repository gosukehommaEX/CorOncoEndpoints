# Figure S4: type I error rate of the 2-in-1 design with PFS and OS tested by
# the graphical procedure G1 of Jin and Zhang (2021), against the expansion
# cutpoint
#
# Reads inst/paper/data/two_in_one_simulation.rds. Writes
# inst/paper/output/figures/figureS4_two_in_one_g1.eps and .pdf and
# inst/paper/output/numbers/figureS4_two_in_one_g1.csv.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
sim <- readRDS(file.path(data_dir, "two_in_one_simulation.rds"))
cuts <- seq(-1, 4, by = 0.05)
ids <- which(sim$scenarios$type %in% c("N1", "N2"))
d <- do.call(rbind, lapply(ids, function(s) {
  st <- sim$results[[s]]$stats
  sc <- sim$scenarios[s, ]
  data.frame(cut = cuts, rate = vapply(cuts, function(cc) mean(two_in_one_rules(st, cc, "g1")),
                                       numeric(1)),
             type = sc$type, resp_cor = sc$resp_cor, kappa = sc$kappa)
}))
d$panel <- factor(ifelse(d$type == "N2", "paste('Response effect only, ', kappa == 1)",
                              paste0("paste('Global null, ', kappa == ", d$kappa, ")")),
                       levels = c("paste('Global null, ', kappa == 1)",
                                  "paste('Global null, ', kappa == 0.6)",
                                  "paste('Response effect only, ', kappa == 1)"))
p <- ggplot(d, aes(cut, rate, colour = factor(resp_cor))) +
  geom_hline(yintercept = 0.025, colour = "grey40", linetype = "dotted") +
  geom_vline(xintercept = stats::qnorm(0.95), colour = "grey40", linetype = "dotted") +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ panel, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = pal[c(2, 3, 4)], name = "Corr(PFS, R)") +
  labs(x = "Expansion cutpoint c for X", y = "Familywise type I error rate") +
  theme_paper + theme(legend.position = "bottom")
save_figure(p, "figureS4_two_in_one_g1", width = fig_width, height = 3.0)

# numbers: largest familywise type I error rate over the cutpoints
grp <- split(d, list(d$type, d$resp_cor, d$kappa), drop = TRUE)
mx <- do.call(rbind, lapply(grp, function(x) {
  data.frame(lab = paste0(x$type[1], "_cor", x$resp_cor[1], "_kappa", x$kappa[1]),
             max_rate = max(x$rate), max_cut = x$cut[which.max(x$rate)],
             stringsAsFactors = FALSE)
}))
write_numbers("figureS4_two_in_one_g1",
              key = c(paste0("max_rate_", mx$lab), paste0("max_cut_", mx$lab)),
              value = c(mx$max_rate, mx$max_cut),
              description = c(paste0("largest familywise type I error rate over c, ", mx$lab),
                              paste0("cutpoint of the largest rate, ", mx$lab)),
              source = "two_in_one_simulation.rds")
