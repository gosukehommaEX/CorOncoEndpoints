# Figure S2: Corr(PFS, OS) and Corr(OS, R) against the proportion of PFS events
# that are deaths (OS median kept at 15 months; Corr(PFS, R) = 0.4). The
# Spearman correlation of PFS and OS (kappa = 1), and the Pearson and Spearman
# correlations for OS medians of 12 and 10 months (kappa = 1, death proportion
# 0.02) and for kappa = 0.25 with Corr(PFS, R) = 0.6 (death proportion 0.02,
# OS median 15 months), are written to the numbers only.
#
# Reads inst/paper/data/model_quantities.rds. Writes
# inst/paper/output/figures/figureS2_death_proportion.eps and .pdf and
# inst/paper/output/numbers/figureS2_death_proportion.csv.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
src <- "model_quantities.rds"
d <- readRDS(file.path(data_dir, src))$figS2
ds <- readRDS(file.path(data_dir, src))$figS2_spearman
dr <- readRDS(file.path(data_dir, src))$figS2_ratio
dk <- readRDS(file.path(data_dir, src))$figS2_kappa
long <- rbind(data.frame(pi = d$death_prop, kappa = d$kappa, v = d$cor.pfs.os,
                         panel = "Corr(PFS, OS)"),
              data.frame(pi = d$death_prop, kappa = d$kappa, v = d$cor.os.resp,
                         panel = "Corr(OS, R)"))
long$panel <- factor(long$panel, levels = c("Corr(PFS, OS)", "Corr(OS, R)"))
long$kappa <- factor(long$kappa, levels = c(1, 0.6))
p <- ggplot(long, aes(pi, v, colour = kappa, linetype = kappa)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = pal[c(2, 3)], labels = expression(kappa == 1, kappa == 0.6),
                      name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed"),
                        labels = expression(kappa == 1, kappa == 0.6), name = NULL) +
  labs(x = "Proportion of PFS events that are deaths", y = "Correlation") +
  theme_paper + theme(legend.position = "bottom")
save_figure(p, "figureS2_death_proportion", width = fig_width, height = 2.9)
pk <- function(k, p0, col) d[[col]][d$kappa == k & abs(d$death_prop - p0) < 1e-9]
sp <- function(p0) ds$spearman_pfs_os[abs(ds$death_prop - p0) < 1e-9]
so <- function(m, col) dr[[col]][dr$os_median == m]
write_numbers("figureS2_death_proportion",
              key = c("cor_pfs_os_kappa1_pi0.02", "cor_pfs_os_kappa1_pi0.4",
                      "cor_pfs_os_kappa06_pi0.02", "cor_pfs_os_kappa06_pi0.4",
                      "spearman_pfs_os_kappa1_pi0.02", "spearman_pfs_os_kappa1_pi0.4",
                      "spearman_pfs_os_kappa1_max", "spearman_n",
                      "spearman_pfs_os_kappa1_pi0.02_os12", "spearman_pfs_os_kappa1_pi0.02_os10",
                      "cor_pfs_os_kappa1_pi0.02_os12", "cor_pfs_os_kappa1_pi0.02_os10",
                      "spearman_pfs_os_kappa0.25_cor0.6_pi0.02", "cor_pfs_os_kappa0.25_cor0.6_pi0.02"),
              value = c(pk(1, 0.02, "cor.pfs.os"), pk(1, 0.4, "cor.pfs.os"),
                        pk(0.6, 0.02, "cor.pfs.os"), pk(0.6, 0.4, "cor.pfs.os"),
                        sp(0.02), sp(0.4), max(ds$spearman_pfs_os), ds$n[1],
                        so(12, "spearman_pfs_os"), so(10, "spearman_pfs_os"),
                        so(12, "cor_pfs_os"), so(10, "cor_pfs_os"),
                        dk$spearman_pfs_os, dk$cor_pfs_os),
              description = c(paste("Corr(PFS, OS), kappa", c(1, 1, 0.6, 0.6),
                                    "death proportion", c(0.02, 0.4, 0.02, 0.4)),
                              paste("Spearman correlation of PFS and OS (simulated), kappa 1",
                                    "death proportion", c(0.02, 0.4)),
                              paste("largest Spearman correlation of PFS and OS (simulated),",
                                    "kappa 1, death proportion 0.02 to 0.40"),
                              "number of simulated patients for each Spearman correlation",
                              paste("Spearman correlation of PFS and OS (simulated), kappa 1,",
                                    "death proportion 0.02, OS median", c(12, 10)),
                              paste("Corr(PFS, OS), kappa 1, death proportion 0.02, OS median",
                                    c(12, 10)),
                              paste("Spearman correlation of PFS and OS (simulated), kappa 0.25,",
                                    "Corr(PFS, R) 0.6, death proportion 0.02, OS median 15"),
                              paste("Corr(PFS, OS), kappa 0.25, Corr(PFS, R) 0.6,",
                                    "death proportion 0.02, OS median 15")),
              source = src)
