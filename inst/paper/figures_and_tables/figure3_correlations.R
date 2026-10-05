# Figure 3: attainable correlations
#   (a) Corr(OS, R) against Corr(PFS, R) for several kappa (OS median kept at
#       15 months by recalibrating gam0)
#   (b) attainable range of Corr(PFS, R) against the response rate
#
# Reads inst/paper/data/model_quantities.rds. Writes
# inst/paper/output/figures/figure3_correlations.eps and .pdf and
# inst/paper/output/numbers/figure3_correlations.csv.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
src <- "model_quantities.rds"
f3 <- readRDS(file.path(data_dir, src))$fig3
d <- f3$cor
kap <- sort(unique(d$kappa))
d$kappa <- factor(d$kappa, levels = kap)
lab_k <- setNames(paste0("\u03ba = ", kap), kap)
p_a <- ggplot(d, aes(resp_cor, cor.os.resp, colour = kappa, linetype = kappa)) +
  geom_abline(intercept = 0, slope = 1, colour = "grey40", linetype = "dotted") +
  geom_line(aes(linewidth = kappa == 1)) +
  scale_linewidth_manual(values = c(0.5, 1.0), guide = "none") +
  scale_colour_manual(values = pal[c(3, 6, 4, 1, 2, 5)], labels = lab_k, name = NULL) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed", "solid", "dotdash", "twodash"),
                        labels = lab_k, name = NULL) +
  coord_cartesian(xlim = c(0, 0.75), ylim = c(-0.2, 0.75)) +
  labs(x = "Corr(PFS, R)", y = "Corr(OS, R)") +
  theme_paper + theme(legend.position = "inside", legend.position.inside = c(0.2, 0.75),
                      legend.text = element_text(size = 7))
b <- f3$bounds
db <- rbind(data.frame(orr = b$orr[b$resp_tau == 0], v = b$upper[b$resp_tau == 0], s = "upper"),
            data.frame(orr = b$orr[b$resp_tau == 0], v = b$lower[b$resp_tau == 0], s = "lower0"),
            data.frame(orr = b$orr[b$resp_tau == 1.5], v = b$lower[b$resp_tau == 1.5],
                       s = "lower15"))
db <- db[!is.na(db$v), ]
lab_b <- c(upper = "Upper bound",
           lower0 = "Lower bound, no landmark",
           lower15 = "Lower bound, landmark 1.5 months (PFS median 6)")
db$s <- factor(db$s, levels = names(lab_b))
p_b <- ggplot(db, aes(orr, v, colour = s, linetype = s)) +
  geom_hline(yintercept = 0, colour = "grey40", linetype = "dotted") +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = pal[c(1, 2, 3)], labels = lab_b, name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed", "longdash"), labels = lab_b, name = NULL) +
  labs(x = "Response rate", y = "Attainable Corr(PFS, R)") +
  theme_paper + theme(legend.position = "inside", legend.position.inside = c(0.55, 0.5),
                      legend.text = element_text(size = 7))
p <- (p_a | p_b) + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
save_figure(p, "figure3_correlations", width = fig_width, height = 3.2)

# numbers
pick <- function(k, r) d$cor.os.resp[d$kappa == k & abs(d$resp_cor - r) < 1e-9]
pfo <- function(k, r) d$cor.pfs.os[d$kappa == k & abs(d$resp_cor - r) < 1e-9]
keys <- c(); vals <- c(); desc <- c()
for (k in kap) for (r in c(0, 0.2, 0.4, 0.6)) {
  keys <- c(keys, sprintf("cor_os_resp_kappa%s_cor%s", k, r), sprintf("cor_pfs_os_kappa%s_cor%s", k, r))
  vals <- c(vals, pick(k, r), pfo(k, r))
  desc <- c(desc, paste("Corr(OS, R), kappa", k, "Corr(PFS, R)", r),
            paste("Corr(PFS, OS), kappa", k, "Corr(PFS, R)", r))
}
bnd <- function(p, tau, w) b[[w]][b$resp_tau == tau & abs(b$orr - p) < 1e-9]
write_numbers("figure3_correlations",
              key = c(keys, "upper_orr0.3", "lower_orr0.3_tau0", "lower_orr0.3_tau1.5"),
              value = c(vals, bnd(0.3, 0, "upper"), bnd(0.3, 0, "lower"), bnd(0.3, 1.5, "lower")),
              description = c(desc, "upper bound of Corr(PFS, R), response rate 0.3",
                              "lower bound, response rate 0.3, no landmark",
                              "lower bound, response rate 0.3, landmark 1.5"),
              source = src)
