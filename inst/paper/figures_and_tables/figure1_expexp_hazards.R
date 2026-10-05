# Figure 1: post-progression death hazards implied by exponential PFS and OS
# (Proposition 1)
#
# Reads inst/paper/data/model_quantities.rds. Writes
# inst/paper/output/figures/figure1_expexp_hazards.eps and .pdf and
# inst/paper/output/numbers/figure1_expexp_hazards.csv.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
src <- "model_quantities.rds"
f1 <- readRDS(file.path(data_dir, src))$fig1
cv <- f1$curves

lab_a <- c(markov = "Post-progression, Markov model with constant pre-progression hazard",
           expexp = "Post-progression, exp-exp model",
           expexp_h02 = "Pre-progression death, exp-exp model",
           const_h02 = "Pre-progression death, constant",
           idm = "Post-progression, illness-death model (\u03ba = 1)")
da <- rbind(data.frame(t = cv$t, h = cv$markov_h12, s = "markov"),
            data.frame(t = cv$t, h = cv$expexp_h12, s = "expexp"),
            data.frame(t = cv$t, h = cv$expexp_h02, s = "expexp_h02"),
            data.frame(t = cv$t, h = cv$markov_h02, s = "const_h02"),
            data.frame(t = cv$t, h = f1$gam0_idm, s = "idm"))
da$s <- factor(da$s, levels = names(lab_a))
ylim <- c(0.004, 30)
brk <- c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)
p_a <- ggplot(da, aes(t, h, colour = s, linetype = s)) +
  geom_line(linewidth = 0.6) +
  geom_hline(yintercept = f1$lam_o, colour = "grey40", linetype = "dotted") +
  scale_y_log10(breaks = brk, labels = format(brk, drop0trailing = TRUE)) +
  coord_cartesian(ylim = ylim) +
  scale_colour_manual(values = pal[c(3, 2, 2, 1, 4)], labels = lab_a, name = NULL) +
  scale_linetype_manual(values = c("solid", "solid", "dashed", "dashed", "dotdash"),
                        labels = lab_a, name = NULL) +
  labs(x = "Months since randomization", y = "Hazard (per month)") +
  guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
  theme_paper + theme(legend.position = "bottom", legend.text = element_text(size = 7))

db <- rbind(data.frame(t = cv$t, h = cv$gumbel_h12, s = "gumbel"),
            data.frame(t = cv$t, h = f1$h02_const, s = "const_h02"))
lab_b <- c(gumbel = "Just after progression, Gumbel latent-time model",
           const_h02 = "Pre-progression death")
db$s <- factor(db$s, levels = names(lab_b))
p_b <- ggplot(db, aes(t, h, colour = s, linetype = s)) +
  geom_line(linewidth = 0.6) +
  geom_hline(yintercept = f1$lam_o, colour = "grey40", linetype = "dotted") +
  scale_y_log10(breaks = brk, labels = format(brk, drop0trailing = TRUE)) +
  coord_cartesian(ylim = ylim) +
  scale_colour_manual(values = pal[c(5, 1)], labels = lab_b, name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed"), labels = lab_b, name = NULL) +
  labs(x = "Month of progression", y = "Hazard (per month)") +
  guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
  theme_paper + theme(legend.position = "bottom", legend.text = element_text(size = 7))

p <- (p_a | p_b) + plot_layout(widths = c(1.15, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
save_figure(p, "figure1_expexp_hazards", width = fig_width, height = 4.2)

tt <- c(0.1, 1, 6, 12, 24)
write_numbers(
  "figure1_expexp_hazards",
  key = c("lam_p", "lam_o", "death_prop", "h02_const", "c_dec", "theta_gumbel", "gam0_idm",
          paste0("markov_h12_t", tt), paste0("expexp_h12_t", tt), paste0("expexp_h02_t", tt),
          paste0("gumbel_h12_s", c(0.5, 1, 6))),
  value = c(f1$lam_p, f1$lam_o, f1$death_prop, f1$h02_const, f1$c_dec, f1$theta_gumbel,
            f1$gam0_idm,
            value_at(cv$t, cv$markov_h12, tt), value_at(cv$t, cv$expexp_h12, tt),
            value_at(cv$t, cv$expexp_h02, tt), value_at(cv$t, cv$gumbel_h12, c(0.5, 1, 6))),
  description = c("PFS hazard", "OS hazard", "proportion of PFS events that are deaths",
                  "constant pre-progression death hazard", "decay rate of h02 in exp-exp model",
                  "Gumbel copula parameter giving the death proportion",
                  "post-progression hazard of the illness-death model (kappa 1, OS median 15)",
                  paste("Markov h12 at month", tt), paste("exp-exp h12 at month", tt),
                  paste("exp-exp h02 at month", tt),
                  paste("Gumbel hazard just after progression at month", c(0.5, 1, 6))),
  source = src)
