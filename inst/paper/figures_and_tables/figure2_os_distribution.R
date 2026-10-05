# Figure 2: OS survival and hazard functions for one group, and the OS hazard
# ratio of two groups
#
# Reads inst/paper/data/model_quantities.rds. Writes
# inst/paper/output/figures/figure2_os_distribution.eps and .pdf and
# inst/paper/output/numbers/figure2_os_distribution.csv.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
src <- "model_quantities.rds"
f2 <- readRDS(file.path(data_dir, src))$fig2
cv <- f2$curves
lev <- c("idm_kappa1", "idm_kappa06", "exponential_os", "pfs")
lab <- expression(paste("OS, illness-death, ", kappa == 1),
                  paste("OS, illness-death, ", kappa == 0.6),
                  "OS, exponential (exp-exp model)",
                  "PFS")
# the exp-exp OS distribution is the exponential distribution with the same
# median, so only the latter is drawn; panel (b) shows the hazards of the same
# four distributions so that (a) and (b) share one legend
d <- cv[cv$model %in% lev, ]
d$model <- factor(d$model, levels = lev)
cols <- pal[c(2, 3, 1, 4)]
ltys <- c("solid", "longdash", "dashed", "dotted")
leg <- theme(legend.position = "bottom", legend.text = element_text(size = 7, hjust = 0),
             legend.key.height = unit(9, "pt"), legend.key.width = unit(18, "pt"),
             legend.key.spacing.y = unit(1, "pt"))
p_a <- ggplot(d, aes(t, surv, colour = model, linetype = model)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = cols, breaks = lev, labels = lab, name = NULL) +
  scale_linetype_manual(values = ltys, breaks = lev, labels = lab, name = NULL) +
  scale_x_continuous(breaks = seq(0, 48, 12)) +
  guides(colour = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2)) +
  labs(x = "Months", y = "Survival probability") +
  theme_paper + leg
p_b <- ggplot(d, aes(t, hazard, colour = model, linetype = model)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = cols, breaks = lev, labels = lab, name = NULL) +
  scale_linetype_manual(values = ltys, breaks = lev, labels = lab, name = NULL) +
  scale_x_continuous(breaks = seq(0, 48, 12)) +
  coord_cartesian(ylim = c(0, NA)) +
  guides(colour = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2)) +
  labs(x = "Months", y = "Hazard (per month)") +
  theme_paper + leg
hr <- f2$hr
lev_k <- c("1", "0.6", "0.3")
hr$kappa <- factor(hr$kappa, levels = lev_k)
lab_k <- expression(kappa == 1, kappa == 0.6, kappa == 0.3)
p_c <- ggplot(hr, aes(t, hr_os, colour = kappa, linetype = kappa)) +
  geom_line(linewidth = 0.6) +
  geom_hline(yintercept = f2$hr_pfs, colour = "grey40", linetype = "dotted") +
  scale_colour_manual(values = pal[c(2, 3, 5)], breaks = lev_k, labels = lab_k, name = NULL) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash"), breaks = lev_k,
                        labels = lab_k, name = NULL) +
  scale_x_continuous(breaks = seq(0, 48, 12)) +
  labs(x = "Months", y = "OS hazard ratio") +
  theme_paper + leg
p <- (((p_a | p_b) + plot_layout(guides = "collect")) | p_c) +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom")
save_figure(p, "figure2_os_distribution", width = fig_width, height = 3.3)

# numbers
tt <- c(0, 6, 12, 24, 36, 48)
get <- function(m, col, x) value_at(cv$t[cv$model == m], cv[[col]][cv$model == m], x)
s_exp <- cv$surv[cv$model == "exponential_os"]
diff_k <- lapply(c("idm_kappa1", "idm_kappa06"), function(m) {
  dlt <- s_exp - cv$surv[cv$model == m]
  c(max = max(dlt), at = cv$t[cv$model == m][which.max(dlt)])
})
hr_keys <- c(1, 6, 12, 24, 36, 48)
hr_vals <- unlist(lapply(c(1, 0.6, 0.3), function(k) {
  value_at(f2$hr$t[f2$hr$kappa == k], f2$hr$hr_os[f2$hr$kappa == k], hr_keys)
}))
med <- unique(f2$hr[, c("kappa", "os_median_ctl", "os_median_trt")])
write_numbers(
  "figure2_os_distribution",
  key = c(paste0("surv_kappa1_t", tt), paste0("surv_kappa06_t", tt), paste0("surv_exp_t", tt),
          paste0("hazard_kappa1_t", tt), paste0("hazard_kappa06_t", tt),
          "maxdiff_exp_minus_kappa1", "maxdiff_at_kappa1",
          "maxdiff_exp_minus_kappa06", "maxdiff_at_kappa06",
          paste0("hr_kappa", rep(c("1", "06", "03"), each = length(hr_keys)), "_t", hr_keys),
          paste0("os_median_ctl_kappa", c("1", "06", "03")),
          paste0("os_median_trt_kappa", c("1", "06", "03"))),
  value = c(get("idm_kappa1", "surv", tt), get("idm_kappa06", "surv", tt),
            get("exponential_os", "surv", tt),
            get("idm_kappa1", "hazard", tt), get("idm_kappa06", "hazard", tt),
            diff_k[[1]], diff_k[[2]], hr_vals,
            med$os_median_ctl[match(c(1, 0.6, 0.3), med$kappa)],
            med$os_median_trt[match(c(1, 0.6, 0.3), med$kappa)]),
  description = c(paste("OS survival, kappa 1, month", tt),
                  paste("OS survival, kappa 0.6, month", tt),
                  paste("exponential OS survival, month", tt),
                  paste("OS hazard, kappa 1, month", tt),
                  paste("OS hazard, kappa 0.6, month", tt),
                  "largest difference exponential minus kappa 1", "month of the largest difference",
                  "largest difference exponential minus kappa 0.6", "month of the largest difference",
                  paste("OS hazard ratio, kappa", rep(c(1, 0.6, 0.3), each = length(hr_keys)),
                        "month", hr_keys),
                  paste("OS median, control, kappa", c(1, 0.6, 0.3)),
                  paste("OS median, experimental, kappa", c(1, 0.6, 0.3))),
  source = src)
