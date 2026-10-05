# Figure S1: PFS and OS survival functions of responders and non-responders
# (control group of Table 1, illness-death model with kappa = 1 and 0.6)
#
# Reads inst/paper/data/model_quantities.rds. Writes
# inst/paper/output/figures/figureS1_survival_by_response.eps and .pdf.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
mq <- readRDS(file.path(data_dir, "model_quantities.rds"))
d <- mq$figS1
long <- rbind(
  data.frame(t = d$t[d$model == "idm_kappa1"], surv = d$pfs[d$model == "idm_kappa1"],
             response = d$response[d$model == "idm_kappa1"], panel = "PFS"),
  data.frame(t = d$t[d$model == "idm_kappa1"], surv = d$os[d$model == "idm_kappa1"],
             response = d$response[d$model == "idm_kappa1"], panel = "paste('OS, ', kappa == 1)"),
  data.frame(t = d$t[d$model == "idm_kappa06"], surv = d$os[d$model == "idm_kappa06"],
             response = d$response[d$model == "idm_kappa06"], panel = "paste('OS, ', kappa == 0.6)"))
long$panel <- factor(long$panel, levels = c("PFS", "paste('OS, ', kappa == 1)",
                                          "paste('OS, ', kappa == 0.6)"))
long$response <- factor(long$response, levels = c("responders", "nonresponders"),
                        labels = c("Responders", "Non-responders"))
p <- ggplot(long, aes(t, surv, colour = response, linetype = response)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ panel, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = pal[c(2, 3)], name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 48, 12)) +
  labs(x = "Months", y = "Survival probability") +
  theme_paper + theme(legend.position = "bottom")
save_figure(p, "figureS1_survival_by_response", width = fig_width, height = 2.9)
