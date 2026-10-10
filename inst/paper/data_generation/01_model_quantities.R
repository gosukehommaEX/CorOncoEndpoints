# Model quantities (Figures 1 to 3, Table 1, Figures S1 and S2), computed
# without simulation except for the Spearman correlation of PFS and OS
#
# Run from the package root:
#   source("inst/paper/data_generation/01_model_quantities.R")
# Writes inst/paper/data/model_quantities.rds. Random numbers are used only for
# the Spearman correlation of PFS and OS, with the seed seed_of(1, i, 1) for the
# i-th proportion of deaths, seed_of(1, 20 + j, 1) for the j-th shorter OS
# median (12 and 10 months) and seed_of(1, 23, 1) for the example with
# kappa = 0.25 and Corr(PFS, R) = 0.6.

source(file.path("inst", "paper", "data_generation", "settings.R"))
t_start <- proc.time()
lam_p <- log(2) / base$pfs_median
lam_o <- log(2) / base$os_median
pi0 <- base$death_prop
ctl_arm <- function(pps.hr.resp, death.prop = pi0, resp.cor = base$resp_cor, ...) {
  OncoArm(pfs.median = base$pfs_median, orr = base$orr, resp.cor = resp.cor,
          death.prop = death.prop, os.median = base$os_median,
          pps.hr.resp = pps.hr.resp, ...)
}
ee <- OncoArm(pfs.median = base$pfs_median, orr = base$orr, resp.cor = base$resp_cor,
              death.prop = pi0, os.model = "expexp", os.median = base$os_median)

# ---- Figure 1: post-progression hazards implied by exponential PFS and OS ----
t1 <- c(seq(0.02, 1, by = 0.02), seq(1.1, 36, by = 0.1))
h02_const <- pi0 * lam_p
theta_gumbel <- log(pi0) / log(lam_o / lam_p)
h02_ee <- lam_o * exp(-ee$c_dec * t1)
fig1 <- list(
  curves = data.frame(
    t = t1,
    markov_h02 = h02_const,
    markov_h12 = lam_o + (lam_o - h02_const) / (exp((lam_p - lam_o) * t1) - 1),
    expexp_h02 = h02_ee,
    expexp_h12 = (lam_o * exp(-lam_o * t1) - h02_ee * exp(-lam_p * t1)) /
      (exp(-lam_o * t1) - exp(-lam_p * t1)),
    gumbel_h12 = h02_const + (theta_gumbel - 1) * pi0 / t1),
  lam_p = lam_p, lam_o = lam_o, death_prop = pi0, h02_const = h02_const,
  c_dec = ee$c_dec, theta_gumbel = theta_gumbel,
  gam0_idm = ctl_arm(1)$gam0)

# ---- Figure 2: OS survival and hazard; OS hazard ratio of two groups ----
t2 <- seq(0, 48, by = 0.25)
one_group <- list(idm_kappa1 = ctl_arm(1), idm_kappa06 = ctl_arm(0.6), expexp = ee)
fig2_curves <- do.call(rbind, lapply(names(one_group), function(nm) {
  a <- one_group[[nm]]
  data.frame(model = nm, t = t2, surv = SurvEndpoint(a, t2, "os"),
             hazard = SurvEndpoint(a, pmax(t2, 1e-6), "os", type = "hazard"))
}))
fig2_curves <- rbind(fig2_curves,
                     data.frame(model = "exponential_os", t = t2, surv = exp(-lam_o * t2),
                                hazard = lam_o),
                     data.frame(model = "pfs", t = t2, surv = exp(-lam_p * t2),
                                hazard = lam_p))
kappas_hr <- c(1, 0.6, 0.3)
fig2_hr <- do.call(rbind, lapply(kappas_hr, function(k) {
  ar <- two_group_arms(pps.hr.resp = k)
  data.frame(kappa = k, t = t2[t2 > 0],
             hr_os = SurvEndpoint(ar$Treatment, t2[t2 > 0], "os", type = "hazard") /
               SurvEndpoint(ar$Control, t2[t2 > 0], "os", type = "hazard"),
             os_median_ctl = QuantileEndpoint(ar$Control, endpoint = "os"),
             os_median_trt = QuantileEndpoint(ar$Treatment, endpoint = "os"))
}))
fig2 <- list(curves = fig2_curves, hr = fig2_hr, hr_pfs = two_group$hr_pfs)

# ---- Figure 3: attainable correlations ----
kappas <- c(0.25, 0.5, 0.75, 1, 1.5, 2)
cors <- seq(0, 0.75, by = 0.05)
grid3 <- expand.grid(resp_cor = cors, kappa = kappas)
cor3 <- t(vapply(seq_len(nrow(grid3)), function(i) {
  a <- ctl_arm(grid3$kappa[i], resp.cor = grid3$resp_cor[i])
  c(CorEndpoints(a), gam0 = a$gam0, theta = a$theta)
}, numeric(6)))
fig3a <- cbind(grid3, cor3)
orrs <- seq(0.05, 0.95, by = 0.01)
fig3b <- do.call(rbind, lapply(c(0, 1.5), function(tau) {
  b <- t(vapply(orrs, function(p) {
    tryCatch(CorBoundPFSResponse(orr = p, pfs.median = base$pfs_median, resp.tau = tau),
             error = function(e) c(lower = NA_real_, upper = NA_real_))
  }, numeric(2)))
  data.frame(resp_tau = tau, orr = orrs, lower = b[, 1], upper = b[, 2])
}))
fig3 <- list(cor = fig3a, bounds = fig3b)

# ---- Table 1: calibration example ----
tab1_arms <- list(
  max_indep = OncoArm(pfs.median = base$pfs_median, orr = base$orr,
                      resp.cor = base$resp_cor, death.prop = lam_o / lam_p,
                      pps.hazard = lam_o, pps.hr.resp = 1),
  idm_kappa1 = ctl_arm(1),
  idm_kappa06 = ctl_arm(0.6),
  expexp = ee)
tab1 <- do.call(rbind, lapply(names(tab1_arms), function(nm) {
  a <- tab1_arms[[nm]]
  ce <- CorEndpoints(a)
  q <- function(ep, r) QuantileEndpoint(a, 0.5, endpoint = ep, response = r)
  idm <- a$os.model == "idm"
  data.frame(model = nm, death_prop = a$death.prop, theta = a$theta, c_resp = a$c_resp,
             h02_0 = if (idm) a$death.prop * a$lam_p else a$lam_o,
             gam0 = if (idm) a$gam0 else NA_real_,
             gam1 = if (idm) a$gam1 else NA_real_,
             pps_median_nonresp = if (idm) log(2) / a$gam0 else NA_real_,
             pps_median_resp = if (idm) log(2) / a$gam1 else NA_real_,
             pfs_median_resp = q("pfs", "responders"),
             pfs_median_nonresp = q("pfs", "nonresponders"),
             os_median = q("os", "all"),
             os_median_resp = q("os", "responders"),
             os_median_nonresp = q("os", "nonresponders"),
             cor_pfs_resp = ce[["cor.pfs.resp"]], cor_os_resp = ce[["cor.os.resp"]],
             cor_pfs_os = ce[["cor.pfs.os"]], pcor_os_resp = ce[["pcor.os.resp"]],
             stringsAsFactors = FALSE)
}))

# ---- Figure S1: PFS and OS by response status ----
figS1 <- do.call(rbind, lapply(c("idm_kappa1", "idm_kappa06"), function(nm) {
  a <- one_group[[nm]]
  do.call(rbind, lapply(c("responders", "nonresponders"), function(r) {
    data.frame(model = nm, response = r, t = t2,
               pfs = SurvEndpoint(a, t2, "pfs", response = r),
               os = SurvEndpoint(a, t2, "os", response = r))
  }))
}))

# ---- Figure S2: correlations as functions of the proportion of deaths ----
pis <- seq(0.02, 0.40, by = 0.02)
gridS2 <- expand.grid(death_prop = pis, kappa = c(1, 0.6))
corS2 <- t(vapply(seq_len(nrow(gridS2)), function(i) {
  a <- ctl_arm(gridS2$kappa[i], death.prop = gridS2$death_prop[i])
  c(CorEndpoints(a), gam0 = a$gam0)
}, numeric(5)))
figS2 <- cbind(gridS2, corS2)
# Spearman correlation of PFS and OS with kappa = 1, estimated from one simulated
# sample of 10^6 patients for each proportion of deaths
n_spearman <- round(1e6 * nsim_scale)
spearman_pfs_os <- vapply(seq_along(pis), function(i) {
  d <- rOncoEndpoints(nsim = 1, n = n_spearman, arms = ctl_arm(1, death.prop = pis[i]),
                      seed = seed_of(1, i, 1))
  stats::cor(d$pfs_time, d$os_time, method = "spearman")
}, numeric(1))
figS2_spearman <- data.frame(death_prop = pis, kappa = 1, spearman_pfs_os = spearman_pfs_os,
                             n = n_spearman)
# Pearson and Spearman correlations of PFS and OS with kappa = 1 and death
# proportion 0.02 for shorter OS medians (ratios of the medians 0.5 and 0.6);
# with kappa = 1 both depend only on the death proportion and the ratio of the
# medians. Seed seed_of(1, 20 + j, 1) for the j-th OS median.
os_medians_short <- c(12, 10)
death_prop_short <- 0.02
figS2_ratio <- do.call(rbind, lapply(seq_along(os_medians_short), function(j) {
  a <- OncoArm(pfs.median = base$pfs_median, orr = base$orr, resp.cor = base$resp_cor,
               death.prop = death_prop_short, os.median = os_medians_short[j],
               pps.hr.resp = 1)
  d <- rOncoEndpoints(nsim = 1, n = n_spearman, arms = a, seed = seed_of(1, 20 + j, 1))
  data.frame(os_median = os_medians_short[j], death_prop = death_prop_short, kappa = 1,
             median_ratio = base$pfs_median / os_medians_short[j],
             cor_pfs_os = CorEndpoints(a)[["cor.pfs.os"]],
             spearman_pfs_os = stats::cor(d$pfs_time, d$os_time, method = "spearman"),
             n = n_spearman)
}))
# Pearson and Spearman correlations of PFS and OS with kappa = 0.25,
# Corr(PFS, R) = 0.6 and death proportion 0.02 (OS median 15 months): an example
# in which the Spearman correlation exceeds that with kappa = 1. Seed
# seed_of(1, 23, 1).
a_k025 <- ctl_arm(0.25, death.prop = death_prop_short, resp.cor = 0.6)
d_k025 <- rOncoEndpoints(nsim = 1, n = n_spearman, arms = a_k025, seed = seed_of(1, 23, 1))
figS2_kappa <- data.frame(kappa = 0.25, resp_cor = 0.6, death_prop = death_prop_short,
                          cor_pfs_os = CorEndpoints(a_k025)[["cor.pfs.os"]],
                          spearman_pfs_os = stats::cor(d_k025$pfs_time, d_k025$os_time,
                                                       method = "spearman"),
                          n = n_spearman)

model_quantities <- list(fig1 = fig1, fig2 = fig2, fig3 = fig3, table1 = tab1,
                         figS1 = figS1, figS2 = figS2, figS2_spearman = figS2_spearman,
                         figS2_ratio = figS2_ratio, figS2_kappa = figS2_kappa,
                         base = base,
                         two_group = two_group, info = run_info(t_start, pilot))
saveRDS(model_quantities, file.path(data_dir, "model_quantities.rds"))
message("01_model_quantities.R: ", round(model_quantities$info$elapsed_sec, 1), " seconds")
