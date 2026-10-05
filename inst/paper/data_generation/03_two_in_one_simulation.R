# Simulation of the 2-in-1 adaptive phase 2/3 design (Figure 4, Table S6,
# Figures S3 and S4)
#
# Run from the package root:
#   source("inst/paper/data_generation/03_two_in_one_simulation.R")
# Writes inst/paper/data/two_in_one_simulation.rds.
#
# Design (based on the hypothetical example of Chen et al. 2018 and the
# multiple-endpoint version of Jin and Zhang 2021):
#   small trial: 120 patients (60 per group) enrolled over 12 months
#   interim analysis: first 90 patients, 3 months after the 90th enrolment;
#     X = pooled Z statistic of the observed response proportions
#   end of the small trial: 70 PFS events in the 120 patients;
#     log-rank statistics Y_PFS and Y_OS
#   expansion: 380 further patients (190 per group) enrolled over 380 / 30
#     months after month 12; end of the large trial at 330 deaths in all 500
#     patients; log-rank statistics Z_PFS and Z_OS
# The statistics are stored for every trial, so that any expansion cutpoint and
# decision rule can be evaluated without new simulations.
#
# Groups: control with PFS median 6, response rate 0.30, Corr(PFS, R) =
# resp.cor, death proportion 0.15, OS median 15, post-progression hazard ratio
# kappa, response timing "ttr" (resp.tau 1.5, ttr.median 2.5).
#   N1 (global null): experimental group identical to the control group
#   N2 (response only): experimental response rate 0.47 with kappa = 1, so that
#      PFS and OS have the same distributions in both groups
#   A  (alternative): experimental PFS hazard ratio 0.55 and response rate 0.47
# Trials: 100,000 per null scenario and 10,000 per alternative scenario.

source(file.path("inst", "paper", "data_generation", "settings.R"))
t_start <- proc.time()
design <- list(n_small = c(60, 60), accrual_small = 12, n_expand = c(190, 190),
               accrual_expand = 380 / 30, m_interim = 90, fu_interim = 3,
               pfs_events_small = 70, os_events_large = 330)
arms_2in1 <- function(type, resp.cor, kappa) {
  switch(type,
         N1 = two_group_arms(resp.cor = resp.cor, pps.hr.resp = kappa, null = TRUE,
                             resp.timing = ttr$resp.timing, resp.tau = ttr$resp.tau,
                             ttr.median = ttr$ttr.median),
         N2 = two_group_arms(resp.cor = resp.cor, pps.hr.resp = kappa, hr.pfs = 1,
                             orr.trt = 0.47, resp.timing = ttr$resp.timing,
                             resp.tau = ttr$resp.tau, ttr.median = ttr$ttr.median),
         A = two_group_arms(resp.cor = resp.cor, pps.hr.resp = kappa, hr.pfs = 0.55,
                            orr.trt = 0.47, resp.timing = ttr$resp.timing,
                            resp.tau = ttr$resp.tau, ttr.median = ttr$ttr.median))
}
scen <- rbind(
  expand.grid(type = "N1", resp_cor = c(0.2, 0.4, 0.6), kappa = c(1, 0.6),
              stringsAsFactors = FALSE),
  expand.grid(type = "N2", resp_cor = c(0.2, 0.4, 0.6), kappa = 1,
              stringsAsFactors = FALSE),
  expand.grid(type = "A", resp_cor = c(0.2, 0.4, 0.6), kappa = c(1, 0.6),
              stringsAsFactors = FALSE))
scen$scenario <- seq_len(nrow(scen))
scen$nsim <- pmax(round(ifelse(scen$type == "A", 10000, 100000) * nsim_scale), 10)

results <- vector("list", nrow(scen))
for (s in seq_len(nrow(scen))) {
  arms <- arms_2in1(scen$type[s], scen$resp_cor[s], scen$kappa[s])
  nsim <- scen$nsim[s]
  nb <- ceiling(nsim / batch_size)
  out <- vector("list", nb)
  for (b in seq_len(nb)) {
    m <- min(batch_size, nsim - (b - 1) * batch_size)
    d <- two_in_one_data(m, arms, design, seed = seed_of(3, s, b))
    st <- two_in_one_statistics(d, design)
    st$sim <- st$sim + (b - 1) * batch_size
    out[[b]] <- st
  }
  ce <- lapply(arms, CorEndpoints)
  results[[s]] <- list(scenario = scen[s, ], stats = do.call(rbind, out),
                       cor_control = ce$Control, cor_treatment = ce$Treatment,
                       os_median = vapply(arms, QuantileEndpoint, numeric(1), endpoint = "os"))
  message("Scenario ", s, " of ", nrow(scen), " done (",
          round((proc.time() - t_start)[["elapsed"]]), " s)")
}

two_in_one_simulation <- list(results = results, scenarios = scen, design = design,
                              ttr = ttr, info = run_info(t_start, pilot))
saveRDS(two_in_one_simulation, file.path(data_dir, "two_in_one_simulation.rds"))
el <- two_in_one_simulation$info$elapsed_sec
message("03_two_in_one_simulation.R: ", round(el, 1), " seconds",
        if (isTRUE(pilot)) paste0("; projected full run: about ", round(el / nsim_scale / 60), " minutes"))
