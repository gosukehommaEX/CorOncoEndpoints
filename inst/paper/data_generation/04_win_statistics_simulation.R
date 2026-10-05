# Simulation of the win ratio and win odds of prioritized OS, PFS and response
# (Table 3 and Table S7)
#
# Run from the package root:
#   source("inst/paper/data_generation/04_win_statistics_simulation.R")
# Writes inst/paper/data/win_statistics_simulation.rds.
#
# Trials: 200 patients per group enrolled uniformly over 24 months, analysed
# at month 30 (fixed calendar time). The sample size is smaller than in the
# design example (350 per group) so that the power of the win ratio test is
# not close to 1 and the differences between the simulated power and the
# formulas remain visible. Endpoints in order of priority: OS, PFS (Gehan
# scores) and response observed by the analysis cutoff (response timing "ttr", resp.tau 1.5,
# ttr.median 2.5). The log-rank statistics of OS and PFS and the pooled Z
# statistic of response at the same cutoff are stored as well.
# Scenarios: Corr(PFS, R) in {0, 0.2, 0.4, 0.6} by kappa in {1, 0.6}, under
# the alternative (experimental PFS hazard ratio 0.7, response rate 0.45) and
# the global null (identical groups); 10,000 trials each.

source(file.path("inst", "paper", "data_generation", "settings.R"))
t_start <- proc.time()
arms_win <- function(resp.cor, kappa, null) {
  two_group_arms(resp.cor = resp.cor, pps.hr.resp = kappa, null = null,
                 hr.pfs = two_group$hr_pfs, orr.trt = two_group$orr_trt,
                 resp.timing = ttr$resp.timing, resp.tau = ttr$resp.tau,
                 ttr.median = ttr$ttr.median)
}
n_win <- c(200, 200)
analysis_time <- 30
scen <- expand.grid(resp_cor = c(0, 0.2, 0.4, 0.6), kappa = c(1, 0.6),
                    hypothesis = c("alternative", "null"), stringsAsFactors = FALSE)
scen$scenario <- seq_len(nrow(scen))
scen$nsim <- max(round(10000 * nsim_scale), 10)

results <- vector("list", nrow(scen))
for (s in seq_len(nrow(scen))) {
  arms <- arms_win(scen$resp_cor[s], scen$kappa[s], scen$hypothesis[s] == "null")
  nsim <- scen$nsim[s]
  nb <- ceiling(nsim / batch_size)
  out <- vector("list", nb)
  for (b in seq_len(nb)) {
    m <- min(batch_size, nsim - (b - 1) * batch_size)
    d <- rOncoEndpoints(nsim = m, n = n_win, arms = arms,
                        a.time = two_group$a_time, seed = seed_of(4, s, b))
    cd <- CutoffData(d, analysis_time)
    trt <- cd$group == "Treatment"
    ws <- win_statistics_batch(cd)
    lr_os <- logrank_z(cd$sim, trt, cd$os_tte, cd$os_event)
    lr_pfs <- logrank_z(cd$sim, trt, cd$pfs_tte, cd$pfs_event)
    oz <- orr_z(cd$sim, trt, cd$response_obs)
    stopifnot(identical(as.numeric(ws$sim), as.numeric(lr_os$sim)),
              identical(as.numeric(ws$sim), as.numeric(oz$sim)))
    ws$z_logrank_os <- lr_os$z
    ws$z_logrank_pfs <- lr_pfs$z
    ws$z_orr <- oz$z
    ws$cutoff <- analysis_time
    ws$deaths <- as.vector(tapply(cd$os_event, cd$sim, sum))
    ws$sim <- ws$sim + (b - 1) * batch_size
    out[[b]] <- ws
  }
  ce <- lapply(arms, CorEndpoints)
  results[[s]] <- list(scenario = scen[s, ], stats = do.call(rbind, out),
                       cor_control = ce$Control, cor_treatment = ce$Treatment)
  message("Scenario ", s, " of ", nrow(scen), " done (",
          round((proc.time() - t_start)[["elapsed"]]), " s)")
}

win_statistics_simulation <- list(results = results, scenarios = scen, n = n_win,
                                  analysis_time = analysis_time, two_group = two_group, ttr = ttr,
                                  info = run_info(t_start, pilot))
saveRDS(win_statistics_simulation, file.path(data_dir, "win_statistics_simulation.rds"))
el <- win_statistics_simulation$info$elapsed_sec
message("04_win_statistics_simulation.R: ", round(el, 1), " seconds",
        if (isTRUE(pilot)) paste0("; projected full run: about ", round(el / nsim_scale / 60), " minutes"))
