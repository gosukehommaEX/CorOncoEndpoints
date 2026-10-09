# Required numbers of events and simulated power of the log-rank test
# (Table 2)
#
# Run from the package root:
#   source("inst/paper/data_generation/02_design_simulation.R")
# Writes inst/paper/data/design_simulation.rds.
#
# Scenarios (350 patients per group unless stated, uniform accrual over 24
# months, no dropout; experimental group with PFS hazard ratio 0.7 and response
# rate 0.45):
#   S1 death proportion 0.15, pps.hr.resp 0.6 (base case)
#   S2 death proportion 0.15, pps.hr.resp 1
#   S3 death proportion 0.30, pps.hr.resp 0.6
#   S4 both groups exp-exp (OS exponential, OS hazard ratio 0.75)
#   S5 death proportion 0.15, pps.hr.resp 0.3
#   S6 as S5 with 2:1 allocation (233 control and 467 experimental patients)
# S5 and S6 were added after the independent review (round 1); the seeds of
# S1 to S4 are unchanged, so their results are the same as before.
# For each scenario the OS analysis is done at the number of events from
# RequiredEvents() (average hazard ratio) and at the number from the Schoenfeld
# formula with the hazard ratio of the two OS medians (exponential assumption),
# and the PFS analysis at the number from RequiredEvents(). A number of events
# that exceeds the sample size cannot be reached and is not simulated. The
# same simulated trials are used for all analyses of a scenario.

source(file.path("inst", "paper", "data_generation", "settings.R"))
t_start <- proc.time()
nsim <- max(round(10000 * nsim_scale), 10)
z2 <- (stats::qnorm(1 - two_group$alpha) + stats::qnorm(two_group$power)) ^ 2
ee_arm <- function(pfs.median, orr, os.median, label) {
  OncoArm(pfs.median = pfs.median, orr = orr, resp.cor = base$resp_cor,
          death.prop = base$death_prop, os.model = "expexp", os.median = os.median,
          label = label)
}
scenarios <- list(
  S1 = list(label = "Death proportion 0.15, kappa 0.6",
            arms = two_group_arms(pps.hr.resp = 0.6)),
  S2 = list(label = "Death proportion 0.15, kappa 1",
            arms = two_group_arms(pps.hr.resp = 1)),
  S3 = list(label = "Death proportion 0.30, kappa 0.6",
            arms = two_group_arms(pps.hr.resp = 0.6, death.prop = 0.30)),
  S4 = list(label = "Exponential-exponential, OS hazard ratio 0.75",
            arms = list(Control = ee_arm(base$pfs_median, base$orr, base$os_median, "Control"),
                        Treatment = ee_arm(base$pfs_median / two_group$hr_pfs,
                                           two_group$orr_trt, base$os_median / 0.75,
                                           "Treatment"))),
  S5 = list(label = "Death proportion 0.15, kappa 0.3",
            arms = two_group_arms(pps.hr.resp = 0.3)),
  S6 = list(label = "Death proportion 0.15, kappa 0.3, allocation 2:1",
            arms = two_group_arms(pps.hr.resp = 0.3), n = c(233, 467)))

results <- list()
for (s in seq_along(scenarios)) {
  sc <- scenarios[[s]]
  arms <- sc$arms
  n_s <- if (is.null(sc$n)) two_group$n else sc$n
  r1 <- n_s[1] / sum(n_s)
  r_os <- RequiredEvents(arms, n = n_s, a.time = two_group$a_time, endpoint = "os",
                         alpha = two_group$alpha, power = two_group$power)
  r_pfs <- RequiredEvents(arms, n = n_s, a.time = two_group$a_time, endpoint = "pfs",
                          alpha = two_group$alpha, power = two_group$power)
  hr_med <- QuantileEndpoint(arms$Control, endpoint = "os") /
    QuantileEndpoint(arms$Treatment, endpoint = "os")
  d_med <- as.integer(ceiling(z2 / (r1 * (1 - r1) * log(hr_med) ^ 2)))
  design <- data.frame(
    analysis = c("os_ahr", "os_median", "pfs_ahr"),
    endpoint = c("os", "os", "pfs"),
    hazard_ratio = c(r_os$ahr, hr_med, r_pfs$ahr),
    events = c(r_os$events, d_med, r_pfs$events),
    expected_time = c(r_os$time, NA, r_pfs$time),
    stringsAsFactors = FALSE)
  if (d_med < sum(n_s)) {
    design$expected_time[2] <- stats::uniroot(function(tc) {
      ExpectedEvents(arms, n = n_s, a.time = two_group$a_time, endpoint = "os",
                     time = tc)$events - d_med
    }, c(1, 500), tol = 1e-8)$root
  }
  design$feasible <- design$events < sum(n_s)
  nb <- ceiling(nsim / batch_size)
  z <- matrix(NA_real_, nsim, nrow(design), dimnames = list(NULL, design$analysis))
  tc <- z
  for (b in seq_len(nb)) {
    idx <- ((b - 1) * batch_size + 1):min(b * batch_size, nsim)
    d <- rOncoEndpoints(nsim = length(idx), n = n_s, arms = arms,
                        a.time = two_group$a_time, seed = seed_of(2, s, b))
    for (k in which(design$feasible)) {
      ep <- design$endpoint[k]
      cutoff <- EventTime(d, design$events[k], ep)
      cd <- CutoffData(d, cutoff)
      lr <- logrank_z(cd$sim, cd$group == "Treatment", cd[[paste0(ep, "_tte")]],
                      cd[[paste0(ep, "_event")]])
      z[idx, k] <- lr$z
      tc[idx, k] <- as.vector(cutoff)
    }
  }
  results[[names(scenarios)[s]]] <- list(label = sc$label, design = design, z = z,
                                         cutoff = tc, nsim = nsim, n = n_s)
  message(names(scenarios)[s], " done (", round((proc.time() - t_start)[["elapsed"]]), " s)")
}

design_simulation <- list(results = results, two_group = two_group, base = base,
                          info = run_info(t_start, pilot))
saveRDS(design_simulation, file.path(data_dir, "design_simulation.rds"))
el <- design_simulation$info$elapsed_sec
message("02_design_simulation.R: ", round(el, 1), " seconds",
        if (isTRUE(pilot)) paste0("; projected full run: about ", round(el / nsim_scale / 60), " minutes"))
