# Computing time of rOncoEndpoints() (Table S5)
#
# Run from the package root:
#   source("inst/paper/data_generation/05_computing_time.R")
# Writes inst/paper/data/computing_time.rds.
#
# Elapsed time of generating nsim trials of 2 x 250 patients (accrual over 24
# months, dropout hazard 0.01), for the OS models "idm" and "expexp" and the
# response timings "none" and "ttr"; median of 5 repetitions. The machine is
# recorded with the results.

source(file.path("inst", "paper", "data_generation", "settings.R"))
t_start <- proc.time()
reps <- 5
conds <- expand.grid(nsim = c(1000, 10000), os_model = c("idm", "expexp"),
                     resp_timing = c("none", "ttr"), stringsAsFactors = FALSE)
if (isTRUE(pilot)) conds <- conds[conds$nsim == 1000, ]
make_arm <- function(os_model, resp_timing, pfs_median, orr) {
  args <- list(pfs.median = pfs_median, orr = orr, resp.cor = base$resp_cor,
               death.prop = base$death_prop, os.model = os_model,
               os.median = base$os_median / if (orr > base$orr) 0.75 else 1)
  if (os_model == "idm") args$pps.hr.resp <- base$pps_hr_resp
  if (resp_timing == "ttr") args <- c(args, ttr)
  do.call(OncoArm, args)
}
conds$median_sec <- NA_real_
conds$min_sec <- NA_real_
for (i in seq_len(nrow(conds))) {
  arms <- list(Control = make_arm(conds$os_model[i], conds$resp_timing[i],
                                  base$pfs_median, base$orr),
               Treatment = make_arm(conds$os_model[i], conds$resp_timing[i],
                                    base$pfs_median / two_group$hr_pfs, two_group$orr_trt))
  el <- vapply(seq_len(reps), function(r) {
    gc()
    unname(system.time(rOncoEndpoints(nsim = conds$nsim[i], n = c(250, 250), arms = arms,
                                      a.time = c(0, 24), d.hazard = 0.01,
                                      seed = seed_of(5, i, r)))[["elapsed"]])
  }, numeric(1))
  conds$median_sec[i] <- stats::median(el)
  conds$min_sec[i] <- min(el)
}
conds$patients <- conds$nsim * 500
machine <- list(sysname = Sys.info()[["sysname"]], release = Sys.info()[["release"]],
                cpu = tryCatch(if (.Platform$OS.type == "windows") {
                  trimws(utils::tail(system("wmic cpu get name", intern = TRUE), -1)[1])
                } else {
                  sub(".*: ", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1])
                }, error = function(e) NA_character_),
                cores = parallel::detectCores(), r_version = R.version.string,
                dqrng_version = as.character(utils::packageVersion("dqrng")))
computing_time <- list(times = conds, reps = reps, machine = machine,
                       info = run_info(t_start, pilot))
saveRDS(computing_time, file.path(data_dir, "computing_time.rds"))
message("05_computing_time.R: ", round(computing_time$info$elapsed_sec, 1), " seconds")
