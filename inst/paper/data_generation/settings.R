# Common settings of the data generation scripts for the article
#
# Sourced at the top of each script in inst/paper/data_generation. Run the
# scripts from the package root, after devtools::load_all() or
# library(CorOncoEndpoints). Set pilot <- TRUE before sourcing a script to run
# 1% of the simulated trials; pilot results are written to inst/paper/data/pilot
# and the projected run time of the full script is printed.

if (!exists("OncoArm")) library(CorOncoEndpoints)
source(file.path("inst", "paper", "load_functions.R"))
if (!exists("pilot")) pilot <- FALSE
data_dir <- if (isTRUE(pilot)) {
  file.path("inst", "paper", "data", "pilot")
} else {
  file.path("inst", "paper", "data")
}
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
nsim_scale <- if (isTRUE(pilot)) 0.01 else 1
batch_size <- 1000

# Base case: control group of Sections 2 and 3 of the article
base <- list(pfs_median = 6, orr = 0.30, resp_cor = 0.40, death_prop = 0.15,
             os_median = 15, pps_hr_resp = 0.6)
# Two-group design: experimental group with PFS hazard ratio 0.7 and response
# rate 0.45; 350 patients per group enrolled uniformly over 24 months;
# one-sided alpha 0.025 and power 0.8
two_group <- list(hr_pfs = 0.7, orr_trt = 0.45, n = c(350, 350), a_time = c(0, 24),
                  alpha = 0.025, power = 0.8)
# Response with time to response (used in the trial examples): first tumour
# assessment at 1.5 months, median time to response 2.5 months
ttr <- list(resp.timing = "ttr", resp.tau = 1.5, ttr.median = 2.5)
