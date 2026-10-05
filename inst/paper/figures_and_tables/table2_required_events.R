# Table 2: required numbers of events from the average hazard ratio and from
# the hazard ratio of the medians, and simulated power of the log-rank test
#
# Reads inst/paper/data/design_simulation.rds. Writes
# inst/paper/output/tables/table2_required_events.tex and
# inst/paper/output/numbers/table2_required_events.csv.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
src <- "design_simulation.rds"
ds <- readRDS(file.path(data_dir, src))
tg <- ds$two_group
lab_sc <- c(S1 = "$\\pi = 0.15$, $\\kappa = 0.6$", S2 = "$\\pi = 0.15$, $\\kappa = 1$",
            S3 = "$\\pi = 0.30$, $\\kappa = 0.6$", S4 = "Exp--exp, OS HR 0.75")
lab_an <- c(os_ahr = "OS, average HR", os_median = "OS, HR of medians",
            pfs_ahr = "PFS")
body <- c()
num <- list()
for (s in names(ds$results)) {
  r <- ds$results[[s]]
  an <- if (s == "S1") c("os_ahr", "os_median", "pfs_ahr") else c("os_ahr", "os_median")
  for (k in an) {
    dz <- r$design[r$design$analysis == k, ]
    feas <- dz$feasible
    pw <- if (feas) mean(r$z[, k] > z_975) else NA_real_
    mcse <- if (feas) sqrt(pw * (1 - pw) / r$nsim) else NA_real_
    mt <- if (feas) mean(r$cutoff[, k]) else NA_real_
    ev <- if (feas) as.character(dz$events) else paste0(dz$events, "\\tnote{a}")
    body <- c(body, paste0(if (k == an[1]) lab_sc[[s]] else "", " & ", lab_an[[k]], " & ",
                           fmt_num(dz$hazard_ratio, 3), " & ", ev, " & ",
                           fmt_num(dz$expected_time, 1), " & ", fmt_num(mt, 1), " & ",
                           if (feas) paste0(fmt_num(pw, 3), " (", fmt_num(mcse, 3), ")") else "--",
                           " \\\\"))
    num[[length(num) + 1L]] <- data.frame(
      key = paste(s, k, c("hr", "events", "expected_time", "mean_time", "power", "mcse"), sep = "_"),
      value = c(dz$hazard_ratio, dz$events, dz$expected_time, mt, pw, mcse),
      description = paste(r$label, lab_an[[k]], c("hazard ratio used", "number of events",
                                                   "expected analysis time (months)",
                                                   "mean simulated analysis time (months)",
                                                   "simulated power", "Monte Carlo SE")))
  }
  if (s != names(ds$results)[length(ds$results)]) body <- c(body, "\\midrule")
}
header <- c(paste0("Scenario & Analysis & HR & Events & \\multicolumn{2}{c}{Analysis time} & ",
                   "Power (MCSE) \\\\"),
            "\\cmidrule(lr){5-6}",
            " & & & & Expected & Simulated & \\\\")
nsim <- ds$results[[1]]$nsim
caption <- paste0("Number of events required by the Schoenfeld formula with the average hazard ",
                  "ratio (average HR) and with the hazard ratio of the two OS medians (HR of ",
                  "medians), and power of the log-rank test in ", format(nsim, big.mark = ","),
                  " simulated trials per scenario.")
notes <- c(paste0(tg$n[1], " patients per group enrolled uniformly over ", tg$a_time[2],
                  " months, one-sided significance level ", tg$alpha, " and target power ",
                  tg$power, ". Control group: PFS median 6 months, response rate 0.30, ",
                  "$\\mathrm{Corr}(\\mathrm{PFS}, R) = 0.40$, OS median 15 months. Experimental ",
                  "group: PFS hazard ratio 0.7 and response rate 0.45, with the same ",
                  "post-progression hazards. Exp--exp: both groups with exponential PFS and OS. ",
                  "Analysis times are in months from the start of enrolment. PFS rows are the ",
                  "same in all scenarios and are shown once. MCSE, Monte Carlo standard error."),
           "[a] More events than patients; not attainable.")
write_table_tex("table2_required_events", caption = caption, label = "tab:events",
                align = "p{2.1cm}p{2.4cm}rrrrr", header = header, body = body, notes = notes,
                size = "\\footnotesize\\setlength{\\tabcolsep}{4pt}")
nm <- do.call(rbind, num)
write_numbers("table2_required_events", key = nm$key, value = nm$value,
              description = nm$description, source = src)
