# Table S5: computing time of rOncoEndpoints()
#
# Reads inst/paper/data/computing_time.rds. Writes
# inst/paper/output/tables/tableS5_computing_time.tex and
# inst/paper/output/numbers/tableS5_computing_time.csv.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
src <- "computing_time.rds"
ct <- readRDS(file.path(data_dir, src))
tm <- ct$times[order(ct$times$os_model != "idm", ct$times$resp_timing, ct$times$nsim), ]
body <- paste0("\\texttt{", tm$os_model, "} & \\texttt{", tm$resp_timing, "} & ",
               format(tm$nsim, big.mark = ","), " & ", format(tm$patients, big.mark = ","), " & ",
               fmt_num(tm$median_sec, 2), " & ",
               fmt_num(tm$patients / tm$median_sec / 1e6, 1), " \\\\")
mc <- ct$machine
write_table_tex("tableS5_computing_time",
                caption = paste0("Elapsed time of \\texttt{rOncoEndpoints()} for two groups of ",
                                 "250 patients (median of ", ct$reps, " repetitions)."),
                label = "tab:time", align = "llrrrr",
                header = c("OS model & Response & Trials & Patients & Seconds & Million patients \\\\",
                           " & timing & & & & per second \\\\"),
                body = body,
                notes = paste0("Accrual over 24 months and dropout hazard 0.01 per month. ",
                               "Machine: ", gsub("_", "\\\\_", mc$cpu), ", ", mc$cores,
                               " logical cores, ", mc$sysname, " ", mc$release, ", ", mc$r_version,
                               ", dqrng ", mc$dqrng_version, ". One core is used."),
                size = "\\small")
write_numbers("tableS5_computing_time",
              key = paste0("sec_", tm$os_model, "_", tm$resp_timing, "_", tm$nsim),
              value = tm$median_sec,
              description = paste("median elapsed seconds:", tm$os_model, tm$resp_timing, tm$nsim, "trials"),
              source = src)
