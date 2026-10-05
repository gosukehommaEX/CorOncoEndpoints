# Table 1: calibrated parameters and implied quantities of one group under four
# OS models
#
# Reads inst/paper/data/model_quantities.rds. Writes
# inst/paper/output/tables/table1_calibration.tex and
# inst/paper/output/numbers/table1_calibration.csv.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
src <- "model_quantities.rds"
mq <- readRDS(file.path(data_dir, src))
tb <- mq$table1
tb <- tb[match(c("max_indep", "idm_kappa1", "idm_kappa06", "expexp"), tb$model), ]
b <- mq$base
rows <- list(
  c("death_prop", "Proportion of PFS events that are deaths, $\\pi$", 2),
  c("theta", "Copula correlation $\\theta$", 3),
  c("h02_0", "Pre-progression death hazard at time 0", 4),
  c("gam0", "Post-progression hazard, non-responders, $\\gamma_0$", 4),
  c("gam1", "Post-progression hazard, responders, $\\gamma_1$", 4),
  c("pps_median_nonresp", "Median post-progression survival, non-responders", 1),
  c("pps_median_resp", "Median post-progression survival, responders", 1),
  c("pfs_median_nonresp", "Median PFS, non-responders", 1),
  c("pfs_median_resp", "Median PFS, responders", 1),
  c("os_median", "Median OS", 1),
  c("os_median_nonresp", "Median OS, non-responders", 1),
  c("os_median_resp", "Median OS, responders", 1),
  c("cor_pfs_resp", "$\\mathrm{Corr}(\\mathrm{PFS}, R)$", 3),
  c("cor_os_resp", "$\\mathrm{Corr}(\\mathrm{OS}, R)$", 3),
  c("cor_pfs_os", "$\\mathrm{Corr}(\\mathrm{PFS}, \\mathrm{OS})$", 3),
  c("pcor_os_resp", "Partial correlation of OS and $R$ given PFS", 3))
body <- vapply(rows, function(r) {
  paste0(r[2], " & ", paste(fmt_num(tb[[r[1]]], as.integer(r[3])), collapse = " & "), " \\\\")
}, character(1))
header <- c(paste0(" & Maximal & \\multicolumn{2}{c}{Illness--death} & Exp--exp \\\\"),
            "\\cmidrule(lr){3-4}",
            "Quantity & independence & $\\kappa = 1$ & $\\kappa = 0.6$ & \\\\")
caption <- paste0("Calibrated parameters and implied quantities of one group with PFS median ",
                  b$pfs_median, " months, response rate ", fmt_num(b$orr, 2),
                  ", $\\mathrm{Corr}(\\mathrm{PFS}, R) = ", fmt_num(b$resp_cor, 2),
                  "$ and OS median ", b$os_median, " months under four OS models.")
notes <- c(paste0("Hazards are per month and medians are in months. Maximal independence: ",
                  "$\\pi$ equals the ratio of the PFS median to the OS median and the ",
                  "post-progression hazard equals the pre-progression death hazard, so that OS ",
                  "is exponential. Illness--death: $\\pi = ", fmt_num(b$death_prop, 2),
                  "$ and $\\gamma_0$ calibrated to the OS median; $\\kappa = \\gamma_1 / \\gamma_0$. ",
                  "Exp--exp: $\\pi = ", fmt_num(b$death_prop, 2), "$ with PFS and OS both ",
                  "exponential; the post-progression hazard depends on the time since ",
                  "randomization and is the same for responders and non-responders (--)."))
write_table_tex("table1_calibration", caption = caption, label = "tab:calibration",
                align = "lrrrr", header = header, body = body, notes = notes, size = "\\small")
keys <- as.vector(outer(vapply(rows, `[`, "", 1), tb$model, paste, sep = "_"))
vals <- as.vector(as.matrix(tb[, vapply(rows, `[`, "", 1)]))
vals <- as.vector(t(matrix(vals, nrow = nrow(tb))))
write_numbers("table1_calibration", key = keys, value = vals,
              description = as.vector(outer(vapply(rows, `[`, "", 2), tb$model, paste, sep = ", ")),
              source = src)
