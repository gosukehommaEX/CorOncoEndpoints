# Table S7: win statistics of all scenarios, including the type I error rates
# under the global null and the contribution of each endpoint
#
# Reads inst/paper/data/win_statistics_simulation.rds. Writes
# inst/paper/output/tables/tableS7_win_statistics.tex and
# inst/paper/output/numbers/tableS7_win_statistics.csv.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
src <- "win_statistics_simulation.rds"
ws <- readRDS(file.path(data_dir, src))
n_tot <- sum(ws$n)
sm <- do.call(rbind, lapply(ws$results, function(r) {
  data.frame(r$scenario[, c("resp_cor", "kappa", "hypothesis")],
             t(win_summary(r$stats, n_tot, alpha = ws$two_group$alpha)))
}))
sm <- sm[order(sm$hypothesis != "null", -sm$kappa, sm$resp_cor), ]
body <- c()
for (h in c("null", "alternative")) {
  v <- sm[sm$hypothesis == h, ]
  body <- c(body, if (length(body) > 0) "\\midrule",
            paste0(ifelse(seq_len(nrow(v)) == 1, ifelse(h == "null", "Null", "Alternative"), ""),
                   " & ", v$kappa, " & ", fmt_num(v$resp_cor, 1), " & ",
                   fmt_num(v$win_os, 3), "/", fmt_num(v$loss_os, 3), " & ",
                   fmt_num(v$win_pfs, 3), "/", fmt_num(v$loss_pfs, 3), " & ",
                   fmt_num(v$win_resp, 3), "/", fmt_num(v$loss_resp, 3), " & ",
                   fmt_num(v$p_tie, 3), " & ", fmt_num(v$rej_wr, 3), " & ",
                   fmt_num(v$rej_wo, 3), " & ", fmt_num(v$power_wo_formula, 3), " & ",
                   fmt_num(v$power_wo_indep, 3), " & ", fmt_num(v$rej_logrank_os, 3), " & ",
                   fmt_num(v$rej_logrank_pfs, 3), " & ", fmt_num(v$rej_orr, 3), " \\\\"))
}
header <- c(paste0("Hypothesis & $\\kappa$ & $\\mathrm{Corr}$ & \\multicolumn{3}{c}{Won/lost on} & ",
                   "Ties & \\multicolumn{4}{c}{WR and WO tests} & ",
                   "\\multicolumn{3}{c}{Other tests} \\\\"),
            "\\cmidrule(lr){4-6} \\cmidrule(lr){8-11} \\cmidrule(lr){12-14}",
            paste0(" & & $(\\mathrm{PFS}, R)$ & OS & PFS & Response & & WR & WO & WO & WO & ",
                   "Log-rank & Log-rank & Response \\\\"),
            paste0(" & & & & & & & & & formula & indep. & OS & PFS & \\\\"))
write_table_tex("tableS7_win_statistics",
                caption = paste0("Win statistics of prioritized OS, PFS and response: mean ",
                                 "proportions of pairs won and lost on each endpoint, ties, and ",
                                 "rejection rates (one-sided 0.025) under the global null and the ",
                                 "alternative (WR and WO: simulated rejection rates and the power ",
                                 "given by the variance formula)."),
                label = "tab:winfull", align = "lrrrrrrrrrrrrr", header = header, body = body,
                notes = paste0("Settings as in Table 3 of the main text. WO formula and WO independence: power of ",
                               "the win odds test from the variance formula with the mean ",
                               "proportions and with the proportions computed as if the ",
                               "endpoints were independent. ", format(sm$nsim[1], big.mark = ","),
                               " trials per scenario."),
                size = "\\scriptsize\\setlength{\\tabcolsep}{3pt}")
vars <- setdiff(names(sm), c("resp_cor", "kappa", "hypothesis"))
lab <- paste0(sm$hypothesis, "_cor", sm$resp_cor, "_kappa", sm$kappa)
write_numbers("tableS7_win_statistics",
              key = as.vector(outer(vars, lab, paste, sep = "_")),
              value = as.vector(t(as.matrix(sm[, vars]))),
              description = as.vector(outer(vars, lab, paste, sep = ", ")),
              source = src)
