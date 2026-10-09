# Table 3: power of the win ratio and win odds tests of prioritized OS, PFS and
# response, simulated and from sample size formulas
#
# Reads inst/paper/data/win_statistics_simulation.rds. Writes
# inst/paper/output/tables/table3_win_statistics.tex and
# inst/paper/output/numbers/table3_win_statistics.csv.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
src <- "win_statistics_simulation.rds"
ws <- readRDS(file.path(data_dir, src))
n_tot <- sum(ws$n)
sm <- do.call(rbind, lapply(ws$results, function(r) {
  data.frame(r$scenario[, c("resp_cor", "kappa", "hypothesis", "scenario")],
             cor_os_resp = r$cor_control[["cor.os.resp"]],
             cor_pfs_os = r$cor_control[["cor.pfs.os"]],
             t(win_summary(r$stats, n_tot, alpha = ws$two_group$alpha)))
}))
alt <- sm[sm$hypothesis == "alternative", ]
alt <- alt[order(-alt$kappa, alt$resp_cor), ]
body <- c()
for (k in unique(alt$kappa)) {
  a <- alt[alt$kappa == k, ]
  body <- c(body, if (length(body) > 0) "\\midrule",
            paste0(ifelse(seq_len(nrow(a)) == 1, paste0("$", k, "$"), ""), " & ",
                   fmt_num(a$resp_cor, 1), " & ", fmt_num(a$cor_os_resp, 2), " & ",
                   fmt_num(a$wr_mean_prob, 3), " & ", fmt_num(a$p_tie, 3), " & ",
                   fmt_num(a$rej_wr, 3), " & ", fmt_num(a$power_wr_formula, 3), " & ",
                   fmt_num(a$power_wr_indep, 3), " & ", fmt_num(a$rej_wo, 3), " & ",
                   fmt_num(a$rej_logrank_os, 3), " & ", fmt_num(a$rej_logrank_pfs, 3), " \\\\"))
}
header <- c(paste0("$\\kappa$ & $\\mathrm{Corr}$ & $\\mathrm{Corr}$ & WR & Ties & ",
                   "\\multicolumn{3}{c}{Power of the WR test} & WO test & \\multicolumn{2}{c}{Log-rank} \\\\"),
            "\\cmidrule(lr){6-8} \\cmidrule(lr){10-11}",
            paste0(" & $(\\mathrm{PFS}, R)$ & $(\\mathrm{OS}, R)$ & & & Simulated & Formula & ",
                   "Independence & Simulated & OS & PFS \\\\"))
nsim <- alt$nsim[1]
caption <- paste0("Win ratio (WR) and win odds (WO) of prioritized OS, PFS and response, and ",
                  "power of their tests and of the log-rank tests in ",
                  format(nsim, big.mark = ","), " simulated trials per scenario.")
notes <- c(paste0(ws$n[1], " patients per group enrolled uniformly over ",
                  ws$two_group$a_time[2], " months; analysis at month ", ws$analysis_time,
                  " (on average ", round(mean(alt$mean_deaths)), " deaths under the alternative); ",
                  "one-sided significance level ", ws$two_group$alpha, ". Pairs are compared on OS, ",
                  "then on PFS (Gehan scores~\\cite{Gehan1965}) and then on response observed by the analysis; ",
                  "WR and Ties are the ratio of the mean proportions of pairs won and lost and ",
                  "the mean proportion of ties. Formula: power from the variance of Yu and Ganju~\\cite{Yu2022} ",
                  "with the mean win, loss and tie proportions. Independence: the same formula ",
                  "with the overall proportions computed from the marginal proportions of each ",
                  "endpoint as if the endpoints were independent (Barnhart et al.~\\cite{Barnhart2025}). ",
                  "$\\mathrm{Corr}(\\mathrm{OS}, R)$ is that of the control group. The Monte Carlo ",
                  "standard error of a power is at most ", fmt_num(sqrt(0.25 / nsim), 3), "."))
write_table_tex("table3_win_statistics", caption = caption, label = "tab:win",
                align = "lrrrrrrrrrr", header = header, body = body, notes = notes,
                size = "\\footnotesize")
vars <- setdiff(names(sm), c("resp_cor", "kappa", "hypothesis", "scenario"))
lab <- paste0(sm$hypothesis, "_cor", sm$resp_cor, "_kappa", sm$kappa)
write_numbers("table3_win_statistics",
              key = as.vector(outer(vars, lab, paste, sep = "_")),
              value = as.vector(t(as.matrix(sm[, vars]))),
              description = as.vector(outer(vars, lab, paste, sep = ", ")),
              source = src)
