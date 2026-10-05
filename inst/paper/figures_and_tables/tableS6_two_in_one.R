# Table S6: operating characteristics of the 2-in-1 design at the expansion
# cutpoint 1.645 under the two decision rules
#
# Reads inst/paper/data/two_in_one_simulation.rds. Writes
# inst/paper/output/tables/tableS6_two_in_one.tex and
# inst/paper/output/numbers/tableS6_two_in_one.csv.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
src <- "two_in_one_simulation.rds"
sim <- readRDS(file.path(data_dir, src))
c0 <- stats::qnorm(0.95)
n_small <- sum(sim$design$n_small)
n_large <- n_small + sum(sim$design$n_expand)
sm <- do.call(rbind, lapply(seq_along(sim$results), function(s) {
  st <- sim$results[[s]]$stats
  sc <- sim$scenarios[s, ]
  go <- mean(st$x >= c0)
  data.frame(scenario = s, type = sc$type, resp_cor = sc$resp_cor, kappa = sc$kappa,
             cor_os_resp = sim$results[[s]]$cor_control[["cor.os.resp"]],
             r_xy = stats::cor(st$x, st$y_pfs), r_xz = stats::cor(st$x, st$z_os),
             r_xy_os = stats::cor(st$x, st$y_os), r_xz_pfs = stats::cor(st$x, st$z_pfs),
             rej_y_pfs = mean(st$y_pfs > z_975), rej_z_os = mean(st$z_os > z_975),
             go = go, ess = n_small + go * (n_large - n_small),
             success_chen = mean(two_in_one_rules(st, c0, "chen")),
             success_g1 = mean(two_in_one_rules(st, c0, "g1")),
             t_x = mean(st$t_x), t_y = mean(st$t_y), t_z = mean(st$t_z), nsim = nrow(st))
}))
lab_type <- c(N1 = "Global null", N2 = "Response only", A = "Alternative")
body <- c()
for (ty in c("N1", "N2", "A")) {
  v <- sm[sm$type == ty, ]
  v <- v[order(-v$kappa, v$resp_cor), ]
  body <- c(body, if (length(body) > 0) "\\midrule",
            paste0(ifelse(seq_len(nrow(v)) == 1, lab_type[[ty]], ""), " & ", v$kappa, " & ",
                   fmt_num(v$resp_cor, 1), " & ", fmt_num(v$cor_os_resp, 2), " & ",
                   fmt_num(v$r_xy, 3), " & ", fmt_num(v$r_xz, 3), " & ",
                   fmt_num(v$rej_y_pfs, 4), " & ", fmt_num(v$rej_z_os, 4), " & ", fmt_num(v$go, 3), " & ",
                   fmt_num(v$ess, 0), " & ", fmt_num(v$success_chen, 4), " & ",
                   fmt_num(v$success_g1, 4), " \\\\"))
}
header <- c(paste0("Scenario & $\\kappa$ & $\\mathrm{Corr}$ & $\\mathrm{Corr}$ & ",
                   "$\\rho_{XY}$ & $\\rho_{XZ}$ & $Y_{\\mathrm{PFS}}$ & $Z_{\\mathrm{OS}}$ & ",
                   "$P(\\mathrm{Go})$ & ESS & ",
                   "\\multicolumn{2}{c}{Probability of success} \\\\"),
            "\\cmidrule(lr){11-12}",
            " & & $(\\mathrm{PFS}, R)$ & $(\\mathrm{OS}, R)$ & & & alone & alone & & & Chen & G1 \\\\")
write_table_tex("tableS6_two_in_one",
                caption = paste0("Operating characteristics of the 2-in-1 design with expansion ",
                                 "when $X \\ge 1.645$."),
                label = "tab:twoinone", align = "lrrrrrrrrrrr", header = header, body = body,
                notes = paste0("$\\rho_{XY}$ and $\\rho_{XZ}$: simulated correlations of $X$ with ",
                               "$Y_{\\mathrm{PFS}}$ and with $Z_{\\mathrm{OS}}$. $Y_{\\mathrm{PFS}}$ alone and ",
                               "$Z_{\\mathrm{OS}}$ alone: probabilities that each statistic exceeds ",
                               "$z_{0.975}$ regardless of the expansion decision. ESS, expected sample ",
                               "size. Chen: success if PFS is significant at the end of the small ",
                               "trial or OS at the end of the large trial (one-sided 0.025). G1: PFS ",
                               "and OS tested at 0.0125 each with the level passed on after a ",
                               "rejection. Global null: identical groups. Response only: response ",
                               "rate 0.47 versus 0.30 with identical PFS and OS distributions. ",
                               "Alternative: PFS hazard ratio 0.55 and response rate 0.47. ",
                               format(sm$nsim[sm$type == "N1"][1], big.mark = ","), " trials per null ",
                               "scenario and ", format(sm$nsim[sm$type == "A"][1], big.mark = ","),
                               " per alternative scenario."),
                size = "\\footnotesize")
vars <- setdiff(names(sm), c("scenario", "type", "resp_cor", "kappa"))
lab <- paste0(sm$type, "_cor", sm$resp_cor, "_kappa", sm$kappa)
write_numbers("tableS6_two_in_one",
              key = as.vector(outer(vars, lab, paste, sep = "_")),
              value = as.vector(t(as.matrix(sm[, vars]))),
              description = as.vector(outer(vars, lab, paste, sep = ", ")),
              source = src)
