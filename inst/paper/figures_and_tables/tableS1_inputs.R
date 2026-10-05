# Table S1: design inputs, model parameters and how they are determined
#
# No data are read. Writes inst/paper/output/tables/tableS1_inputs.tex.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
body <- c(
  "PFS median (\\texttt{pfs.median}) or hazard & $\\lambda_P$ & $\\lambda_P = \\log 2 / \\mathrm{median}$ \\\\",
  "Response rate (\\texttt{orr}) & $p$ & $c$ solves $P(Z_1 > z_\\tau, Z_2 > c) = p$ \\\\",
  paste0("$\\mathrm{Corr}(\\mathrm{PFS}, R)$ (\\texttt{resp.cor}) or median PFS of responders ",
         "(\\texttt{resp.pfs.median}) & $\\theta$ & One-dimensional root finding \\\\"),
  "Landmark (\\texttt{resp.tau}) & $z_\\tau$ & $z_\\tau = \\Phi^{-1}(1 - e^{-\\lambda_P \\tau})$ \\\\",
  "Median time to response (\\texttt{ttr.median}) & $\\rho$ & $\\rho = \\log 2 / (\\mathrm{median} - \\tau)$ \\\\",
  "Proportion of PFS events that are deaths (\\texttt{death.prop}) & $\\pi$ & Used as given \\\\",
  "\\midrule",
  "\\multicolumn{3}{l}{Illness--death model (\\texttt{os.model = \"idm\"})} \\\\",
  paste0("OS median (\\texttt{os.median}), or median or hazard of post-progression survival ",
         "& $\\gamma_0$ & Root finding of $S_{\\mathrm{OS}}(\\mathrm{median}) = 1/2$, or given \\\\"),
  "Post-progression hazard ratio of responders (\\texttt{pps.hr.resp}) & $\\kappa$ & $\\gamma_1 = \\kappa \\gamma_0$ \\\\",
  "\\midrule",
  "\\multicolumn{3}{l}{Exp--exp model (\\texttt{os.model = \"expexp\"})} \\\\",
  paste0("OS median (\\texttt{os.median}) & $\\lambda_O$, $c$ & $\\lambda_O = \\log 2 / \\mathrm{median}$, ",
         "$c = \\lambda_O / \\pi - \\lambda_P$ \\\\"))
write_table_tex("tableS1_inputs",
                caption = paste0("Design inputs of \\texttt{OncoArm()}, the model parameters they ",
                                 "determine and how the parameters are obtained."),
                label = "tab:inputs", align = "p{5.2cm}p{1.6cm}p{5.4cm}",
                header = "Input (argument) & Parameter & Determination \\\\",
                body = body,
                notes = paste0("$(Z_1, Z_2)$ is standard bivariate normal with correlation $\\theta$; ",
                               "$z_\\tau = -\\infty$ without a landmark. In the exp--exp model the ",
                               "pre-progression death hazard is $\\lambda_O e^{-c t}$, which ",
                               "requires $\\pi \\le \\lambda_O / \\lambda_P$."),
                size = "\\footnotesize")
