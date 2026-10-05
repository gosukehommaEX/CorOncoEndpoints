# Table S4: validation of the functions used for the trial examples (summary
# by check; the full list is in
# inst/paper/validation/output/validate_example_functions.md)
#
# Reads inst/paper/validation/output/validate_example_functions.csv. Writes
# inst/paper/output/tables/tableS4_function_validation.tex and
# inst/paper/output/numbers/tableS4_function_validation.csv.

source(file.path("inst", "paper", "figures_and_tables", "settings.R"))
src <- file.path("inst", "paper", "validation", "output", "validate_example_functions.csv")
vf <- utils::read.csv(src, stringsAsFactors = FALSE)
grp <- sub(" \\(trial [0-9]+\\)$", "", vf$check)
grp <- factor(grp, levels = unique(grp))
sm <- do.call(rbind, lapply(split(vf, grp), function(v) {
  data.frame(check = sub(" \\(trial [0-9]+\\)$", "", v$check[1]), n = nrow(v),
             max_diff = max(v$abs_diff), max_tol = max(v$tolerance),
             pass = sum(v$judgment == "PASS"), explained = sum(v$judgment == "EXPLAINED"),
             fail = sum(v$judgment == "FAIL"), stringsAsFactors = FALSE)
}))
esc <- function(x) gsub("_", "\\\\_", gsub("::", "::", x))
judg <- ifelse(sm$fail == 0 & sm$explained == 0, "all PASS",
               paste0(sm$pass, " PASS, ", sm$explained, " EXPLAINED, ", sm$fail, " FAIL"))
body <- paste0(esc(sm$check), " & ", sm$n, " & ",
               formatC(sm$max_diff, format = "e", digits = 1), " & ",
               formatC(sm$max_tol, format = "g", digits = 2), " & ", judg, " \\\\")
write_table_tex("tableS4_function_validation",
                caption = paste0("Validation of the functions used for the trial examples and ",
                                 "of the hazard formulas of Figure 1."),
                label = "tab:functions", align = "p{5.6cm}rrrl",
                header = c("Check & Comparisons & Largest & Tolerance & Judgment \\\\",
                           " & & difference & & \\\\"),
                body = body,
                notes = paste0("Expected values of the hand-made data set were computed ",
                               "independently with loop-based Python code. Simulated data sets ",
                               "were compared with \\texttt{survival::survdiff()}, ",
                               "\\texttt{stats::prop.test()} and the \\texttt{BuyseTest} package."),
                size = "\\scriptsize\\setlength{\\tabcolsep}{3pt}")
write_numbers("tableS4_function_validation",
              key = c(paste0("n_", seq_len(nrow(sm))), paste0("fail_", seq_len(nrow(sm))),
                      "total", "total_fail"),
              value = c(sm$n, sm$fail, sum(sm$n), sum(sm$fail)),
              description = c(paste("comparisons:", sm$check), paste("FAIL:", sm$check),
                              "total comparisons", "total FAIL"),
              source = src)
