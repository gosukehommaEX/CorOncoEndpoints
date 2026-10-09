# Reproduction of published numbers with CorOncoEndpoints
#
# Run from the package root, after devtools::load_all() or
# library(CorOncoEndpoints):
#   source("inst/reproduce/reproduce_published.R")
# Writes inst/reproduce/output/reproduce_results.csv and
# inst/reproduce/output/reproduce_summary.md (generated, do not edit by hand).
#
# Judgments: PASS (agrees within the stated tolerance), EXPLAINED (differs, with
# a documented reason), FAIL (differs without a reason), INFO (not reproducible
# from the published information; reason given).

if (!exists("OncoArm")) library(CorOncoEndpoints)
out_dir <- file.path("inst", "reproduce", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rows <- list()
add <- function(source, quantity, published, computed, tol, reason = "") {
  judg <- if (is.na(published) || is.na(computed)) {
    "INFO"
  } else if (abs(computed - published) <= tol) {
    "PASS"
  } else if (nzchar(reason)) {
    "EXPLAINED"
  } else {
    "FAIL"
  }
  if (judg == "PASS") reason <- ""
  rows[[length(rows) + 1L]] <<- data.frame(
    source = source, quantity = quantity, published = published,
    computed = computed, tolerance = tol, judgment = judg, reason = reason,
    stringsAsFactors = FALSE)
}

# Fleischer et al. (2009) parameterization: TTP hazard l1, death hazard l2,
# post-progression death hazard l3; response is not involved (resp.cor = 0)
fl_arm <- function(l1, l2, l3) {
  OncoArm(pfs.hazard = l1 + l2, orr = 0.3, resp.cor = 0, death.prop = l2 / (l1 + l2),
          pps.hazard = l3)
}

# ---- Fleischer et al. (2009), Statistics in Medicine 28:2669-2686 ----
src <- "Fleischer et al. (2009)"
add(src, "Corr(PFS, OS), lambda = (2, 0.5, 1.5), Section 3", 0.52,
    CorEndpoints(fl_arm(2, 0.5, 1.5))[["cor.pfs.os"]], 0.005)
ex1 <- fl_arm(0.284, 0.075, 0.128)
add(src, "Example 1: Corr(PFS, OS)", 0.34, CorEndpoints(ex1)[["cor.pfs.os"]], 0.005)
d1 <- rOncoEndpoints(nsim = 1, n = 1e6, arms = ex1, seed = 2009)
add(src, "Example 1: proportion dying without progression (simulated, n = 1e6)", 0.209,
    mean(d1$progression == 0), 0.002)
ex3c <- fl_arm(0.342, 0.0654, 0.0797)
ex3t <- fl_arm(0.139, 0.0648, 0.0876)
ex3s <- fl_arm(0.342, 0.0654, 0.0876)
add(src, "Example 3: median PFS, control", 1.7, QuantileEndpoint(ex3c, endpoint = "pfs"), 0.05)
add(src, "Example 3: median PFS, treatment", 3.4, QuantileEndpoint(ex3t, endpoint = "pfs"), 0.05)
add(src, "Example 3: median OS, control", 9.2, QuantileEndpoint(ex3c, endpoint = "os"), 0.05)
add(src, "Example 3: median OS, treatment", 9.3, QuantileEndpoint(ex3t, endpoint = "os"), 0.05)
add(src, "Example 3: 25th percentile of OS, control", 4.0,
    QuantileEndpoint(ex3c, prob = 0.25, endpoint = "os"), 0.05)
add(src, "Example 3: 25th percentile of OS, treatment", 4.1,
    QuantileEndpoint(ex3t, prob = 0.25, endpoint = "os"), 0.05)
add(src, "Example 3: median OS, control with lambda3 = 0.0876", 8.6,
    QuantileEndpoint(ex3s, endpoint = "os"), 0.05)

# Example 2: 1300 patients, 45 per month, medians PFS/OS 5.1/10.6 (treatment)
# and 4/9 (control); probability of death before progression 1/3
acc_end <- 1300 / 45
ex2 <- function(pfs, os, model) {
  if (model == "exp") {
    OncoArm(pfs.median = pfs, orr = 0.3, resp.cor = 0, death.prop = pfs / os,
            os.model = "expexp", os.median = os)
  } else {
    OncoArm(pfs.median = pfs, orr = 0.3, resp.cor = 0, death.prop = 1 / 3,
            os.median = os)
  }
}
arms_exp <- list(Control = ex2(4, 9, "exp"), Treatment = ex2(5.1, 10.6, "exp"))
arms_mod <- list(Control = ex2(4, 9, "idm"), Treatment = ex2(5.1, 10.6, "idm"))
e_exp <- ExpectedEvents(arms_exp, n = c(650, 650), a.time = c(0, acc_end),
                        endpoint = "os", time = 47)$events
e_mod <- ExpectedEvents(arms_mod, n = c(650, 650), a.time = c(0, acc_end),
                        endpoint = "os", time = 47)$events
e_mod29 <- ExpectedEvents(arms_mod, n = c(652.5, 652.5), a.time = c(0, 29),
                          endpoint = "os", time = 47)$events
add(src, "Example 2: expected OS events at 47 months, exponential OS", 1147, e_exp, 1)
add(src, "Example 2: expected OS events at 47 months, model with death proportion 1/3",
    1201, e_mod, 1,
    reason = paste0("The paper gives the accrual as 45 patients per month with a total ",
                    "accrual time of approx. 29 months. With exactly 1300 / 45 months the ",
                    "expected number is ", sprintf("%.1f", e_mod), "; with 29 months and ",
                    "1305 patients it is ", sprintf("%.1f", e_mod29), ". The published ",
                    "1201 lies between, and the exponential case above differs in the ",
                    "same direction."))
add(src, "Example 2: power of the log-rank test at 1201 events (88.7%)", 0.887, NA, NA,
    reason = paste0("The paper does not state how the power was computed. The Schoenfeld ",
                    "formula with the hazard ratio 9/10.6 gives ",
                    sprintf("%.3f", pnorm(sqrt(1201) * abs(log(9 / 10.6)) / 2 -
                                            qnorm(0.975))),
                    " at 1201 events."))

# ---- TrialSimulator documentation, example of solveThreeStateModel() ----
# The published values are comments of the example, computed from 1e6 patients
# generated by CorrelatedPfsAndOs3() (Monte Carlo values).
src <- "TrialSimulator (solveThreeStateModel example)"
ts <- fl_arm(0.1, 0.05, 0.12)
add(src, "Corr(PFS, OS), h01 = 0.1, h02 = 0.05, h12 = 0.12 (Monte Carlo, n = 1e6)", 0.65,
    CorEndpoints(ts)[["cor.pfs.os"]], 0.005)
add(src, "median PFS (Monte Carlo, n = 1e6)", 4.62, QuantileEndpoint(ts, endpoint = "pfs"), 0.005)
add(src, "median OS (Monte Carlo, n = 1e6)", 9.61, QuantileEndpoint(ts, endpoint = "os"), 0.005)

# ---- Erdmann et al. (2025), Biometrical Journal 67:e70017, Table 1 ----
src <- "Erdmann et al. (2025)"
erd <- list(c(0.06, 0.30, 0.10, 0.40, 0.720), c(0.30, 0.28, 0.50, 0.30, 0.725),
            c(0.140, 0.112, 0.180, 0.150, 0.764), c(0.18, 0.06, 0.23, 0.07, 0.800))
for (i in seq_along(erd)) {
  v <- erd[[i]]
  trt <- OncoArm(pfs.hazard = v[1] + v[2], orr = 0.3, resp.cor = 0,
                 death.prop = v[2] / (v[1] + v[2]), pps.hazard = 0.3)
  ctl <- OncoArm(pfs.hazard = v[3] + v[4], orr = 0.3, resp.cor = 0,
                 death.prop = v[4] / (v[3] + v[4]), pps.hazard = 0.3)
  hr <- SurvEndpoint(trt, 1, "pfs", type = "hazard") / SurvEndpoint(ctl, 1, "pfs", type = "hazard")
  add(src, paste0("Table 1, scenario ", i, ": PFS hazard ratio"), v[5], hr, 0.0005)
}

res <- do.call(rbind, rows)
utils::write.csv(res, file.path(out_dir, "reproduce_results.csv"), row.names = FALSE)
tab <- table(factor(res$judgment, levels = c("PASS", "EXPLAINED", "FAIL", "INFO")))
md <- c("# Reproduction of published numbers",
        "",
        "Generated by `inst/reproduce/reproduce_published.R`. Do not edit by hand.",
        "",
        paste0("PASS: ", tab[["PASS"]], ", EXPLAINED: ", tab[["EXPLAINED"]],
               ", FAIL: ", tab[["FAIL"]], ", INFO: ", tab[["INFO"]]),
        "",
        "| Source | Quantity | Published | Computed | Tolerance | Judgment |",
        "|---|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %s | %s | %s |", res$source, res$quantity,
                format(res$published), ifelse(is.na(res$computed), "-",
                                              sprintf("%.4f", res$computed)),
                ifelse(is.na(res$tolerance), "-", format(res$tolerance)), res$judgment),
        "")
ex <- res[res$judgment %in% c("EXPLAINED", "INFO"), ]
if (nrow(ex) > 0) {
  md <- c(md, "## Reasons", "", sprintf("- %s, %s: %s", ex$source, ex$quantity, ex$reason), "")
}
writeLines(md, file.path(out_dir, "reproduce_summary.md"))
print(tab)
