# Validation of the functions in inst/paper/functions and of the hazard
# formulas used for Figure 1
#
# Run from the package root, after devtools::load_all() or
# library(CorOncoEndpoints):
#   source("inst/paper/validation/validate_example_functions.R")
# Requires the survival and BuyseTest packages. Writes
# inst/paper/validation/output/validate_example_functions.csv and
# inst/paper/validation/output/validate_example_functions.md (generated, do not
# edit by hand).
#
# Judgments: PASS (agreement within the tolerance), EXPLAINED (difference with
# a registered reason), FAIL, INFO. Expected values of the hand-made data set
# were computed independently with loop-based Python code
# (inst/paper/validation/python/expected_values.py).

if (!exists("OncoArm")) library(CorOncoEndpoints)
source(file.path("inst", "paper", "load_functions.R"))
for (p in c("survival", "BuyseTest")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("The ", p, " package is required for this script.", call. = FALSE)
  }
}
# BuyseTest is attached so that its S4 methods for coef() and confint() are used
suppressPackageStartupMessages(library(BuyseTest))
message("Using survival ", utils::packageVersion("survival"), " and BuyseTest ",
        utils::packageVersion("BuyseTest"))
out_dir <- file.path("inst", "paper", "validation", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rows <- list()
add <- function(check, quantity, expected, computed, tol, reason = "") {
  diff <- abs(computed - expected)
  judg <- if (is.na(diff)) "FAIL" else if (diff <= tol) "PASS" else if (nzchar(reason)) "EXPLAINED" else "FAIL"
  rows[[length(rows) + 1L]] <<- data.frame(check = check, quantity = quantity,
                                           expected = expected, computed = computed,
                                           abs_diff = diff, tolerance = tol,
                                           judgment = judg, reason = reason,
                                           stringsAsFactors = FALSE)
}

# ---- 1. hand-made data set: expected values from expected_values.py ----
hm <- data.frame(
  trt = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0),
  os_tte = c(5.2, 8.1, 3.3, 12.0, 7.5, 2.0, 4.4, 8.1, 6.0, 1.5, 9.9, 3.3, 10.4),
  os_event = c(1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1),
  pfs_tte = c(3.1, 8.1, 2.0, 6.5, 7.5, 2.0, 2.2, 5.0, 6.0, 1.5, 4.0, 3.3, 6.6),
  pfs_event = c(1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1),
  resp = c(1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0))
lr_os <- logrank_z(rep(1, 13), hm$trt, hm$os_tte, hm$os_event)
lr_pfs <- logrank_z(rep(1, 13), hm$trt, hm$pfs_tte, hm$pfs_event)
add("hand-made data vs Python", "log-rank OS: O - E", -0.215073815073815, lr_os$o_minus_e, 1e-12)
add("hand-made data vs Python", "log-rank OS: variance", 2.18171528204162, lr_os$var, 1e-12)
add("hand-made data vs Python", "log-rank OS: Z", 0.145609094801767, lr_os$z, 1e-12)
add("hand-made data vs Python", "log-rank PFS: O - E", -1.026221001221, lr_pfs$o_minus_e, 1e-12)
add("hand-made data vs Python", "log-rank PFS: variance", 2.32897492625744, lr_pfs$var, 1e-12)
add("hand-made data vs Python", "log-rank PFS: Z", 0.672447667668736, lr_pfs$z, 1e-12)
add("hand-made data vs Python", "response: pooled Z",
    0.791697994367621, orr_z(rep(1, 13), hm$trt, hm$resp)$z, 1e-12)
ws <- win_statistics(hm$trt == 1, hm$os_tte, hm$os_event, hm$pfs_tte, hm$pfs_event, hm$resp)
py_win <- c(p_win = 0.523809523809524, p_loss = 0.452380952380952,
            p_tie = 0.0238095238095238, wr = 1.15789473684211, wo = 1.15384615384615,
            nb = 0.0714285714285715, se_log_wr = 0.698654710055471,
            se_log_wo = 0.682306145214951, se_nb = 0.339412495706417,
            win_os = 0.380952380952381, loss_os = 0.428571428571429,
            win_pfs = 0.0952380952380952, loss_pfs = 0.0238095238095238,
            win_resp = 0.0476190476190476, loss_resp = 0,
            mwin_os = 0.380952380952381, mloss_os = 0.428571428571429,
            mwin_pfs = 0.5, mloss_pfs = 0.428571428571429,
            mwin_resp = 0.357142857142857, mloss_resp = 0.142857142857143)
for (nm in names(py_win)) {
  add("hand-made data vs Python", paste("win statistics:", nm), py_win[[nm]], ws[[nm]], 1e-12)
}
ind <- win_independence(c(0.30, 0.20, 0.10), c(0.20, 0.15, 0.05))
add("formulas vs Python", "independence: p_win", 0.4325, ind[["p_win"]], 1e-12)
add("formulas vs Python", "independence: p_loss", 0.29125, ind[["p_loss"]], 1e-12)
add("formulas vs Python", "independence: p_tie", 0.27625, ind[["p_tie"]], 1e-12)
add("formulas vs Python", "power of the win ratio test",
    0.652052312293839, win_power(0.45, 0.35, 0.20, 700, measure = "wr"), 1e-12)
add("formulas vs Python", "power of the win odds test",
    0.650404998767263, win_power(0.45, 0.35, 0.20, 700, measure = "wo"), 1e-12)
add("published value", "Yu and Ganju (2022): total sample size for WR 1.5, ties 0.1",
    417, win_sample_size(1.5, 0.1, alpha = 0.025, power = 0.9), 0)

# ---- 2. log-rank and response statistics vs survival and stats ----
arms <- two_group_arms(resp.timing = "ttr", resp.tau = 1.5, ttr.median = 2.5)
sim_d <- rOncoEndpoints(nsim = 20, n = c(60, 60), arms = arms, a.time = c(0, 12),
                        d.hazard = 0.01, seed = 9001)
cut <- CutoffData(sim_d, cutoff = 18)
cut_tie <- cut
cut_tie$os_tte <- round(cut_tie$os_tte, 0)
cut_tie$pfs_tte <- round(cut_tie$pfs_tte, 0)
for (case in c("continuous times", "times rounded to months (ties)")) {
  dd <- if (case == "continuous times") cut else cut_tie
  for (ep in c("os", "pfs")) {
    tte <- dd[[paste0(ep, "_tte")]]
    ev <- dd[[paste0(ep, "_event")]]
    lr <- logrank_z(dd$sim, dd$group == "Treatment", tte, ev)
    ref <- t(vapply(split(seq_len(nrow(dd)), dd$sim), function(i) {
      sd <- survival::survdiff(survival::Surv(tte[i], ev[i]) ~
                                 factor(dd$group[i], levels = c("Control", "Treatment")))
      c(sd$obs[2] - sd$exp[2], sd$var[2, 2])
    }, numeric(2)))
    add("log-rank vs survival::survdiff (20 trials)",
        paste0(toupper(ep), ", ", case, ": max |difference| of O - E"),
        0, max(abs(lr$o_minus_e - ref[, 1])), 1e-8)
    add("log-rank vs survival::survdiff (20 trials)",
        paste0(toupper(ep), ", ", case, ": max |difference| of variance"),
        0, max(abs(lr$var - ref[, 2])), 1e-8)
  }
}
oz <- orr_z(cut$sim, cut$group == "Treatment", cut$response_obs)
ref_z <- vapply(split(seq_len(nrow(cut)), cut$sim), function(i) {
  x <- c(sum(cut$response_obs[i] * (cut$group[i] == "Treatment")),
         sum(cut$response_obs[i] * (cut$group[i] == "Control")))
  n <- c(sum(cut$group[i] == "Treatment"), sum(cut$group[i] == "Control"))
  pt <- suppressWarnings(stats::prop.test(x, n, correct = FALSE))
  sign(x[1] / n[1] - x[2] / n[2]) * sqrt(unname(pt$statistic))
}, numeric(1))
add("response Z vs stats::prop.test (20 trials)", "max |difference| of Z", 0,
    max(abs(oz$z - ref_z)), 1e-10)

# ---- 3. win statistics vs BuyseTest ----
for (s in 1:5) {
  dd <- cut[cut$sim == s, ]
  ws <- win_statistics(dd$group == "Treatment", dd$os_tte, dd$os_event,
                       dd$pfs_tte, dd$pfs_event, dd$response_obs)
  df <- data.frame(arm = factor(dd$group, levels = c("Control", "Treatment")),
                   os_tte = dd$os_tte, os_event = dd$os_event,
                   pfs_tte = dd$pfs_tte, pfs_event = dd$pfs_event,
                   resp = dd$response_obs)
  bt <- BuyseTest(arm ~ tte(os_tte, status = os_event) +
                               tte(pfs_tte, status = pfs_event) + bin(resp),
                             data = df, scoring.rule = "Gehan",
                             method.inference = "u-statistic", trace = 0)
  last <- function(x) unname(x[length(x)])
  chk <- paste0("win statistics vs BuyseTest (trial ", s, ")")
  add(chk, "proportion of pairs won", last(coef(bt, statistic = "favorable")),
      ws[["p_win"]], 1e-10)
  add(chk, "proportion of pairs lost", last(coef(bt, statistic = "unfavorable")),
      ws[["p_loss"]], 1e-10)
  add(chk, "net benefit", last(coef(bt, statistic = "netBenefit")), ws[["nb"]], 1e-10)
  add(chk, "win ratio", last(coef(bt, statistic = "winRatio")), ws[["wr"]], 1e-10)
  ci <- confint(bt, statistic = "netBenefit")
  add(chk, "standard error of the net benefit", ci[nrow(ci), "se"], ws[["se_nb"]], 1e-6)
}

# ---- 4. 2-in-1 statistics recomputed with survival and stats ----
design <- list(n_small = c(60, 60), accrual_small = 12, n_expand = c(190, 190),
               accrual_expand = 380 / 30, m_interim = 90, fu_interim = 3,
               pfs_events_small = 70, os_events_large = 330)
d2 <- two_in_one_data(5, arms, design, seed = 9101)
st <- two_in_one_statistics(d2, design)
add("2-in-1 data", "patients per trial in the small part", 120,
    max(abs(table(d2$sim[d2$part == 1]) - 120)) + 120, 0)
add("2-in-1 data", "small part enrolled before the expansion part (all trials)", 1,
    as.numeric(all(tapply(d2$accrual_time[d2$part == 1], d2$sim[d2$part == 1], max) <=
                     tapply(d2$accrual_time[d2$part == 2], d2$sim[d2$part == 2], min))), 0)
ref_stat <- t(vapply(1:5, function(s) {
  ds <- d2[d2$sim == s, ]
  sm <- ds[ds$part == 1, ]
  sm <- sm[order(sm$accrual_time), ]
  di <- sm[1:90, ]
  ci <- CutoffData(di, max(di$accrual_time) + 3)
  x <- c(sum(ci$response_obs[ci$group == "Treatment"]), sum(ci$response_obs[ci$group == "Control"]))
  n <- c(sum(ci$group == "Treatment"), sum(ci$group == "Control"))
  pt <- suppressWarnings(stats::prop.test(x, n, correct = FALSE))
  zx <- sign(x[1] / n[1] - x[2] / n[2]) * sqrt(unname(pt$statistic))
  lr <- function(dd, ep) {
    sd <- survival::survdiff(survival::Surv(dd[[paste0(ep, "_tte")]], dd[[paste0(ep, "_event")]]) ~
                               factor(dd$group, levels = c("Control", "Treatment")))
    -(sd$obs[2] - sd$exp[2]) / sqrt(sd$var[2, 2])
  }
  ty <- sort(sm$pfs_calendar_time[sm$pfs_event == 1])[70]
  cy <- CutoffData(sm, ty)
  tz <- sort(ds$os_calendar_time[ds$os_event == 1])[330]
  cz <- CutoffData(ds, tz)
  c(zx, lr(cy, "pfs"), lr(cy, "os"), lr(cz, "pfs"), lr(cz, "os"))
}, numeric(5)))
labs <- c("x", "y_pfs", "y_os", "z_pfs", "z_os")
for (k in seq_along(labs)) {
  add("2-in-1 statistics vs survival and stats (5 trials)",
      paste("max |difference| of", labs[k]), 0, max(abs(st[[labs[k]]] - ref_stat[, k])), 1e-8)
}

# ---- 5. hazard formulas of Figure 1 ----
lam_p <- log(2) / 6
lam_o <- log(2) / 15
pi0 <- 0.15
tt <- c(0.1, 0.5, 1, 3, 6, 12, 24, 36)
h02 <- pi0 * lam_p
h12_markov <- lam_o + (lam_o - h02) / (exp((lam_p - lam_o) * tt) - 1)
bal <- h02 * exp(-lam_p * tt) + h12_markov * (exp(-lam_o * tt) - exp(-lam_p * tt)) -
  lam_o * exp(-lam_o * tt)
add("Figure 1 formulas", "Markov model, constant h02: max |death outflow - lam_o S_OS|",
    0, max(abs(bal)), 1e-14)
ee <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, death.prop = pi0,
              os.model = "expexp", os.median = 15)
h02_t <- lam_o * exp(-ee$c_dec * tt)
h12_ee <- (lam_o * exp(-lam_o * tt) - h02_t * exp(-lam_p * tt)) /
  (exp(-lam_o * tt) - exp(-lam_p * tt))
add("Figure 1 formulas", "expexp model: max relative |closed form - package .expexp_h12()|",
    0, max(abs(h12_ee / CorOncoEndpoints:::.expexp_h12(tt, lam_p, lam_o, ee$c_dec) - 1)), 1e-10)
add("Figure 1 formulas", "expexp model: proportion of deaths lam_o / (lam_p + c)", pi0,
    lam_o / (lam_p + ee$c_dec), 1e-12)
th <- log(pi0) / log(lam_o / lam_p)
a_ttp <- (lam_p ^ th - lam_o ^ th) ^ (1 / th)
log_f <- function(s, y) {
  # log of -dS(x, y)/dx at x = s for the Gumbel copula of exponential latent times
  A <- ((a_ttp * s) ^ th + (lam_o * y) ^ th) ^ (1 / th)
  (th - 1) * log(a_ttp * s) + (1 - th) * log(A) - A
}
for (s in c(0.5, 6)) {
  eps <- 1e-6 * s
  h_num <- -(log_f(s, s + 2 * eps) - log_f(s, s)) / (2 * eps)
  add("Figure 1 formulas",
      paste0("Gumbel model: hazard just after progression at s = ", s, " (finite difference)"),
      h02 + (th - 1) * pi0 / s, h_num, 1e-4 * (h02 + (th - 1) * pi0 / s))
}

# ---- output ----
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(out_dir, "validate_example_functions.csv"), row.names = FALSE)
tab <- table(factor(res$judgment, levels = c("PASS", "EXPLAINED", "FAIL", "INFO")))
fmt <- function(x) ifelse(is.na(x), "-", formatC(x, format = "g", digits = 8))
md <- c("# Validation of the functions used for the article",
        "",
        "Generated by `inst/paper/validation/validate_example_functions.R`. Do not edit by hand.",
        "",
        paste0("PASS: ", tab[["PASS"]], ", EXPLAINED: ", tab[["EXPLAINED"]],
               ", FAIL: ", tab[["FAIL"]], ", INFO: ", tab[["INFO"]]),
        "",
        "| Check | Quantity | Expected | Computed | Tolerance | Judgment |",
        "|---|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %s | %s | %s |", gsub("|", "\\|", res$check, fixed = TRUE),
                gsub("|", "\\|", res$quantity, fixed = TRUE),
                fmt(res$expected), fmt(res$computed), fmt(res$tolerance), res$judgment),
        "")
if (any(res$judgment == "EXPLAINED")) {
  md <- c(md, "## Reasons", "",
          paste0("- ", res$check[res$judgment == "EXPLAINED"], ", ",
                 res$quantity[res$judgment == "EXPLAINED"], ": ",
                 res$reason[res$judgment == "EXPLAINED"]), "")
}
writeLines(md, file.path(out_dir, "validate_example_functions.md"))
print(tab)
