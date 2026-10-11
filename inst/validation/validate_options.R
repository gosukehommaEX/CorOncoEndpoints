# Validation of the options added in version 0.3.0 against the theory
#
# Run from the package root, after devtools::load_all() or
# library(CorOncoEndpoints):
#   source("inst/validation/validate_options.R")
# Writes inst/validation/output/validate_options.csv and
# inst/validation/output/validate_options.md (generated, do not edit by hand).
#
# Groups with two endpoints, the correlation between PFS and OS as input,
# accrual by rates with an open final interval, piecewise exponential dropout
# and the censoring of PFS. The scenarios of validate_generator.R (Table S3 of
# the article) are not changed. For each scenario, 40 batches are generated
# (50,000 patients per batch, or 25 trials per batch for the trial scenario).
# The Monte Carlo estimate is the mean over batches and its standard error is
# the standard deviation over batches divided by sqrt(40). PASS when
# |z| = |estimate - theory| / SE < 4; FAIL otherwise. Structural checks are
# PASS or FAIL. Batch b of scenario s uses the seed 10^6 * 9 + 1000 * b + s.

if (!exists("OncoArm")) library(CorOncoEndpoints)
out_dir <- file.path("inst", "validation", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

nb <- 40
nper <- 50000
rows <- list()
add <- function(scen, quantity, theory, est, se) {
  z <- if (is.na(se) || se == 0) NA_real_ else (est - theory) / se
  judg <- if (is.na(z)) {
    if (isTRUE(all.equal(est, theory))) "PASS" else "FAIL"
  } else if (abs(z) < 4) "PASS" else "FAIL"
  rows[[length(rows) + 1L]] <<- data.frame(scenario = scen, quantity = quantity,
                                           theory = theory, estimate = est, se = se,
                                           z = z, judgment = judg,
                                           stringsAsFactors = FALSE)
}
seed_of <- function(s, b) 9e6 + 1000 * b + s
# mean and standard error over batches of the statistics returned by stat(b)
mc <- function(stat) {
  m <- do.call(rbind, lapply(seq_len(nb), stat))
  list(est = colMeans(m), se = apply(m, 2, stats::sd) / sqrt(nb))
}

# ---- 1. PFS and OS without response ----
po <- OncoArm(pfs.median = 6, death.prop = 0.15, os.median = 15, endpoints = c("pfs", "os"))
struct <- TRUE
r <- mc(function(b) {
  d <- rOncoEndpoints(nsim = 1, n = nper, arms = po, seed = seed_of(1, b))
  struct <<- struct && all(is.na(d$response)) && all(d$pfs_time <= d$os_time)
  c(mean(d$progression == 0), cor(d$pfs_time, d$os_time), mean(d$os_time > 12),
    mean(d$os_time))
})
th <- c(0.15, CorEndpoints(po)[["cor.pfs.os"]], SurvEndpoint(po, 12, "os"),
        1 / po$lam_p + (1 - po$death.prop) / po$gam0)
labs <- c("proportion of deaths before progression", "Corr(PFS, OS)", "S_OS(12)", "mean OS")
for (i in seq_along(labs)) add("pfs_os", labs[i], th[i], r$est[i], r$se[i])
add("pfs_os", "structural checks (response NA, PFS <= OS)", 1, as.numeric(struct), NA)

# ---- 2. Correlation between PFS and OS as input ----
pc <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, pfs.os.cor = 0.55,
              os.median = 15, pps.hr.resp = 0.6)
r <- mc(function(b) {
  d <- rOncoEndpoints(nsim = 1, n = nper, arms = pc, seed = seed_of(2, b))
  c(cor(d$pfs_time, d$os_time), mean(d$os_time > 15), mean(d$progression == 0))
})
th <- c(0.55, 0.5, pc$death.prop)
labs <- c("Corr(PFS, OS) (target)", "S_OS(15) (OS median)",
          "proportion of deaths before progression (found)")
for (i in seq_along(labs)) add("pfs_os_cor", labs[i], th[i], r$est[i], r$se[i])

# ---- 3. PFS and response without OS ----
pr <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, resp.timing = "landmark",
              resp.tau = 1.5, endpoints = c("pfs", "response"))
struct <- TRUE
r <- mc(function(b) {
  d <- rOncoEndpoints(nsim = 1, n = nper, arms = pr, seed = seed_of(3, b))
  rr <- d$response == 1
  struct <<- struct && all(is.na(d$os_time)) && all(d$pfs_time[rr] > 1.5)
  c(mean(rr), cor(d$pfs_time, d$response), mean(d$pfs_time[rr] > 6))
})
th <- c(0.30, CorEndpoints(pr)[["cor.pfs.resp"]], SurvEndpoint(pr, 6, "pfs", "responders"))
labs <- c("response rate", "Corr(PFS, R)", "S_PFS(6), responders")
for (i in seq_along(labs)) add("pfs_resp", labs[i], th[i], r$est[i], r$se[i])
add("pfs_resp", "structural checks (OS NA, landmark)", 1, as.numeric(struct), NA)

# ---- 4. OS and response without PFS ----
osr <- OncoArm(os.median = 15, orr = 0.30, resp.cor = 0.35, resp.timing = "ttr",
               resp.tau = 1.5, ttr.median = 3, endpoints = c("os", "response"))
struct <- TRUE
r <- mc(function(b) {
  d <- rOncoEndpoints(nsim = 1, n = nper, arms = osr, seed = seed_of(4, b))
  rr <- d$response == 1
  struct <<- struct && all(is.na(d$pfs_time)) && all(d$os_time[rr] > 1.5) &&
    all(d$ttr[rr] >= 1.5 & d$ttr[rr] < d$os_time[rr])
  c(mean(rr), cor(d$os_time, d$response), mean(d$os_time > 12),
    mean(d$os_time[rr] > 12), mean(d$os_time[!rr] > 12))
})
th <- c(0.30, 0.35, SurvEndpoint(osr, 12, "os"), SurvEndpoint(osr, 12, "os", "responders"),
        SurvEndpoint(osr, 12, "os", "nonresponders"))
labs <- c("response rate", "Corr(OS, R) (target)", "S_OS(12)", "S_OS(12), responders",
          "S_OS(12), non-responders")
for (i in seq_along(labs)) add("os_resp", labs[i], th[i], r$est[i], r$se[i])
add("os_resp", "structural checks (PFS NA, landmark, TTR range)", 1, as.numeric(struct), NA)

# ---- 5. Trials with accrual rates, piecewise dropout and censoring of PFS ----
ctl <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, death.prop = 0.15,
               os.median = 15, pps.hr.resp = 0.6)
trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.40, death.prop = 0.15,
               pps.hazard = ctl$gam0, pps.hr.resp = 0.6)
arms <- list(Control = ctl, Treatment = trt)
# 10 patients per month on [0, 6), then 40 per month until 700 patients (month 22)
design <- list(n = c(350, 350), a.time = c(0, 6), a.rate = c(10, 40),
               d.hazard = c(0.02, 0.005), d.time = c(0, 12, Inf))
pcens <- c(0.01, 0.03)
ntr <- 25
struct <- TRUE
r <- mc(function(b) {
  d <- do.call(rOncoEndpoints, c(list(nsim = ntr, arms = arms, pfs.cens.hazard = pcens,
                                      seed = seed_of(5, b)), design))
  struct <<- struct && all(tapply(d$accrual_time < 6, list(d$sim, d$group), sum) == 30) &&
    max(d$accrual_time) <= 22
  cut <- CutoffData(d, cutoff = 24)
  c(sum(cut$os_event) / ntr, sum(cut$pfs_event) / ntr, mean(d$dropout_time <= 12),
    mean(d$dropout_time <= 24),
    mean(d$pfs_cens_time[d$group == "Treatment"] <= 12))
})
ee_os <- do.call(ExpectedEvents, c(list(arms = arms, endpoint = "os", time = 24), design))
ee_pfs <- do.call(ExpectedEvents, c(list(arms = arms, endpoint = "pfs", time = 24,
                                         pfs.cens.hazard = pcens), design))
th <- c(ee_os$events, ee_pfs$events, 1 - exp(-0.24), 1 - exp(-0.24 - 0.06),
        1 - exp(-0.03 * 12))
labs <- c("OS events by month 24 (ExpectedEvents)", "PFS events by month 24 (ExpectedEvents)",
          "P(dropout <= 12)", "P(dropout <= 24)", "P(censoring of PFS <= 12), treatment")
for (i in seq_along(labs)) add("trial", labs[i], th[i], r$est[i], r$se[i])
add("trial", "structural checks (30 patients per group before month 6, end at 22)", 1,
    as.numeric(struct), NA)

res <- do.call(rbind, rows)
utils::write.csv(res, file.path(out_dir, "validate_options.csv"), row.names = FALSE)
tab <- table(factor(res$judgment, levels = c("PASS", "FAIL")))
md <- c("# Validation of the options of version 0.3.0",
        "",
        "Generated by `inst/validation/validate_options.R`. Do not edit by hand.",
        "",
        paste0(nb, " batches per scenario; PASS when |z| < 4."),
        "",
        paste0("PASS: ", tab[["PASS"]], ", FAIL: ", tab[["FAIL"]],
               ", largest |z|: ", sprintf("%.2f", max(abs(res$z), na.rm = TRUE))),
        "",
        "| Scenario | Quantity | Theory | Estimate | SE | z | Judgment |",
        "|---|---|---|---|---|---|---|",
        sprintf("| %s | %s | %.5f | %.5f | %s | %s | %s |", res$scenario, res$quantity,
                res$theory, res$estimate,
                ifelse(is.na(res$se), "-", sprintf("%.5f", res$se)),
                ifelse(is.na(res$z), "-", sprintf("%.2f", res$z)), res$judgment),
        "")
writeLines(md, file.path(out_dir, "validate_options.md"))
print(tab)
