# 2-in-1 adaptive design

``` r

library(CorOncoEndpoints)
```

In the 2-in-1 design of Chen et al. (2018), a statistic `X` at an
interim analysis decides whether a phase 2 trial is kept as it is or
expanded into a phase 3 trial. With standardized statistics `Y` at the
end of the small trial and `Z` at the end of the large trial, the study
is positive if `X < c` and `Y > w`, or if `X >= c` and `Z > w`, where
`c` is the expansion cutpoint and `w` the critical value of a one-sided
test. Chen et al. showed that the type I error is controlled for every
`c` when the statistics are bivariate normal and `X` is at least as
strongly correlated with `Y` as with `Z`. When `X` is based on response,
`Y` on PFS and `Z` on OS, these correlations follow from the joint
distribution of the endpoints and the timing of the analyses, and the
generator allows them to be evaluated. This vignette uses a smaller
number of simulated trials than the article on the package.

## Design and groups

The small trial has 120 patients enrolled over 12 months. `X` is the
standardized difference in response rates among the first 90 patients, 3
months after the last of them has entered, and counts only responses
observed by then. The small trial ends at 70 PFS events, and `Y` is the
log-rank statistic for PFS. In an expansion, 380 further patients are
enrolled at 30 per month from month 12, and the large trial ends at 330
deaths, with `Z` the log-rank statistic for OS. Under the global null
hypothesis, both groups have the same distribution:

``` r

ctl <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, death.prop = 0.15,
               os.median = 15, pps.hr.resp = 0.6, resp.timing = "ttr",
               resp.tau = 1.5, ttr.median = 2.5)
arms <- list(Control = ctl, Treatment = ctl)
```

## Simulated trials

Both parts of every trial are generated, and whether the expansion takes
place is decided afterwards from `X`, so that every cutpoint is
evaluated on the same trials. The expansion starts at month 12
(`a.time = 12`) and enrolls 30 patients per month until its 380 patients
are enrolled; it uses another random-number stream than the small part.

``` r

nsim <- 1000
small <- rOncoEndpoints(nsim = nsim, n = c(60, 60), arms = arms,
                        a.time = c(0, 12), seed = 2018)
expansion <- rOncoEndpoints(nsim = nsim, n = c(190, 190), arms = arms,
                            a.time = 12, a.rate = 30, seed = 2018, stream = 1)
large <- rbind(small, expansion)
range(expansion$accrual_time)
#> [1] 12.00000 24.66661
```

The statistics are computed with
[`CutoffData()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md)
and
[`EventTime()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/EventTime.md),
the pooled test of two proportions for `X`, and
[`survival::survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html)
for the log-rank statistics, all oriented so that positive values favor
the treatment group.

``` r

# pooled Z statistic for the difference in observed response proportions
orr_z <- function(d) {
  x <- tapply(d$response_obs, list(d$sim, d$group), sum)
  m <- tapply(d$response_obs, list(d$sim, d$group), length)
  p <- (x[, "Control"] + x[, "Treatment"]) / (m[, "Control"] + m[, "Treatment"])
  se <- sqrt(p * (1 - p) * (1 / m[, "Control"] + 1 / m[, "Treatment"]))
  ifelse(se > 0, (x[, "Treatment"] / m[, "Treatment"] -
                    x[, "Control"] / m[, "Control"]) / se, 0)
}
# log-rank statistic per simulated trial
logrank_z <- function(d, tte, event) {
  vapply(split(d, d$sim), function(x) {
    s <- survival::survdiff(survival::Surv(time, status) ~ group,
                            data = data.frame(time = x[[tte]], status = x[[event]],
                                              group = x$group))
    -(s$obs[2] - s$exp[2]) / sqrt(s$var[2, 2])
  }, numeric(1))
}
# interim analysis: first 90 patients of the small part
small <- small[order(small$sim, small$accrual_time), ]
first90 <- small[ave(small$accrual_time, small$sim, FUN = seq_along) <= 90, ]
t_x <- tapply(first90$accrual_time, first90$sim, max) + 3
X <- orr_z(CutoffData(first90, t_x))
# end of the small trial (70 PFS events) and of the large trial (330 deaths)
Y <- logrank_z(CutoffData(small, EventTime(small, 70, "pfs")), "pfs_tte", "pfs_event")
Z <- logrank_z(CutoffData(large, EventTime(large, 330, "os")), "os_tte", "os_event")
round(cor(cbind(X, Y, Z)), 3)
#>       X     Y     Z
#> X 1.000 0.368 0.148
#> Y 0.368 1.000 0.248
#> Z 0.148 0.248 1.000
```

The correlation of `X` with `Y` is 0.37 and that with `Z` is 0.15. The
interim population overlaps more with the small trial (90 of 120
patients) than with the large trial (90 of 500), and response is more
strongly associated with PFS than with OS.

## Type I error

``` r

w <- qnorm(0.975)
type1 <- function(cc) mean(ifelse(X >= cc, Z > w, Y > w))
cutpoints <- c(0, 0.5, 1, 1.645, 2, 3)
data.frame(cutpoint = cutpoints,
           expanded = sapply(cutpoints, function(cc) mean(X >= cc)),
           type1 = sapply(cutpoints, type1))
#>   cutpoint expanded type1
#> 1    0.000    0.514 0.018
#> 2    0.500    0.286 0.016
#> 3    1.000    0.151 0.022
#> 4    1.645    0.055 0.027
#> 5    2.000    0.027 0.032
#> 6    3.000    0.003 0.034
```

With 1000 simulated trials, the Monte Carlo standard error of a type I
error rate near 0.025 is about 0.005. For a large cutpoint almost no
trial is expanded, and the rate approaches the proportion of trials with
`Y > w`, here 0.034, which differs from 0.025 by 1.8 standard errors;
for a small cutpoint it approaches the proportion with `Z > w`, here
0.028. At the cutpoint 1.645, 6% of the trials were expanded, and the
estimated type I error rate was 0.027. A stronger association of
response with OS, for example through a longer survival of responders
after progression (`pps.hr.resp` further below 1), increases the
correlation of `X` with `Z`, and the generator shows how the operating
characteristics of a planned design depend on it.

## Reference

Chen, C., Anderson, K., Mehrotra, D. V., Rubin, E. H. and Tse, A.
(2018). A 2-in-1 adaptive phase 2/3 design for expedited oncology drug
development. *Contemporary Clinical Trials*, 64, 238-242.
<https://doi.org/10.1016/j.cct.2017.09.006>
