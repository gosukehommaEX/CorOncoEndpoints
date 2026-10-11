# Accrual, dropout and the choice of endpoints

``` r

library(CorOncoEndpoints)
ctl <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, death.prop = 0.15,
               os.median = 15, pps.hr.resp = 0.6)
trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.40,
               death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6)
arms <- list(Control = ctl, Treatment = trt)
```

This vignette shows how the settings of a trial are specified: accrual,
dropout and the censoring of PFS, the observation of response, groups
with two of the three endpoints, and random-number streams for
simulations run in batches. Time is in months throughout.

## Accrual

Accrual is specified as in the FastSurvival package. `a.time` gives the
breakpoints of the accrual intervals, and the number of patients in each
interval is given by rates in patients per month (`a.rate`) or by
proportions (`a.prop`). The rates refer to all groups together. The
number of patients of each group in each interval is fixed, and the
patients are placed uniformly within their interval.

With one more rate than intervals, the last rate applies to an open
final interval that ends when all patients are enrolled. Below, 10
patients per month are enrolled for 6 months and then 30 per month until
the 300 patients are enrolled.

``` r

d1 <- rOncoEndpoints(nsim = 1, n = c(150, 150), arms = arms,
                     a.time = c(0, 6), a.rate = c(10, 30), seed = 1)
table(interval = cut(d1$accrual_time, c(0, 6, Inf), right = FALSE), d1$group)
#>          
#> interval  Control Treatment
#>   [0,6)        30        30
#>   [6,Inf)     120       120
max(d1$accrual_time)
#> [1] 13.99021
```

The first 6 months enroll 60 patients, 30 per group, and accrual ends at
month 14. With as many rates as intervals, the rates must accrue the
total sample size; otherwise an error suggests the alternatives.
Proportions describe the same pattern without rates:

``` r

d2 <- rOncoEndpoints(nsim = 1, n = c(150, 150), arms = arms,
                     a.time = c(0, 6, 14), a.prop = c(0.2, 0.8), seed = 2)
table(interval = cut(d2$accrual_time, c(0, 6, 14), right = FALSE), d2$group)
#>         
#> interval Control Treatment
#>   [0,6)       30        30
#>   [6,14)     120       120
```

Accrual may start after time 0, for example for the patients of the
expansion of a trial. Here 20 patients per month are enrolled from month
12:

``` r

d3 <- rOncoEndpoints(nsim = 1, n = c(50, 50), arms = arms, a.time = 12,
                     a.rate = 20, seed = 3)
range(d3$accrual_time)
#> [1] 12.06227 16.95176
```

[`ExpectedEvents()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/ExpectedEvents.md),
[`AverageHR()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/AverageHR.md)
and
[`RequiredEvents()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/RequiredEvents.md)
take the same arguments:

``` r

ExpectedEvents(arms, n = c(150, 150), a.time = c(0, 6), a.rate = c(10, 30),
               endpoint = "os", time = c(12, 24, 36))
#>   time    events   Control Treatment
#> 1   12  28.49156  16.26906  12.22251
#> 2   24 136.73264  75.49925  61.23340
#> 3   36 215.87670 115.12934 100.74737
```

## Dropout and censoring of PFS

Dropout is exponential with hazard `d.hazard` (or median `d.median`) and
censors PFS and OS. A proportion `q` lost to follow-up within a year
corresponds to the hazard `-log(1 - q) / 12` per month. Piecewise
exponential dropout is given by hazards in `d.hazard` and breakpoints in
`d.time`.

`pfs.cens.hazard` adds a censoring time that applies to PFS only, for
example the start of another anticancer therapy before progression,
after which OS is still followed. Both censoring mechanisms are
independent of the endpoints.

``` r

q <- 0.05
d4 <- rOncoEndpoints(nsim = 1, n = c(2000, 2000), arms = arms, a.time = c(0, 12),
                     d.hazard = -log(1 - q) / 12, pfs.cens.hazard = 0.02, seed = 4)
c(dropout_within_a_year = mean(d4$dropout_time <= 12),
  pfs_censored_by_therapy = mean(d4$pfs_cens_time < pmin(d4$pfs_time, d4$dropout_time)),
  os_events = mean(d4$os_event))
#>   dropout_within_a_year pfs_censored_by_therapy               os_events 
#>                 0.04425                 0.16425                 0.90975
d5 <- rOncoEndpoints(nsim = 1, n = c(2000, 2000), arms = arms, a.time = c(0, 12),
                     d.hazard = c(0.02, 0.005), d.time = c(0, 12, Inf), seed = 5)
c(dropout_before_12 = mean(d5$dropout_time <= 12), expected = 1 - exp(-0.02 * 12))
#> dropout_before_12          expected 
#>         0.2160000         0.2133721
```

## Observation of response

The column `response` of the generated data is the response of the
patient, whether or not it would be observed.
[`CutoffData()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md)
decides which responses are observed at an analysis. With
`resp.timing = "landmark"`, the response is determined at the first
tumor assessment `resp.tau`, and it is observed only if the patient is
still followed then; with `resp.timing = "ttr"` the time to response
plays this role. In `response_obs`, a patient without an observed
response counts as a non-responder, as in an intention-to-treat
analysis, and `response_eval` indicates the patients who were followed
up to the first assessment, for an analysis restricted to them.

``` r

lm_arm <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, death.prop = 0.15,
                  os.median = 15, resp.timing = "landmark", resp.tau = 1.5)
d6 <- rOncoEndpoints(nsim = 1, n = 20000, arms = lm_arm, a.time = c(0, 12),
                     d.hazard = 0.05, seed = 6)
c6 <- CutoffData(d6, cutoff = 24)
c(response = mean(c6$response), observed_itt = mean(c6$response_obs),
  observed_among_assessed = mean(c6$response_obs[c6$response_eval == 1]))
#>                response            observed_itt observed_among_assessed 
#>               0.3050000               0.2837000               0.3058925
```

Because dropout is independent of response, the proportion among the
patients who were assessed estimates the response rate of 0.3, whereas
the intention-to-treat proportion is expected to be 0.3 times the
probability of being followed up to month 1.5, that is, 0.278. The input
`orr` is therefore the response rate without dropout.

## Two of the three endpoints

A group can be described by two endpoints through `endpoints` in
[`OncoArm()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md);
the columns of the third endpoint are then `NA`. With PFS and OS only,
the response inputs are not needed:

``` r

po <- OncoArm(pfs.median = 6, death.prop = 0.15, os.median = 15,
              endpoints = c("pfs", "os"))
CorEndpoints(po)
#> cor.pfs.resp  cor.os.resp   cor.pfs.os pcor.os.resp 
#>           NA           NA    0.6088695           NA
```

With OS and response only, OS is exponential and is linked to response
directly, so that `resp.cor` is the correlation between OS and response:

``` r

osr <- OncoArm(os.median = 15, orr = 0.30, resp.cor = 0.35,
               endpoints = c("os", "response"))
CorEndpoints(osr)
#> cor.pfs.resp  cor.os.resp   cor.pfs.os pcor.os.resp 
#>           NA         0.35           NA           NA
head(rOncoEndpoints(nsim = 1, n = 5, arms = osr, seed = 7)[
  , c("pfs_time", "os_time", "response")])
#>   pfs_time   os_time response
#> 1       NA  6.102199        0
#> 2       NA 18.407025        1
#> 3       NA 32.557469        1
#> 4       NA 78.446937        1
#> 5       NA 10.191542        0
```

Alternatively, the three endpoints can be generated and PFS ignored. The
correlation between OS and response is then a result of the other
inputs, here 0.381 for the control group, and OS is not exponential.

## Correlation between PFS and OS as input

The proportion of PFS events that are deaths can be replaced by a target
correlation between PFS and OS:

``` r

pc <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, pfs.os.cor = 0.55,
              os.median = 15)
c(death.prop = pc$death.prop, CorEndpoints(pc)["cor.pfs.os"])
#> death.prop cor.pfs.os 
#>  0.2405415  0.5500000
```

The attainable correlations are limited by the medians: with medians of
6 and 15 months, values above about 0.68 cannot be reached, and
[`OncoArm()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md)
then reports the attainable range.

## Simulations in batches

A large simulation can be run in batches with the same `seed` and
different values of `stream`. Each stream is an independent sequence of
random numbers, so a batch gives the same data whether the batches are
run one after another or in parallel.

``` r

batch <- function(b) {
  rOncoEndpoints(nsim = 10, n = c(150, 150), arms = arms, a.time = c(0, 12),
                 seed = 2026, stream = b)
}
identical(batch(1), batch(1))
#> [1] TRUE
identical(batch(1)$pfs_time, batch(2)$pfs_time)
#> [1] FALSE
```
