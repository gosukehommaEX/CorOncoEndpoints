# Getting started with CorOncoEndpoints

``` r

library(CorOncoEndpoints)
```

## Describing a treatment group

A treatment group is described by
[`OncoArm()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md).
The inputs are the quantities usually discussed at the design stage: the
PFS median, the objective response rate, the association between PFS and
response, the proportion of PFS events that are deaths, and the OS
median.

``` r

ctl <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40,
               death.prop = 0.15, os.median = 15, pps.hr.resp = 0.6,
               label = "Control")
ctl
#> OncoArm (Control)
#>   OS model           : idm
#>   Response timing    : none
#>   PFS hazard         : 0.1155
#>   Response rate      : 0.3
#>   Copula correlation : 0.5376
#>   Death proportion   : 0.15
#>   Post-progression hazard (non-responders, responders): 0.09702, 0.05821
#> Implied medians (all, responders, non-responders)
#>   PFS: 6, 11.39, 4.393
#>   OS : 15, 24.07, 12.23
#> Implied correlations
#>   Corr(PFS, R) = 0.4, Corr(OS, R) = 0.3808, Corr(PFS, OS) = 0.6034
#>   Partial correlation of OS and R given PFS = 0.1908
```

The post-progression hazard of non-responders was calibrated to the OS
median: it is 0.097 per month, so the median post-progression survival
of non-responders is 7.1 months, and responders have 0.6 times this
hazard. The association with response was given as a Pearson
correlation; it can instead be given as the median PFS among responders
through `resp.pfs.median`.

The implied correlations among the three endpoints are

``` r

CorEndpoints(ctl)
#> cor.pfs.resp  cor.os.resp   cor.pfs.os pcor.os.resp 
#>    0.4000000    0.3808329    0.6034093    0.1908293
```

## A second group

The treatment group below has a PFS hazard ratio of 0.7 and a higher
response rate. Its post-progression hazard is set equal to that of the
control group, so the treatment acts on OS only through PFS and
response.

``` r

trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.40,
               death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6,
               label = "Treatment")
QuantileEndpoint(trt, endpoint = "os")
#> [1] 18.76724
```

The implied OS median of the treatment group is 18.8 months.

## Generating data

[`rOncoEndpoints()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md)
generates patient-level data for many simulated trials. The random
numbers come from the `dqrng` generator; set `seed` for reproducible
results.

``` r

dat <- rOncoEndpoints(nsim = 2, n = c(200, 200),
                      arms = list(Control = ctl, Treatment = trt),
                      a.time = c(0, 24), d.hazard = 0.005, seed = 1)
head(dat)
#>   sim   group accrual_time  pfs_time   os_time progression response ttr
#> 1   1 Control   7.25918546  6.312920 29.503975           1        0  NA
#> 2   1 Control  14.80959863  5.339498  8.219314           1        0  NA
#> 3   1 Control  19.16633898  4.687368  4.687368           0        0  NA
#> 4   1 Control  17.07641780  4.445474  9.800482           1        0  NA
#> 5   1 Control   4.91513070  0.341015 14.370445           1        0  NA
#> 6   1 Control   0.01352835 26.618913 32.477194           1        0  NA
#>   dropout_time   pfs_tte pfs_event    os_tte os_event pfs_calendar_time
#> 1     32.19099  6.312920         1 29.503975        1         13.572105
#> 2    112.77597  5.339498         1  8.219314        1         20.149097
#> 3    140.66633  4.687368         1  4.687368        1         23.853707
#> 4    432.44744  4.445474         1  9.800482        1         21.521892
#> 5     32.67168  0.341015         1 14.370445        1          5.256146
#> 6    315.94705 26.618913         1 32.477194        1         26.632441
#>   os_calendar_time
#> 1         36.76316
#> 2         23.02891
#> 3         23.85371
#> 4         26.87690
#> 5         19.28558
#> 6         32.49072
```

Every patient has `pfs_time <= os_time`, with equality when the PFS
event is a death (`progression = 0`). In this example 0 patients violate
the ordering.

## Analysis at a cutoff

[`EventTime()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/EventTime.md)
gives the calendar time at which a number of events is reached in each
simulated trial, and
[`CutoffData()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md)
returns the data observed at a cutoff.

``` r

cut <- CutoffData(dat, EventTime(dat, events = 150, endpoint = "os"))
table(cut$sim, cut$os_event)
#>    
#>       0   1
#>   1 250 150
#>   2 250 150
```

## Response observed at an interim analysis

With `resp.timing = "ttr"` each responder has a time to response, so the
response observed by an early cutoff can be derived. Responders also
have PFS longer than `resp.tau`, here the time of the first tumor
assessment.

``` r

arm_ttr <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40,
                   death.prop = 0.15, os.median = 15,
                   resp.timing = "ttr", resp.tau = 1.5, ttr.median = 2.5)
d_ttr <- rOncoEndpoints(nsim = 1, n = 300, arms = arm_ttr, a.time = c(0, 12),
                        seed = 2)
early <- CutoffData(d_ttr, cutoff = 9)
c(enrolled = nrow(early), responders_eventually = sum(early$response),
  responses_observed = sum(early$response_obs))
#>              enrolled responders_eventually    responses_observed 
#>                   209                    66                    42
```
