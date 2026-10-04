# Number of events for PFS and OS

## Two groups

``` r

ctl <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40, death.prop = 0.15,
               os.median = 15, pps.hr.resp = 0.6)
trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.40,
               death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6)
arms <- list(Control = ctl, Treatment = trt)
```

PFS is exponential in both groups, so its hazard ratio is constant at
0.7. The OS hazard ratio changes over time:

``` r

tt <- c(1, 6, 12, 24, 48)
hr_os <- SurvEndpoint(trt, tt, "os", type = "hazard") /
  SurvEndpoint(ctl, tt, "os", type = "hazard")
round(setNames(hr_os, tt), 3)
#>     1     6    12    24    48 
#> 0.712 0.736 0.754 0.782 0.843
```

## Average hazard ratio and required events

The log-rank test does not require proportional hazards, but the
Schoenfeld formula needs a single hazard ratio.
[`AverageHR()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/AverageHR.md)
averages the log hazard ratio with weights equal to the expected events,
and
[`RequiredEvents()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/RequiredEvents.md)
iterates between the number of events and the analysis time.

``` r

AverageHR(arms, n = c(350, 350), a.time = c(0, 24), endpoint = "os",
          time = c(24, 36, 48))
#>   time   events       ahr
#> 1   24 242.2590 0.7404792
#> 2   36 449.6579 0.7505123
#> 3   48 573.5266 0.7587009
r_os <- RequiredEvents(arms, n = c(350, 350), a.time = c(0, 24), endpoint = "os")
r_pfs <- RequiredEvents(arms, n = c(350, 350), a.time = c(0, 24), endpoint = "pfs")
r_os[c("events", "time", "ahr")]
#> $events
#> [1] 366
#> 
#> $time
#> [1] 30.61522
#> 
#> $ahr
#> [1] 0.7460352
r_pfs[c("events", "time", "ahr")]
#> $events
#> [1] 247
#> 
#> $time
#> [1] 16.76118
#> 
#> $ahr
#> [1] 0.7
```

For OS, 366 events are required, expected at month 30.6. Using instead
the hazard ratio implied by the two OS medians, 0.799, as if OS were
exponential would give 626 events. The script
`inst/validation/validate_design.R` checks by simulation that the number
from
[`RequiredEvents()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/RequiredEvents.md)
gives the target power.
