# Number of events required by the log-rank test

Computes the number of PFS or OS events required for a one-sided
log-rank test to have a given power, using the Schoenfeld formula with
the average hazard ratio of
[`AverageHR`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/AverageHR.md),
and the calendar time at which this number of events is expected.

## Usage

``` r
RequiredEvents(
  arms,
  n,
  a.time = NULL,
  a.rate = NULL,
  d.hazard = 0,
  endpoint = c("os", "pfs"),
  alpha = 0.025,
  power = 0.8,
  tmax = NULL
)
```

## Arguments

- arms:

  A list of two `OncoArm` objects, control first.

- n:

  Sample sizes of the two groups.

- a.time, a.rate:

  Accrual as in
  [`rOncoEndpoints`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md).

- d.hazard:

  Dropout hazard, a scalar or one value per group.

- endpoint:

  `"os"` or `"pfs"`.

- alpha:

  One-sided significance level.

- power:

  Target power.

- tmax:

  Upper limit of calendar time searched (default: end of accrual plus 10
  times the largest median of the endpoint).

## Value

A list with elements `events` (rounded up), `events.exact`, `time`
(calendar time at which `events` are expected), `ahr`, `max.events`
(expected events by `tmax`) and `iterations`.

## Details

With allocation proportions `r1` and `r2`, the required number of events
is \$\$D = (z\_{1 - alpha} + z\_{power})^2 / (r1 r2 log(AHR)^2).\$\$
Because the average hazard ratio depends on the analysis time when the
hazard ratio varies over time (as it does for OS under
`os.model = "idm"`), the number of events and the analysis time are
found by iteration: the analysis time is the calendar time at which `D`
events are expected, and the average hazard ratio is recomputed at that
time until the rounded number of events no longer changes. For
exponential PFS in both groups the hazard ratio is constant and the
result is the usual Schoenfeld number.

## Examples

``` r
ctl <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
               os.median = 15, pps.hr.resp = 0.6)
trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.4,
               death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6)
RequiredEvents(list(Control = ctl, Treatment = trt), n = c(350, 350),
               a.time = c(0, 24), endpoint = "os")
#> $events
#> [1] 366
#> 
#> $events.exact
#> [1] 365.75
#> 
#> $time
#> [1] 30.61522
#> 
#> $ahr
#> [1] 0.7460352
#> 
#> $max.events
#> [1] 699.9859
#> 
#> $iterations
#> [1] 3
#> 
```
