# Average hazard ratio between two groups

Computes an average of the time-varying hazard ratio of the second group
versus the first, weighted by the expected number of events, for an
analysis at given calendar times.

## Usage

``` r
AverageHR(
  arms,
  n,
  a.time = NULL,
  a.rate = NULL,
  d.hazard = 0,
  endpoint = c("os", "pfs"),
  time
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

- time:

  Numeric vector of calendar times.

## Value

A data frame with columns `time`, `events` (expected, total) and `ahr`.

## Details

The average hazard ratio at calendar time `T` is \$\$AHR(T) = exp( int
log HR(x) w(x) dx / int w(x) dx ),\$\$ where `HR(x)` is the ratio of the
hazards of the second and first groups at follow-up time `x` and `w(x)`
is the expected number of events at `x` in both groups among patients
followed for at least `x` by `T` (see
[`ExpectedEvents`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/ExpectedEvents.md)).
When the hazard ratio is constant, `AHR(T)` equals it. With this average
hazard ratio, the Schoenfeld formula gives the number of events required
by the log-rank test (see
[`RequiredEvents`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/RequiredEvents.md)).

## Examples

``` r
ctl <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
               os.median = 15, pps.hr.resp = 0.6)
trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.4,
               death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6)
AverageHR(list(Control = ctl, Treatment = trt), n = c(350, 350),
          a.time = c(0, 24), endpoint = "os", time = c(24, 36))
#>   time   events       ahr
#> 1   24 242.2591 0.7404792
#> 2   36 449.6581 0.7505123
```
