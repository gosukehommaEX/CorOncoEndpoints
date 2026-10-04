# Expected number of events by calendar time

Computes the expected number of observed PFS or OS events by given
calendar times, for groups described by
[`OncoArm`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md)
objects, with piecewise-uniform accrual and exponential dropout.

## Usage

``` r
ExpectedEvents(
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

  An `OncoArm` object or a list of them, one per group.

- n:

  Sample sizes, one per group.

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

A data frame with columns `time`, `events` (total) and one column per
group.

## Details

For group `j` with `n_j` patients, event density `f_j`, dropout hazard
`d_j` and accrual distribution function `A`, the expected number of
events by calendar time `T` is the integral of
`n_j f_j(x) exp(-d_j x) A(T - x)` over follow-up time `x` from 0 to `T`.
The event densities are tabulated on a grid of 2001 points over
`[0, max(time)]` and the integral is computed by the midpoint rule.

## Examples

``` r
ctl <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
               os.median = 15, pps.hr.resp = 0.6)
trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.4,
               death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6)
ExpectedEvents(list(Control = ctl, Treatment = trt), n = c(350, 350),
               a.time = c(0, 24), endpoint = "os", time = c(12, 24, 36))
#>   time   events  Control Treatment
#> 1   12  58.1432  33.0397   25.1035
#> 2   24 242.2591 133.8728  108.3863
#> 3   36 449.6581 241.5483  208.1098
```
