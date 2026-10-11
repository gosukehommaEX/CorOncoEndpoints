# Expected number of events by calendar time

Computes the expected number of observed PFS or OS events by given
calendar times, for groups described by
[`OncoArm`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md)
objects, with piecewise-uniform accrual, exponential or piecewise
exponential dropout and, for PFS, censoring of PFS.

## Usage

``` r
ExpectedEvents(
  arms,
  n,
  a.time = NULL,
  a.rate = NULL,
  a.prop = NULL,
  d.hazard = 0,
  d.median = NULL,
  d.time = NULL,
  pfs.cens.hazard = 0,
  endpoint = c("os", "pfs"),
  time
)
```

## Arguments

- arms:

  An `OncoArm` object or a list of them, one per group.

- n:

  Sample sizes, one per group.

- a.time, a.rate, a.prop:

  Accrual as in
  [`rOncoEndpoints`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md).

- d.hazard, d.median, d.time:

  Dropout as in
  [`rOncoEndpoints`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md).

- pfs.cens.hazard:

  Hazard of the censoring of PFS, a scalar or one value per group (used
  for `endpoint = "pfs"` only).

- endpoint:

  `"os"` or `"pfs"`.

- time:

  Numeric vector of calendar times.

## Value

A data frame with columns `time`, `events` (total) and one column per
group.

## Details

For group `j` with `n_j` patients, event density `f_j`, cumulative
dropout hazard `D_j` and accrual distribution function `A`, the expected
number of events by calendar time `T` is the integral of
`n_j f_j(x) exp(-D_j(x)) A(T - x)` over follow-up time `x` from 0 to
`T`; for PFS the integrand is also multiplied by `exp(-c_j x)`, where
`c_j` is `pfs.cens.hazard`. `A` is piecewise linear, as specified by
`a.time`, `a.rate` and `a.prop` in
[`rOncoEndpoints`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md);
the generated data have a fixed number of patients in each accrual
interval, which follows `A` up to rounding. The event densities are
tabulated on a grid of 2001 points over `[0, max(time)]` and the
integral is computed by the midpoint rule.

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
