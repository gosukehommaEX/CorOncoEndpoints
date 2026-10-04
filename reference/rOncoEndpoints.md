# Generate correlated PFS, OS and response

Generates patient-level progression-free survival (PFS), overall
survival (OS), objective response and, optionally, time to response for
one or more treatment groups and many simulated trials, with accrual and
dropout.

## Usage

``` r
rOncoEndpoints(
  nsim = 1,
  n,
  arms,
  a.time = NULL,
  a.rate = NULL,
  d.hazard = 0,
  seed = NULL
)
```

## Arguments

- nsim:

  Number of simulated trials.

- n:

  Integer vector of sample sizes, one per group.

- arms:

  An `OncoArm` object or a list of them, one per group. Group labels are
  taken from the names of the list, otherwise from the `label` of each
  object, otherwise `"Group1"`, `"Group2"`, and so on.

- a.time:

  Accrual breakpoints starting at 0, or `NULL`.

- a.rate:

  Relative accrual intensities on the intervals of `a.time` (default:
  equal).

- d.hazard:

  Dropout hazard, a scalar or one value per group.

- seed:

  Optional integer seed for
  [`dqrng::dqset.seed()`](https://daqana.github.io/dqrng/reference/dqrng-functions.html).

## Value

A data frame with one row per patient, ordered by `sim` and group, with
columns `sim`, `group`, `accrual_time`, `pfs_time`, `os_time`,
`progression` (1 if the PFS event is a progression, 0 if a death),
`response`, `ttr` (time to response, `NA` unless `resp.timing = "ttr"`
and the patient responded), `dropout_time`, `pfs_tte`, `pfs_event`,
`os_tte`, `os_event`, `pfs_calendar_time` and `os_calendar_time`.

## Details

Each group is described by an
[`OncoArm`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md)
object; see its Details for the model. Random numbers come from the
`dqrng` generator and the transformation runs in C++. For reproducible
results set `seed`; the `dqrng` generator is not affected by
[`set.seed()`](https://rdrr.io/r/base/Random.html).

Accrual is piecewise uniform: patients enter independently, with density
proportional to `a.rate` on the intervals defined by `a.time`. Without
`a.time` all patients enter at time 0. Dropout is exponential with
hazard `d.hazard` and censors PFS, OS and the observation of response.

The latent times `pfs_time` and `os_time` always satisfy
`pfs_time <= os_time`, with equality when the PFS event is a death
(`progression = 0`). The observed columns `pfs_tte`, `pfs_event`,
`os_tte` and `os_event` account for dropout but not for an analysis
cutoff; use
[`CutoffData`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md)
for that.

## See also

[`OncoArm`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md),
[`CutoffData`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md)

## Examples

``` r
ctl <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40,
               death.prop = 0.15, os.median = 15, pps.hr.resp = 0.6)
trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.40,
               death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6)
dat <- rOncoEndpoints(nsim = 2, n = c(100, 100),
                      arms = list(Control = ctl, Treatment = trt),
                      a.time = c(0, 24), seed = 1)
head(dat)
#>   sim   group accrual_time  pfs_time   os_time progression response ttr
#> 1   1 Control     1.293824  6.312920  6.312920           0        1  NA
#> 2   1 Control     6.826035  5.339498  5.339498           0        0  NA
#> 3   1 Control    12.854298  4.687368 29.327064           1        1  NA
#> 4   1 Control    14.767060  4.445474  7.748702           1        0  NA
#> 5   1 Control    16.533835  0.341015 47.863850           1        0  NA
#> 6   1 Control     9.079699 26.618913 44.321493           1        1  NA
#>   dropout_time   pfs_tte pfs_event    os_tte os_event pfs_calendar_time
#> 1          Inf  6.312920         1  6.312920        1          7.606744
#> 2          Inf  5.339498         1  5.339498        1         12.165534
#> 3          Inf  4.687368         1 29.327064        1         17.541666
#> 4          Inf  4.445474         1  7.748702        1         19.212534
#> 5          Inf  0.341015         1 47.863850        1         16.874850
#> 6          Inf 26.618913         1 44.321493        1         35.698612
#>   os_calendar_time
#> 1         7.606744
#> 2        12.165534
#> 3        42.181362
#> 4        22.515762
#> 5        64.397685
#> 6        53.401192
```
