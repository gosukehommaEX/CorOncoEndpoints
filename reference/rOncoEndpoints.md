# Generate correlated PFS, OS and response

Generates patient-level progression-free survival (PFS), overall
survival (OS), objective response and, optionally, time to response for
one or more treatment groups and many simulated trials, with accrual,
dropout and censoring of PFS.

## Usage

``` r
rOncoEndpoints(
  nsim = 1,
  n,
  arms,
  a.time = NULL,
  a.rate = NULL,
  a.prop = NULL,
  d.hazard = 0,
  d.median = NULL,
  d.time = NULL,
  pfs.cens.hazard = 0,
  seed = NULL,
  stream = NULL
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

  Accrual breakpoints (non-negative and increasing), or `NULL` for entry
  of all patients at time 0.

- a.rate:

  Absolute accrual rates (patients per unit time over all groups), with
  `length(a.time) - 1` or `length(a.time)` elements; see Details.

- a.prop:

  Proportions of patients enrolled in the intervals of `a.time`, with
  `length(a.time) - 1` elements. Give at most one of `a.rate` and
  `a.prop`.

- d.hazard:

  Dropout hazard: a scalar or one value per group, or, with `d.time`,
  the hazards of the pieces (a numeric vector for all groups or a list
  with one vector per group).

- d.median:

  Dropout median, an alternative to `d.hazard` with the same structure
  (`Inf` for no dropout).

- d.time:

  Breakpoints of a piecewise exponential dropout distribution, starting
  at 0 and ending with `Inf`: a numeric vector for all groups or a list
  with one vector per group.

- pfs.cens.hazard:

  Hazard of an exponential censoring time that censors PFS only, a
  scalar or one value per group (0 for none).

- seed:

  Optional integer seed for
  [`dqrng::dqset.seed()`](https://daqana.github.io/dqrng/reference/dqrng-functions.html).

- stream:

  Optional non-negative whole number that selects an independent stream
  of the generator for `seed`; requires `seed`.

## Value

A data frame with one row per patient, ordered by `sim` and group, with
columns `sim`, `group`, `accrual_time`, `pfs_time`, `os_time`,
`progression` (1 if the PFS event is a progression, 0 if a death),
`response`, `ttr` (time to response, `NA` unless `resp.timing = "ttr"`
and the patient responded), `dropout_time`, `pfs_tte`, `pfs_event`,
`os_tte`, `os_event`, `pfs_calendar_time`, `os_calendar_time`,
`pfs_cens_time` (censoring time of PFS, `Inf` without `pfs.cens.hazard`)
and `resp_tau` (landmark of the group, used by
[`CutoffData`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md)).

## Details

Each group is described by an
[`OncoArm`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md)
object; see its Details for the model. Random numbers come from the
`dqrng` generator and the transformation runs in C++. For reproducible
results set `seed`; the `dqrng` generator is not affected by
[`set.seed()`](https://rdrr.io/r/base/Random.html).

**Accrual.** Accrual is piecewise uniform over the intervals defined by
`a.time`, as in the FastSurvival package. Each group has a fixed number
of patients in each interval: the group size is divided among the
intervals in proportion to their expected numbers of patients, with the
remainder going to the intervals with the largest fractional parts, and
the patients are placed uniformly within their interval. The expected
numbers are given in one of three ways.

- `a.rate`: absolute rates, in patients per unit time over all groups.
  With `length(a.time) - 1` rates, the rates must accrue the total of
  `n`. With `length(a.time)` rates, the last rate applies to an open
  final interval whose end is computed so that the total of `n` is
  reached; for example, `a.time = 0` with `a.rate = 20` enrolls 20
  patients per unit time until the total is reached.

- `a.prop`: the proportions of patients in the intervals.

- Neither: proportional to the lengths of the intervals, that is,
  uniform accrual from the first to the last element of `a.time`.

`a.time` may start after 0, for example for patients enrolled in the
expansion of a trial. Without `a.time` all patients enter at time 0.

**Dropout and censoring of PFS.** Dropout is exponential with hazard
`d.hazard` (or median `d.median`), or piecewise exponential with
breakpoints `d.time`, and censors PFS and OS. `pfs.cens.hazard` adds an
exponential censoring time that censors PFS only, for example the start
of another anticancer therapy before progression, while OS remains under
follow-up. Both are independent of the endpoints. Whether a response is
observed before dropout, before the censoring of PFS and before an
analysis cutoff is determined by
[`CutoffData`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md).

**Random numbers.** The groups draw their random numbers in turn. The
censoring time of PFS is drawn only for groups with a positive
`pfs.cens.hazard`, so that data generated without it do not change. With
`stream`, the generator is seeded with `seed` and moved to an
independent stream, so that batches of simulated trials can be generated
with the same `seed` and different streams, in sequence or in parallel.

The latent times `pfs_time` and `os_time` always satisfy
`pfs_time <= os_time`, with equality when the PFS event is a death
(`progression = 0`). The observed columns `pfs_tte`, `pfs_event`,
`os_tte` and `os_event` account for dropout and the censoring of PFS but
not for an analysis cutoff, and `response` is the response of the
patient whether or not it would be observed; use
[`CutoffData`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md)
for the data observed at an analysis. For groups with two endpoints
(`endpoints` in
[`OncoArm`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md)),
the columns of the missing endpoint are `NA`; all groups must have the
same endpoints.

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
#>   os_calendar_time pfs_cens_time resp_tau
#> 1         7.606744           Inf        0
#> 2        12.165534           Inf        0
#> 3        42.181362           Inf        0
#> 4        22.515762           Inf        0
#> 5        64.397685           Inf        0
#> 6        53.401192           Inf        0

# 10 patients per month for 6 months, then 30 per month until 200 patients
# are enrolled (month 6 + 140 / 30); dropout of 5% per year (months)
dat2 <- rOncoEndpoints(nsim = 2, n = c(100, 100),
                       arms = list(Control = ctl, Treatment = trt),
                       a.time = c(0, 6), a.rate = c(10, 30),
                       d.hazard = -log(0.95) / 12, seed = 1)
range(dat2$accrual_time)
#> [1]  0.06413425 10.64764795
```
