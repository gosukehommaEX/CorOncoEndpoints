# Calendar time of a given number of events

Returns, for each simulated trial, the calendar time at which the given
number of observed PFS or OS events is reached. Events censored by
dropout are not counted. The result can be passed to
[`CutoffData`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md)
for event-driven analyses.

## Usage

``` r
EventTime(data, events, endpoint = c("os", "pfs"))
```

## Arguments

- data:

  A data frame returned by
  [`rOncoEndpoints`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md).

- events:

  Number of events.

- endpoint:

  `"os"` or `"pfs"`.

## Value

A numeric vector with one value per simulated trial (`NA` when fewer
events occur), named by `sim`.

## Examples

``` r
arm <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
               os.median = 15)
dat <- rOncoEndpoints(nsim = 3, n = 200, arms = arm, a.time = c(0, 12),
                      seed = 1)
EventTime(dat, events = 100, endpoint = "os")
#>        1        2        3 
#> 22.26770 20.44160 21.58353 
```
