# CorOncoEndpoints

CorOncoEndpoints generates correlated progression-free survival (PFS),
overall survival (OS) and binary objective response for oncology
clinical trial simulations.

## Model

- PFS is exactly exponential and the response probability equals the
  specified objective response rate. PFS and response are linked through
  a Gaussian copula, calibrated to the Pearson correlation between PFS
  and response or to the median PFS of responders.
- PFS is never longer than OS. A PFS event is a death with a specified
  probability and otherwise a progression, after which the patient
  survives for a post-progression time.
- With `os.model = "idm"` (illness-death model), the death hazard after
  progression can differ from the death hazard before progression and
  can depend on response. OS is then a mixture of exponential
  distributions.
- With `os.model = "expexp"`, both PFS and OS are exactly exponential.
  The pre-progression death hazard decreases over time and the
  post-progression hazard is determined by the exponential OS
  distribution.
- Response can have a landmark (responders have PFS longer than the
  first tumor assessment) and a time to response, so that the response
  observed at an interim analysis can be derived.

## Installation

``` r

# install.packages("devtools")
devtools::install_github("gosukehommaEX/CorOncoEndpoints")
```

## Example

``` r

library(CorOncoEndpoints)

ctl <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40,
               death.prop = 0.15, os.median = 15, pps.hr.resp = 0.6)
trt <- OncoArm(pfs.median = 6 / 0.7, orr = 0.45, resp.cor = 0.40,
               death.prop = 0.15, pps.hazard = ctl$gam0, pps.hr.resp = 0.6)
ctl                  # calibrated parameters, implied medians and correlations
CorEndpoints(trt)    # Corr(PFS, R), Corr(OS, R), Corr(PFS, OS), partial correlation

arms <- list(Control = ctl, Treatment = trt)
dat <- rOncoEndpoints(nsim = 1000, n = c(350, 350), arms = arms,
                      a.time = c(0, 24), seed = 1)

# number of OS events for 80% power of the one-sided log-rank test at 2.5%
RequiredEvents(arms, n = c(350, 350), a.time = c(0, 24), endpoint = "os")

# data at the calendar time of the 366th OS event in each simulated trial
cut <- CutoffData(dat, EventTime(dat, events = 366, endpoint = "os"))
```

## Functions

| Function | Purpose |
|----|----|
| [`OncoArm()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md) | Define a treatment group and calibrate the generator |
| [`rOncoEndpoints()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md) | Generate patient-level PFS, OS, response and time to response |
| [`CorEndpoints()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CorEndpoints.md) | Implied correlations among PFS, OS and response |
| [`SurvEndpoint()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/SurvEndpoint.md), [`QuantileEndpoint()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/QuantileEndpoint.md) | Survival, density, hazard and quantiles, overall or by response |
| [`CorBoundPFSResponse()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CorBoundPFSResponse.md) | Attainable range of the PFS and response correlation |
| [`ExpectedEvents()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/ExpectedEvents.md), [`AverageHR()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/AverageHR.md), [`RequiredEvents()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/RequiredEvents.md) | Expected events, average hazard ratio and required events |
| [`CutoffData()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md), [`EventTime()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/EventTime.md) | Data at an analysis cutoff and event-driven cutoffs |

## Verification

- `tests/testthat`: expected values computed independently in Python
  (`inst/validation/python`), and published values of Fleischer et
  al. (2009) and of the TrialSimulator documentation.
- `inst/reproduce/reproduce_published.R`: reproduction of published
  numbers.
- `inst/validation/validate_generator.R`: Monte Carlo check of the
  generator against the theory.
- `inst/validation/validate_design.R`: simulated power of the log-rank
  test at the number of events from
  [`RequiredEvents()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/RequiredEvents.md).

## References

Fleischer, F., Gaschler-Markefski, B. and Bluhmki, E. (2009). A
statistical model for the dependence between progression-free survival
and overall survival. *Statistics in Medicine*, 28, 2669-2686.
<https://doi.org/10.1002/sim.3637>

Meller, M., Beyersmann, J. and Rufibach, K. (2019). Joint modeling of
progression-free and overall survival and computation of correlation
measures. *Statistics in Medicine*, 38, 4270-4289.
<https://doi.org/10.1002/sim.8295>

## License

MIT
