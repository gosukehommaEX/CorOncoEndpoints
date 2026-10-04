# CorOncoEndpoints: Correlated Progression-Free Survival, Overall Survival and Objective Response

The package generates progression-free survival (PFS), overall survival
(OS) and binary objective response for oncology trial simulations so
that the three endpoints of a patient are dependent in a clinically
interpretable way.

Each treatment group is described by an object created with
[`OncoArm`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md).
PFS follows an exponential distribution exactly and the response
probability equals the specified objective response rate exactly. PFS
and response are linked through a Gaussian copula. OS is obtained from
an illness-death structure: a PFS event is a death with a specified
probability and otherwise a progression, after which the patient
survives for a post-progression time. Two models for OS are available.

- `os.model = "idm"`:

  Time-constant pre-progression hazards and an exponential
  post-progression survival whose hazard may differ from the
  pre-progression death hazard and may depend on response. OS is then a
  mixture of exponential distributions with a closed-form or
  one-dimensional integral survival function.

- `os.model = "expexp"`:

  Both PFS and OS are exactly exponential. The pre-progression death
  hazard decreases over time and the post-progression hazard is
  determined by the exponential OS distribution, so it cannot be
  specified separately and does not depend on response.

Main functions:

- [`OncoArm`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md):

  Define one treatment group from design inputs (medians, response rate,
  association) and calibrate the generator.

- [`rOncoEndpoints`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md):

  Generate patient-level data for one or more groups and many simulated
  trials, with accrual and dropout.

- [`CorEndpoints`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CorEndpoints.md):

  Pearson correlations among PFS, OS and response implied by a group.

- [`SurvEndpoint`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/SurvEndpoint.md),
  [`QuantileEndpoint`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/QuantileEndpoint.md):

  Survival, density, hazard and quantiles of PFS and OS, overall or by
  response.

- [`CorBoundPFSResponse`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CorBoundPFSResponse.md):

  Attainable range of the PFS and response correlation.

- [`ExpectedEvents`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/ExpectedEvents.md),
  [`AverageHR`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/AverageHR.md),
  [`RequiredEvents`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/RequiredEvents.md):

  Expected events, average hazard ratio and the number of events
  required by the log-rank test.

- [`CutoffData`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md),
  [`EventTime`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/EventTime.md):

  Data at an analysis cutoff and event-driven cutoff times.

## References

Fleischer, F., Gaschler-Markefski, B. and Bluhmki, E. (2009). A
statistical model for the dependence between progression-free survival
and overall survival. *Statistics in Medicine*, 28, 2669–2686.
[doi:10.1002/sim.3637](https://doi.org/10.1002/sim.3637)

Meller, M., Beyersmann, J. and Rufibach, K. (2019). Joint modeling of
progression-free and overall survival and computation of correlation
measures. *Statistics in Medicine*, 38, 4270–4289.
[doi:10.1002/sim.8295](https://doi.org/10.1002/sim.8295)

## See also

Useful links:

- <https://github.com/gosukehommaEX/CorOncoEndpoints>

- <https://gosukehommaEX.github.io/CorOncoEndpoints/>

- Report bugs at
  <https://github.com/gosukehommaEX/CorOncoEndpoints/issues>

## Author

**Maintainer**: Gosuke Homma <my.name.is.gosuke@gmail.com>

Authors:

- Gosuke Homma <my.name.is.gosuke@gmail.com>
