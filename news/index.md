# Changelog

## CorOncoEndpoints 0.3.0

### Breaking changes

- `a.rate` is now an absolute accrual rate (patients per unit time over
  all groups), as in the FastSurvival package. With `length(a.time) - 1`
  rates they must accrue the total sample size; with `length(a.time)`
  rates the last one applies to an open final interval. The relative
  intensities of version 0.2.0 are given by the new `a.prop`.
- Each accrual interval now receives a fixed number of patients of each
  group (largest-remainder rounding), placed uniformly within the
  interval, instead of independently drawn entry times. With one accrual
  interval, or none, the data are the same as in version 0.2.0.
- [`CutoffData()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md)
  now counts a response with `resp.timing = "landmark"` as observed only
  if the landmark `resp.tau` was reached before the cutoff and before
  dropout. Responses with `resp.timing = "none"` and `"ttr"` are treated
  as before, unless the new censoring of PFS applies.
- In
  [`ExpectedEvents()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/ExpectedEvents.md),
  [`AverageHR()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/AverageHR.md)
  and
  [`RequiredEvents()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/RequiredEvents.md),
  the new arguments follow `a.rate`, so the arguments after `a.rate`
  must be named.
- `orr` and `death.prop` of
  [`OncoArm()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md)
  have the default `NULL`, because they are not needed for every
  combination of endpoints.

### New features

- `OncoArm(endpoints = )` describes a group by two of the three
  endpoints: PFS and OS, PFS and response, or OS and response. With OS
  and response, OS is exactly exponential and is linked to response by
  the Gaussian copula, so that `resp.cor` is the correlation between OS
  and response and `resp.os.median` the median OS of responders. The
  columns of the missing endpoint are `NA`.
- `OncoArm(pfs.os.cor = )` gives the correlation between PFS and OS
  instead of `death.prop`, which is then found by a one-dimensional
  search.
- Accrual may start after time 0 (`a.time[1] > 0`).
- Dropout can be given by its median (`d.median`) and can be piecewise
  exponential (`d.time`, with group-specific hazards through lists).
- `pfs.cens.hazard` censors PFS (and the assessment of response) but not
  OS, for example at the start of another anticancer therapy.
- [`CutoffData()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md)
  returns `response_eval`, whether the patient was followed up to the
  first tumor assessment.
- `rOncoEndpoints(stream = )` selects an independent random-number
  stream for simulations in batches.
- [`rOncoEndpoints()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md)
  returns the columns `pfs_cens_time` and `resp_tau`.
- The design functions take the same accrual and dropout arguments as
  [`rOncoEndpoints()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md)
  and `pfs.cens.hazard` for PFS.

### Documentation

- New vignettes: “Accrual, dropout and the choice of endpoints”, “2-in-1
  adaptive design” and “Win ratio and win odds of prioritized
  endpoints”.
- The help of
  [`rOncoEndpoints()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md)
  no longer states that dropout censors the observation of response for
  every timing of response;
  [`CutoffData()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md)
  documents when a response counts as observed.

### Verification

- `tests/testthat/fixtures/reference_v0.2.0.rds` stores outputs of
  version 0.2.0, and the tests check that the same calls give the same
  results.
- New tests with expected values computed independently in Python
  (`inst/validation/python`), and `inst/validation/validate_options.R`,
  a Monte Carlo check of the new options.

## CorOncoEndpoints 0.2.0

### Breaking changes

- The package has been rewritten. A treatment group is now defined with
  [`OncoArm()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md),
  and
  [`rOncoEndpoints()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md)
  takes a list of `OncoArm` objects (`arms`), sample sizes, accrual,
  dropout and a `seed`.
- Response is now linked to PFS (through a Gaussian copula) instead of
  to the death time. In version 0.1.0, response was linked to OS by a
  Clayton or Frank copula and time to progression was independent of
  response, so that `Corr(PFS, R)` was always smaller than
  `Corr(OS, R)`.
- Removed `CheckSimResults()`, `CopulaParamResponseTTE()`,
  `CorBoundResponsePFS()`, `CorBoundResponseTTE()` and
  `CorResponsePFS()`.

### New features

- [`OncoArm()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md)
  calibrates the generator to design inputs: PFS median or hazard,
  objective response rate, the correlation between PFS and response
  (`resp.cor`) or the median PFS of responders (`resp.pfs.median`), the
  proportion of PFS events that are deaths (`death.prop`), and the OS
  median or the post-progression survival.
- Two OS models: `os.model = "idm"` (illness-death model in which the
  death hazard after progression can differ from the death hazard before
  progression and can depend on response through `pps.hr.resp`) and
  `os.model = "expexp"` (PFS and OS both exactly exponential).
- Timing of response: `resp.timing = "landmark"` (responders have PFS
  longer than `resp.tau`) and `resp.timing = "ttr"` (time to response
  between `resp.tau` and PFS). PFS stays exactly exponential and the
  response probability stays exactly `orr`.
- [`CorEndpoints()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CorEndpoints.md),
  [`SurvEndpoint()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/SurvEndpoint.md),
  [`QuantileEndpoint()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/QuantileEndpoint.md)
  and
  [`CorBoundPFSResponse()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CorBoundPFSResponse.md)
  give the implied correlations, survival functions, quantiles and the
  attainable range of the PFS and response correlation.
- [`ExpectedEvents()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/ExpectedEvents.md),
  [`AverageHR()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/AverageHR.md)
  and
  [`RequiredEvents()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/RequiredEvents.md)
  compute expected events, the average hazard ratio and the number of
  events required by the log-rank test when OS hazards are not
  proportional.
- [`CutoffData()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md)
  and
  [`EventTime()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/EventTime.md)
  give data at an analysis cutoff, including the response observed by
  the cutoff, and event-driven cutoff times.
- Random numbers are generated with ‘dqrng’ and transformed in C++ via
  ‘Rcpp’.

### Verification

- Tests compare with values computed independently in Python
  (`inst/validation/python`) and with published values of Fleischer et
  al.
  2009. and of the TrialSimulator documentation.
- `inst/reproduce/reproduce_published.R` reproduces published numbers;
  `inst/validation/validate_generator.R` and
  `inst/validation/validate_design.R` check the generator and
  [`RequiredEvents()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/RequiredEvents.md)
  by simulation.

### Dependencies

- Imports ‘Rcpp’ and ‘dqrng’; ‘tibble’ is no longer used. R (\>= 4.1.0).

## CorOncoEndpoints 0.1.0

### Initial Release

This is the first release of the CorOncoEndpoints package, providing
comprehensive tools for generating correlated oncology endpoints in
clinical trial simulations.

#### Main Features

- **Data Generation**
  - [`rOncoEndpoints()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md):
    Generate correlated OS, PFS, and Response endpoints
  - Support for 7 different endpoint patterns (OS only, PFS only,
    Response only, OS+Response, PFS+Response, OS+PFS, OS+PFS+Response)
  - Multi-group support for comparing multiple treatment arms
  - Automatic enforcement of PFS ≤ OS constraint
- **Copula Support**
  - Clayton copula (lower tail dependence, positive correlations only)
  - Frank copula (flexible symmetric dependence, positive and negative
    correlations)
- **Validation Tools**
  - `CheckSimResults()`: Validate simulation results against theoretical
    values
  - Calculate bias, relative bias, SE, MSE, and RMSE
  - Compare empirical estimates with theoretical expectations
- **Correlation Analysis**
  - `CorResponsePFS()`: Calculate PFS-Response correlation in
    OS-PFS-Response framework
  - `CorBoundResponsePFS()`: Compute correlation bounds for PFS-Response
  - `CorBoundResponseTTE()`: Compute general correlation bounds for
    TTE-Response
  - `CopulaParamResponseTTE()`: Calculate copula parameters for given
    correlations

#### Theoretical Framework

- Implementation of Fleischer model (2009) for OS-PFS dependency
- Fréchet-Hoeffding bounds for correlation constraints
- Copula-based dependence modeling

#### Documentation

- Comprehensive function documentation with roxygen2
- Multiple examples for each function
- Detailed usage guides in README

#### Dependencies

- R (\>= 3.5.0)
- stats (base R)
- tibble

### Future Plans

- Add support for additional copula families (Gumbel, Gaussian)
- Implement time-to-event endpoints with censoring
- Add visualization functions for correlation structures
- Extend to handle more complex trial designs
