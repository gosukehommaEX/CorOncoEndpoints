# CorOncoEndpoints 0.2.0

## Breaking changes

* The package has been rewritten. A treatment group is now defined with
  `OncoArm()`, and `rOncoEndpoints()` takes a list of `OncoArm` objects
  (`arms`), sample sizes, accrual, dropout and a `seed`.
* Response is now linked to PFS (through a Gaussian copula) instead of to the
  death time. In version 0.1.0, response was linked to OS by a Clayton or Frank
  copula and time to progression was independent of response, so that
  `Corr(PFS, R)` was always smaller than `Corr(OS, R)`.
* Removed `CheckSimResults()`, `CopulaParamResponseTTE()`,
  `CorBoundResponsePFS()`, `CorBoundResponseTTE()` and `CorResponsePFS()`.

## New features

* `OncoArm()` calibrates the generator to design inputs: PFS median or
  hazard, objective response rate, the correlation between PFS and response
  (`resp.cor`) or the median PFS of responders (`resp.pfs.median`), the
  proportion of PFS events that are deaths (`death.prop`), and the OS median
  or the post-progression survival.
* Two OS models: `os.model = "idm"` (illness-death model in which the death
  hazard after progression can differ from the death hazard before progression
  and can depend on response through `pps.hr.resp`) and
  `os.model = "expexp"` (PFS and OS both exactly exponential).
* Timing of response: `resp.timing = "landmark"` (responders have PFS longer
  than `resp.tau`) and `resp.timing = "ttr"` (time to response between
  `resp.tau` and PFS). PFS stays exactly exponential and the response
  probability stays exactly `orr`.
* `CorEndpoints()`, `SurvEndpoint()`, `QuantileEndpoint()` and
  `CorBoundPFSResponse()` give the implied correlations, survival functions,
  quantiles and the attainable range of the PFS and response correlation.
* `ExpectedEvents()`, `AverageHR()` and `RequiredEvents()` compute expected
  events, the average hazard ratio and the number of events required by the
  log-rank test when OS hazards are not proportional.
* `CutoffData()` and `EventTime()` give data at an analysis cutoff, including
  the response observed by the cutoff, and event-driven cutoff times.
* Random numbers are generated with 'dqrng' and transformed in C++ via 'Rcpp'.

## Verification

* Tests compare with values computed independently in Python
  (`inst/validation/python`) and with published values of Fleischer et al.
  (2009) and of the TrialSimulator documentation.
* `inst/reproduce/reproduce_published.R` reproduces published numbers;
  `inst/validation/validate_generator.R` and
  `inst/validation/validate_design.R` check the generator and
  `RequiredEvents()` by simulation.

## Dependencies

* Imports 'Rcpp' and 'dqrng'; 'tibble' is no longer used. R (>= 4.1.0).

# CorOncoEndpoints 0.1.0

## Initial Release

This is the first release of the CorOncoEndpoints package, providing comprehensive tools for generating correlated oncology endpoints in clinical trial simulations.

### Main Features

* **Data Generation**
  - `rOncoEndpoints()`: Generate correlated OS, PFS, and Response endpoints
  - Support for 7 different endpoint patterns (OS only, PFS only, Response only, OS+Response, PFS+Response, OS+PFS, OS+PFS+Response)
  - Multi-group support for comparing multiple treatment arms
  - Automatic enforcement of PFS ≤ OS constraint

* **Copula Support**
  - Clayton copula (lower tail dependence, positive correlations only)
  - Frank copula (flexible symmetric dependence, positive and negative correlations)

* **Validation Tools**
  - `CheckSimResults()`: Validate simulation results against theoretical values
  - Calculate bias, relative bias, SE, MSE, and RMSE
  - Compare empirical estimates with theoretical expectations

* **Correlation Analysis**
  - `CorResponsePFS()`: Calculate PFS-Response correlation in OS-PFS-Response framework
  - `CorBoundResponsePFS()`: Compute correlation bounds for PFS-Response
  - `CorBoundResponseTTE()`: Compute general correlation bounds for TTE-Response
  - `CopulaParamResponseTTE()`: Calculate copula parameters for given correlations

### Theoretical Framework

* Implementation of Fleischer model (2009) for OS-PFS dependency
* Fréchet-Hoeffding bounds for correlation constraints
* Copula-based dependence modeling

### Documentation

* Comprehensive function documentation with roxygen2
* Multiple examples for each function
* Detailed usage guides in README

### Dependencies

* R (>= 3.5.0)
* stats (base R)
* tibble

## Future Plans

* Add support for additional copula families (Gumbel, Gaussian)
* Implement time-to-event endpoints with censoring
* Add visualization functions for correlation structures
* Extend to handle more complex trial designs
