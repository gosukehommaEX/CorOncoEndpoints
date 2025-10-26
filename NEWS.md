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
