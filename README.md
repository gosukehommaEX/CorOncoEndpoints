# CorOncoEndpoints

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/CorOncoEndpoints)](http://cran.r-project.org/package=CorOncoEndpoints)
<!-- badges: end -->

An R package for generating correlated oncology endpoints (OS, PFS, Response) in clinical trial simulations.

## Overview

The `CorOncoEndpoints` package provides functionality to generate correlated oncology endpoints with appropriate constraints and correlation structures using copula models. This package is particularly useful for simulation studies in oncology clinical trials.

## Installation

You can install the development version from GitHub:
```r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/CorOncoEndpoints")
```

## Key Functions

- `rOncoEndpoints()`: Generate correlated endpoints (OS, PFS, Response)
- `CheckSimResults()`: Validate simulation results against theoretical values
- `CorResponsePFS()`: Calculate PFS-Response correlation in OS-PFS-Response framework
- `CorBoundResponsePFS()`: Compute correlation bounds for PFS-Response
- `CorBoundResponseTTE()`: Compute correlation bounds for TTE-Response
- `CopulaParamResponseTTE()`: Calculate copula parameters for given correlation

## Usage Example
```r
library(CorOncoEndpoints)

# Generate correlated OS, PFS, and Response data
set.seed(123)
data <- rOncoEndpoints(
  nsim = 100,
  group = c("Treatment", "Control"),
  n = c(150, 150),
  p = c(0.4, 0.3),
  hazard_OS = c(0.05, 0.07),
  hazard_PFS = c(0.08, 0.10),
  rho_tte_resp = c(0.3, 0.2),
  copula = "Clayton"
)

head(data)
```

## References

Fleischer, F., Gaschler-Markefski, B., & Bluhmki, E. (2009). A statistical model for the dependence between progression-free survival and overall survival. Statistics in Medicine, 28(21), 2669-2686.

Trivedi, P. K., & Zimmer, D. M. (2005). Copula modeling: an introduction for practitioners. Foundations and Trends in Econometrics, 1(1), 1-111.
