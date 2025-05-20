# CorrSurvBinom
<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/CorrSurvBinom)](http://cran.r-project.org/package=CorrSurvBinom)
<!-- badges: end -->
An R package for analyzing of correlated time-to-event and binary outcomes.

## Overview

The `CorrSurvBinom` package provides functionality to generate paired outcomes of time-to-event and binary variables with specified correlation structures using copula models. This package is particularly useful for simulation studies in biostatistics, survival analysis, and clinical trial design.

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/CorrSurvBinom")
```

## Usage

```r
library(CorrSurvBinom)

# Generate 1000 observations with median survival time of 12 months,
# response probability of 0.2, and correlation of 0.5 using Clayton copula
sim_data <- rCorrExpBinom(hazard = log(2) / 12, p = 0.2, rho = 0.5, n = 1000, copula.fun = 'Clayton')

# Access the simulated dataset
head(sim_data$simdata[[1]])

# View summary statistics
sim_data[c("hat.MST", "hat.p", "hat.rho")]
```

## Function Arguments

- `hazard`: Hazard rate (= log(2) / median survival time)
- `p`: Probability of patients being a responder
- `rho`: Correlation between time-to-event and binary outcomes
- `n`: Sample size
- `copula.fun`: Copula function to generate pseudorandom numbers (options: 'Clayton' or 'Frank')

## Return Values

The function returns a tibble with the following columns:

- `hazard`: True hazard rate (= log(2) / median survival time)
- `MST`: True median survival time
- `p`: True probability of responder
- `rho`: True correlation between two outcomes
- `UL.rho`: Upper bound of the correlation
- `n`: Sample size
- `copula.fun`: Specified copula function
- `theta`: Estimated association parameter
- `simdata`: Simulated dataset including correlated time-to-event and binary outcomes
- `hat.MST`: Estimated median survival time for simulated dataset
- `hat.p`: Estimated probability of patients being a responder for simulated dataset
- `hat.rho`: Estimated correlation between two outcomes for simulated dataset

## Notes

Since correlation of paired time-to-event and binary outcomes is not free to range over [-1, 1], the boundary of correlation can be calculated by Frechet-Hoeffding bounds. This function currently only covers positive correlation cases.

## References

1. Trivedi PK and Zimmer DM. Copula modeling: an introduction for practitioners. Found Trends Econ 2007; 1: 1-111.
2. Hofert M, Kojadinovic I, Machler M and et al. Elements of Copula Modeling with R. Use R! 2018; 1: 1-267.
3. Bender R, Augustin T, and Blettner M. Generating survival times to simulate Cox proportional hazards models. Statistics in Medicine 2005; 24: 1713-1723
