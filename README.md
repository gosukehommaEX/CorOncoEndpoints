# CorOncoEndpoints

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/CorOncoEndpoints)](http://cran.r-project.org/package=CorOncoEndpoints)
[![R-CMD-check](https://github.com/gosukehommaEX/CorOncoEndpoints/workflows/R-CMD-check/badge.svg)](https://github.com/gosukehommaEX/CorOncoEndpoints/actions)
[![downloads](http://cranlogs.r-pkg.org/badges/grand-total/CorOncoEndpoints)](https://cran.r-project.org/package=CorOncoEndpoints)
[![downloads](http://cranlogs.r-pkg.org/badges/CorOncoEndpoints)](https://cran.r-project.org/package=CorOncoEndpoints)
[![pkgdown](https://img.shields.io/badge/pkgdown-documentation-blue.svg)](https://gosukehommaEX.github.io/CorOncoEndpoints/)
[![R-CMD-check](https://github.com/gosukehommaEX/CorOncoEndpoints/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gosukehommaEX/CorOncoEndpoints/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

An R package for generating correlated oncology endpoints (OS, PFS, Response) in clinical trial simulations.

**Documentation**: Visit our [pkgdown website](https://gosukehommaEX.github.io/CorOncoEndpoints/) for comprehensive documentation.

## Overview

The `CorOncoEndpoints` package provides a comprehensive set of functions to generate correlated oncology endpoints with appropriate constraints and correlation structures using copula models. This package is particularly useful for simulation studies in oncology clinical trials where multiple endpoints need to be modeled with realistic dependencies.

**📚 [Visit our documentation website](https://gosukehommaEX.github.io/CorOncoEndpoints/)** for detailed guides and examples.

### Key Features

- **Flexible Endpoint Generation**: Generate OS only, PFS only, Response only, or any combination of these endpoints
- **Copula-Based Correlation**: Support for Clayton and Frank copulas to model dependence structures
- **Appropriate Constraints**: Automatically ensures PFS ≤ OS when both endpoints are generated
- **Multi-Group Support**: Generate data for multiple treatment groups simultaneously
- **Validation Tools**: Built-in functions to validate simulation results against theoretical values
- **Correlation Bounds**: Calculate feasible correlation ranges based on Fréchet-Hoeffding bounds

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/CorOncoEndpoints")
```

## Documentation and Vignettes

📖 **Comprehensive guides available on our [pkgdown website](https://gosukehommaEX.github.io/CorOncoEndpoints/)**

Explore our detailed vignettes:

- **[Introduction to CorOncoEndpoints](https://gosukehommaEX.github.io/CorOncoEndpoints/articles/introduction.html)** - Get started with basic usage and examples
- **[Advanced Usage and Examples](https://gosukehommaEX.github.io/CorOncoEndpoints/articles/advanced-usage.html)** - Complex scenarios and sensitivity analyses  
- **[Theoretical Background](https://gosukehommaEX.github.io/CorOncoEndpoints/articles/theoretical-background.html)** - Mathematical foundations and derivations

You can also access vignettes locally after installation:

```r
# View vignette in R
vignette("introduction", package = "CorOncoEndpoints")
vignette("advanced-usage", package = "CorOncoEndpoints")
vignette("theoretical-background", package = "CorOncoEndpoints")

# List all available vignettes
browseVignettes("CorOncoEndpoints")
```

## Main Functions

### Data Generation

- **`rOncoEndpoints()`**: Generate correlated endpoints (OS, PFS, Response) for multiple groups
  - Supports 7 different endpoint patterns
  - Handles multiple treatment groups
  - Ensures PFS ≤ OS constraint automatically

### Validation and Analysis

- **`CheckSimResults()`**: Validate simulation results against theoretical values
  - Calculates bias, relative bias, SE, MSE, and RMSE
  - Compares empirical estimates with theoretical expectations
  
### Correlation Tools

- **`CorResponsePFS()`**: Calculate PFS-Response correlation in the OS-PFS-Response framework
  - Useful for understanding implied correlations when all three endpoints are generated
  
- **`CorBoundResponsePFS()`**: Compute correlation bounds for PFS-Response
  - Accounts for hazard rates in the three-endpoint model
  
- **`CorBoundResponseTTE()`**: Compute general correlation bounds for TTE-Response
  - Based on Fréchet-Hoeffding bounds
  - Depends only on response probability

### Copula Parameters

- **`CopulaParamResponseTTE()`**: Calculate copula parameters for a given correlation
  - Converts correlation to copula parameter (theta)
  - Supports Clayton and Frank copulas

## Usage Examples

### Example 1: Generate OS and Response with Correlation

```r
library(CorOncoEndpoints)

set.seed(123)
data1 <- rOncoEndpoints(
  nsim = 100,
  group = c("Treatment", "Control"),
  n = c(150, 150),
  p = c(0.4, 0.3),
  hazard_OS = c(0.05, 0.07),
  rho_tte_resp = c(0.3, 0.2),
  copula = "Clayton"
)

head(data1)
```

### Example 2: Generate All Three Endpoints (OS, PFS, Response)

```r
set.seed(456)
data2 <- rOncoEndpoints(
  nsim = 100,
  group = c("Experimental", "Standard"),
  n = c(200, 200),
  p = c(0.5, 0.35),
  hazard_OS = c(0.04, 0.06),
  hazard_PFS = c(0.08, 0.10),
  rho_tte_resp = c(0.4, 0.25),
  copula = "Frank"
)

# Check correlation between OS and Response
with(subset(data2, Group == "Experimental"), cor(OS, Response))

# Check correlation between PFS and Response (implied)
with(subset(data2, Group == "Experimental"), cor(PFS, Response))

# Verify PFS <= OS constraint
with(data2, all(PFS <= OS))  # Should be TRUE
```

### Example 3: Validate Simulation Results

```r
set.seed(789)
sim_data <- rOncoEndpoints(
  nsim = 1000,
  group = c("Treatment", "Control"),
  n = c(100, 100),
  p = c(0.4, 0.3),
  hazard_OS = c(0.05, 0.07),
  hazard_PFS = c(0.08, 0.10),
  rho_tte_resp = c(0.3, 0.2),
  copula = "Clayton"
)

validation <- CheckSimResults(
  dataset = sim_data,
  p = c(Treatment = 0.4, Control = 0.3),
  hazard_OS = c(Treatment = 0.05, Control = 0.07),
  hazard_PFS = c(Treatment = 0.08, Control = 0.10),
  rho_tte_resp = c(Treatment = 0.3, Control = 0.2),
  copula = "Clayton"
)

print(validation, n = Inf)
```

### Example 4: Calculate Correlation Bounds

```r
# General TTE-Response bounds (depends only on p)
CorBoundResponseTTE(p = 0.4)

# PFS-Response bounds in OS-PFS-Response framework
CorBoundResponsePFS(
  p = 0.4,
  hazard_OS = 0.05,
  hazard_PFS = 0.08
)

# Calculate implied PFS-Response correlation
CorResponsePFS(
  p = 0.4,
  hazard_OS = 0.05,
  hazard_PFS = 0.08,
  rho_OS_Response = 0.3,
  copula = "Clayton"
)
```

### Example 5: Calculate Copula Parameters

```r
# Calculate theta for Clayton copula
theta_clayton <- CopulaParamResponseTTE(
  p = 0.4,
  rho = 0.3,
  copula = "Clayton"
)

# Calculate theta for Frank copula
theta_frank <- CopulaParamResponseTTE(
  p = 0.4,
  rho = 0.3,
  copula = "Frank"
)

cat("Clayton theta:", theta_clayton, "\n")
cat("Frank theta:", theta_frank, "\n")
```

## Supported Endpoint Patterns

The `rOncoEndpoints()` function supports seven different patterns:

1. **OS only**: Generate overall survival data
2. **PFS only**: Generate progression-free survival data
3. **Response only**: Generate binary response data
4. **OS + Response**: Correlated OS and Response
5. **PFS + Response**: Correlated PFS and Response
6. **OS + PFS**: Correlated OS and PFS (with PFS ≤ OS)
7. **OS + PFS + Response**: All three endpoints correlated (with PFS ≤ OS)

## Copula Families

### Clayton Copula
- Exhibits lower tail dependence
- Cannot model negative dependence (rho > 0 only)
- Good for modeling positive associations in survival data

### Frank Copula
- Flexible for both positive and negative dependence
- Symmetric tail behavior
- Suitable for a wider range of correlation structures

## Theoretical Background

This package implements methods based on:

- **Fleischer Model** (2009): Framework for modeling dependence between OS and PFS
  - PFS = min(OS, TTP) where TTP is time to progression
  - Ensures PFS ≤ OS constraint automatically

- **Copula Theory**: Fréchet-Hoeffding bounds and copula-based dependence modeling
  - Clayton copula for lower tail dependence
  - Frank copula for flexible symmetric dependence

- **Correlation Bounds**: Based on Fréchet-Hoeffding copula bounds
  - Determines feasible correlation ranges
  - Accounts for marginal distributions

## References

1. Fleischer, F., Gaschler-Markefski, B., & Bluhmki, E. (2009). A statistical model for the dependence between progression-free survival and overall survival. *Statistics in Medicine*, 28(21), 2669-2686.

2. Trivedi, P. K., & Zimmer, D. M. (2005). Copula modeling: an introduction for practitioners. *Foundations and Trends in Econometrics*, 1(1), 1-111.

3. Hofert, M., Kojadinovic, I., Maechler, M., & Yan, J. (2018). *Elements of copula modeling with R*. Springer.

4. Nelsen, R. B. (2006). *An introduction to copulas* (2nd ed.). Springer.

## Citation

If you use this package in your research, please cite:

```
Homma, G. (2025). CorOncoEndpoints: Generate Correlated Oncology Endpoints for Clinical 
Trial Simulations. R package version 0.1.0. 
https://github.com/gosukehommaEX/CorOncoEndpoints
```

## License

MIT License

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

**Before contributing, please check:**
- [Our documentation](https://gosukehommaEX.github.io/CorOncoEndpoints/)
- [Open issues](https://github.com/gosukehommaEX/CorOncoEndpoints/issues)

## Issues

If you encounter any problems or have suggestions, please file an issue at:
https://github.com/gosukehommaEX/CorOncoEndpoints/issues

## Author

Gosuke Homma
