## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## Test environments

* local: Windows 11, R 4.3.2
* GitHub Actions:
  - ubuntu-latest (devel, release, oldrel-1)
  - windows-latest (release)
  - macos-latest (release)

## Downstream dependencies

There are currently no downstream dependencies for this package.

## Comments

This is the initial CRAN submission of CorOncoEndpoints.

The package provides tools for generating correlated oncology endpoints 
(overall survival, progression-free survival, and objective response) 
using copula-based methods. It is particularly useful for simulation 
studies in oncology clinical trials.

Key features:
- Flexible endpoint generation (7 different patterns)
- Support for Clayton and Frank copulas
- Validation tools for simulation results
- Comprehensive documentation with three vignettes

All examples run successfully and all tests pass.
