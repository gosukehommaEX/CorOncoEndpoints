## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local: Windows 11
* GitHub Actions (R-CMD-check workflow): ubuntu-latest (devel, release,
  oldrel-1), windows-latest (release), macos-latest (release)

## Downstream dependencies

There are currently no downstream dependencies for this package.

## Comments

This is a new submission of CorOncoEndpoints (version 0.3.0).

The package generates correlated progression-free survival, overall survival
and objective response for oncology trial simulations, calibrated to design
inputs (medians, the response rate and the association between progression-free
survival and response). It also computes the implied correlations, survival
functions, expected numbers of events, average hazard ratios and the number of
events required by the log-rank test. Random numbers are generated with
'dqrng' and transformed in C++ via 'Rcpp'.
