# Survival, density and hazard of PFS or OS

Evaluates the survival function, density or hazard of PFS or OS implied
by an `OncoArm` object, for all patients or separately for responders
and non-responders.

## Usage

``` r
SurvEndpoint(
  arm,
  t,
  endpoint = c("os", "pfs"),
  response = c("all", "responders", "nonresponders"),
  type = c("survival", "density", "hazard")
)
```

## Arguments

- arm:

  An object of class `OncoArm`.

- t:

  Numeric vector of non-negative times.

- endpoint:

  `"os"` or `"pfs"`.

- response:

  `"all"`, `"responders"` or `"nonresponders"`.

- type:

  `"survival"`, `"density"` or `"hazard"`.

## Value

Numeric vector of the same length as `t`.

## Details

PFS is exponential for all patients; among responders its survival
function is a bivariate normal orthant probability divided by the
response rate. For OS with `os.model = "idm"` and `pps.hr.resp = 1`, the
survival function of all patients has the closed form of Fleischer et
al. (2009, Theorem 5); otherwise the contribution of patients who
progressed is a one-dimensional integral over the progression time. For
`os.model = "expexp"`, OS of all patients is exponential.

## References

Fleischer, F., Gaschler-Markefski, B. and Bluhmki, E. (2009). A
statistical model for the dependence between progression-free survival
and overall survival. *Statistics in Medicine*, 28, 2669–2686.
[doi:10.1002/sim.3637](https://doi.org/10.1002/sim.3637)

## Examples

``` r
arm <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
               os.median = 15, pps.hr.resp = 0.6)
SurvEndpoint(arm, c(6, 12, 24), endpoint = "os")
#> [1] 0.8161914 0.5966966 0.2830903
SurvEndpoint(arm, c(6, 12, 24), endpoint = "os", response = "responders")
#> [1] 0.9363714 0.8022802 0.5014853
SurvEndpoint(arm, c(6, 12, 24), endpoint = "os", type = "hazard")
#> [1] 0.04565910 0.05723889 0.06506027
```
