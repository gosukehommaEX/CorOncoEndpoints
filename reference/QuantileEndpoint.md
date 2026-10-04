# Quantiles of PFS or OS

Computes the time at which the survival function of PFS or OS, for all
patients or by response, equals `1 - prob`. The median is `prob = 0.5`.

## Usage

``` r
QuantileEndpoint(
  arm,
  prob = 0.5,
  endpoint = c("os", "pfs"),
  response = c("all", "responders", "nonresponders")
)
```

## Arguments

- arm:

  An object of class `OncoArm`.

- prob:

  Numeric vector of probabilities in (0, 1).

- endpoint:

  `"os"` or `"pfs"`.

- response:

  `"all"`, `"responders"` or `"nonresponders"`.

## Value

Numeric vector of times.

## Examples

``` r
arm <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
               os.median = 15)
QuantileEndpoint(arm, endpoint = "os")
#> [1] 15
QuantileEndpoint(arm, endpoint = "pfs", response = "responders")
#> [1] 11.38721
```
