# Attainable range of the correlation between PFS and response

Returns the smallest and largest Pearson correlation between an
exponential PFS and a binary response with probability `orr`, optionally
under the restriction that responders have PFS longer than a landmark
time `resp.tau`.

## Usage

``` r
CorBoundPFSResponse(orr, pfs.median = NULL, pfs.hazard = NULL, resp.tau = 0)
```

## Arguments

- orr:

  Objective response rate, in (0, 1).

- pfs.median, pfs.hazard:

  PFS median or hazard; needed only when `resp.tau > 0`.

- resp.tau:

  Landmark time (0 for none).

## Value

A numeric vector `c(lower = , upper = )`.

## Details

Without a landmark the bounds are the Frechet bounds
`sqrt((1 - p) / p) log(1 - p)` and `-sqrt(p / (1 - p)) log(p)`, which do
not depend on the PFS hazard. With a landmark the upper bound is
unchanged, and the lower bound is attained when the responders are the
patients with the shortest PFS above `resp.tau`. The Gaussian copula
used by
[`OncoArm`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md)
approaches both bounds as the copula correlation tends to -1 or 1;
[`OncoArm()`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/OncoArm.md)
restricts the copula correlation to (-0.9999, 0.9999), so correlations
very close to a bound may not be attainable.

## Examples

``` r
CorBoundPFSResponse(orr = 0.3)
#>      lower      upper 
#> -0.5448300  0.7881852 
CorBoundPFSResponse(orr = 0.3, pfs.median = 6, resp.tau = 1.5)
#>      lower      upper 
#> -0.4073680  0.7881852 
```
