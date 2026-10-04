# Pearson correlations among PFS, OS and response

Computes the Pearson correlations `Corr(PFS, R)`, `Corr(OS, R)` and
`Corr(PFS, OS)` implied by an `OncoArm` object, together with the
partial correlation of OS and response given PFS.

## Usage

``` r
CorEndpoints(arm)
```

## Arguments

- arm:

  An object of class `OncoArm`.

## Value

A named numeric vector with elements `cor.pfs.resp`, `cor.os.resp`,
`cor.pfs.os` and `pcor.os.resp` (partial correlation of OS and response
given PFS).

## Details

Write `OS = PFS + Delta * D`, where `Delta` indicates progression as the
first event and `D` is the post-progression survival. For
`os.model = "idm"`, `Delta` and `D` are independent of PFS given
response, and with `m_r = (1 - death.prop) / gam_r` and
`delta = m_1 - m_0`, \$\$Cov(OS, R) = Cov(PFS, R) + p (1 - p) delta,\$\$
\$\$Cov(PFS, OS) = Var(PFS) + delta Cov(PFS, R),\$\$ \$\$Var(OS) =
Var(PFS) + Var(Delta D) + 2 delta Cov(PFS, R),\$\$ where `p` is the
response rate. The partial covariance of OS and response given PFS is
`delta p (1 - p) (1 - Corr(PFS, R)^2)`, so the partial correlation has
the sign of `delta`; it is zero when `pps.hr.resp = 1`, in which case
`Corr(OS, R) = Corr(PFS, R) Corr(PFS, OS)`. `Cov(PFS, R)` is a
one-dimensional integral over the latent normal score of PFS.

For `os.model = "expexp"`, the variances are those of the exponential
distributions, and the covariances are computed from the expected
post-progression survival given the progression time, which is
tabulated.

## References

Fleischer, F., Gaschler-Markefski, B. and Bluhmki, E. (2009). A
statistical model for the dependence between progression-free survival
and overall survival. *Statistics in Medicine*, 28, 2669–2686.
[doi:10.1002/sim.3637](https://doi.org/10.1002/sim.3637)

## Examples

``` r
arm <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
               os.median = 15, pps.hr.resp = 0.6)
CorEndpoints(arm)
#> cor.pfs.resp  cor.os.resp   cor.pfs.os pcor.os.resp 
#>    0.4000000    0.3808329    0.6034093    0.1908293 
```
