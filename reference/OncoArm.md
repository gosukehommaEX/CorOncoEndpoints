# Define a treatment group for correlated PFS, OS and response

Creates an `OncoArm` object that describes one treatment group and
calibrates the parameters of the generator to design inputs: the PFS
median (or hazard), the objective response rate, the association between
PFS and response, the proportion of PFS events that are deaths, and the
OS median or the post-progression survival.

## Usage

``` r
OncoArm(
  pfs.median = NULL,
  pfs.hazard = NULL,
  orr,
  resp.cor = NULL,
  resp.pfs.median = NULL,
  death.prop,
  os.model = c("idm", "expexp"),
  os.median = NULL,
  pps.median = NULL,
  pps.hazard = NULL,
  pps.hr.resp = 1,
  resp.timing = c("none", "landmark", "ttr"),
  resp.tau = 0,
  ttr.median = NULL,
  label = NULL
)
```

## Arguments

- pfs.median, pfs.hazard:

  PFS median or hazard; give exactly one.

- orr:

  Objective response rate, in (0, 1).

- resp.cor:

  Target Pearson correlation between PFS and response.

- resp.pfs.median:

  Target median PFS among responders. Give exactly one of `resp.cor` and
  `resp.pfs.median`.

- death.prop:

  Proportion of PFS events that are deaths without prior progression, in
  \[0, 1) for `"idm"` and in (0, pfs.median / os.median\] for
  `"expexp"`.

- os.model:

  `"idm"` (default) or `"expexp"`; see Details.

- os.median:

  Target OS median. Required for `"expexp"`; for `"idm"` give exactly
  one of `os.median`, `pps.median` and `pps.hazard`.

- pps.median, pps.hazard:

  Median or hazard of post-progression survival of non-responders
  (`"idm"` only).

- pps.hr.resp:

  Post-progression hazard ratio of responders versus non-responders
  (`"idm"` only; must be 1 for `"expexp"`).

- resp.timing:

  `"none"` (default), `"landmark"` or `"ttr"`; see Details.

- resp.tau:

  Landmark time: responders have PFS longer than `resp.tau`. Must be 0
  for `"none"` and positive for `"landmark"`; for `"ttr"` it may be 0
  (no landmark).

- ttr.median:

  Median of the untruncated time to response, larger than `resp.tau`
  (`"ttr"` only).

- label:

  Optional group label used by
  [`rOncoEndpoints`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md).

## Value

An object of class `OncoArm`: a list with the inputs and the calibrated
parameters `lam_p` (PFS hazard), `theta` (copula correlation), `c_resp`
(response threshold), `z_tau` (latent score of the landmark), and for
`"idm"` `gam0` and `gam1`, or for `"expexp"` `lam_o`, `c_dec` and a
tabulated cumulative post-progression hazard.

## Details

**PFS and response.** Let `(Z1, Z2)` be standard bivariate normal with
correlation `theta`. PFS is `-log(Phi(-Z1)) / lam_p`, so it is exactly
exponential with hazard `lam_p`. Response is `R = 1` when
`PFS > resp.tau` and `Z2 > c_resp`, where `c_resp` is chosen so that
`P(R = 1) = orr` exactly. The copula correlation `theta` is calibrated
either to the Pearson correlation `resp.cor` between PFS and response or
to the median PFS among responders `resp.pfs.median`. Exactly one of
these two must be given.

**Timing of response.** With `resp.timing = "none"` response is a
patient-level binary variable without a time. With
`resp.timing = "landmark"` responders must have PFS longer than
`resp.tau` (for example, the time of the first tumor assessment). With
`resp.timing = "ttr"` the same restriction applies and, in addition,
each responder gets a time to response that lies between `resp.tau` and
the responder's PFS: `resp.tau` plus an exponential time with median
`ttr.median - resp.tau`, truncated at `PFS - resp.tau`. The time to
response allows the response observed by an interim analysis to be
derived (see
[`CutoffData`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CutoffData.md)).
In all three modes, PFS remains exactly exponential and the response
probability remains exactly `orr`.

**OS, illness-death model (`os.model = "idm"`).** A PFS event is a death
with probability `death.prop`, independently of PFS and response; the
pre-progression hazards of progression and death are therefore constant,
`(1 - death.prop) * lam_p` and `death.prop * lam_p`. After progression,
the patient survives for an exponential time with hazard `gam0`
(non-responders) or `gam1 = pps.hr.resp * gam0` (responders), measured
from progression. `gam0` is calibrated to `os.median` or given through
`pps.median` or `pps.hazard`. OS is not exponential unless
`pps.hr.resp = 1` and `gam0 = death.prop * lam_p`.

**OS, exponential-exponential model (`os.model = "expexp"`).** OS is
exactly exponential with median `os.median`. This requires the
pre-progression death hazard to start at the OS hazard `lam_o`; it is
set to `lam_o * exp(-c_dec * t)` with `c_dec` chosen so that the
proportion of PFS events that are deaths equals `death.prop`, which must
not exceed `pfs.median / os.median`. The post-progression hazard is then
determined by the exponential OS distribution (it depends on time since
randomization) and cannot depend on response, so `pps.hr.resp` must
be 1. With `death.prop = pfs.median / os.median` the model reduces to
the maximal independence model of Fleischer et al. (2009).

## References

Fleischer, F., Gaschler-Markefski, B. and Bluhmki, E. (2009). A
statistical model for the dependence between progression-free survival
and overall survival. *Statistics in Medicine*, 28, 2669–2686.
[doi:10.1002/sim.3637](https://doi.org/10.1002/sim.3637)

## See also

[`rOncoEndpoints`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md),
[`CorEndpoints`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/CorEndpoints.md)

## Examples

``` r
ctl <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40,
               death.prop = 0.15, os.median = 15, pps.hr.resp = 0.6,
               label = "Control")
ctl
#> OncoArm (Control)
#>   OS model           : idm
#>   Response timing    : none
#>   PFS hazard         : 0.1155
#>   Response rate      : 0.3
#>   Copula correlation : 0.5376
#>   Death proportion   : 0.15
#>   Post-progression hazard (non-responders, responders): 0.09702, 0.05821
#> Implied medians (all, responders, non-responders)
#>   PFS: 6, 11.39, 4.393
#>   OS : 15, 24.07, 12.23
#> Implied correlations
#>   Corr(PFS, R) = 0.4, Corr(OS, R) = 0.3808, Corr(PFS, OS) = 0.6034
#>   Partial correlation of OS and R given PFS = 0.1908

# Both PFS and OS exactly exponential
ee <- OncoArm(pfs.median = 6, orr = 0.30, resp.cor = 0.40,
              death.prop = 0.15, os.model = "expexp", os.median = 15)
QuantileEndpoint(ee, endpoint = "os")
#> [1] 15
```
