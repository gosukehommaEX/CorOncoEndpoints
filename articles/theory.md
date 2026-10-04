# Model and theory

## The model for one group

Let $`(Z_1, Z_2)`$ be standard bivariate normal with correlation
$`\theta`$.

- PFS is $`-\log\{\Phi(-Z_1)\}/\lambda_P`$, so it is exactly exponential
  with hazard $`\lambda_P`$.
- Response is $`R = 1`$ when PFS exceeds the landmark $`\tau`$ (0 by
  default) and $`Z_2 > c`$, where $`c`$ gives $`P(R = 1) = p`$ exactly.
  The copula correlation $`\theta`$ is calibrated to the Pearson
  correlation between PFS and $`R`$ or to the median PFS of responders.
- A PFS event is a death with probability $`\pi`$ and otherwise a
  progression. After progression the patient survives for a time $`D`$,
  and $`\mathrm{OS} = \mathrm{PFS} + \Delta D`$, where $`\Delta`$
  indicates progression as the first event. Hence
  $`\mathrm{PFS} \le \mathrm{OS}`$ for every patient.

Two models are available for $`D`$. In the illness-death model
(`os.model = "idm"`), $`D`$ is exponential with hazard $`\gamma_0`$ for
non-responders and $`\gamma_1 = \kappa \gamma_0`$ for responders
(`pps.hr.resp` $`= \kappa`$). The pre-progression death hazard is
$`\pi \lambda_P`$, so the post-progression hazard can be chosen
separately from it. In the exponential-exponential model
(`os.model = "expexp"`), OS is exactly exponential.

## Why OS is not exponential in the illness-death model

When the pre-progression hazards are constant and the post-progression
death hazard stays bounded right after progression, both PFS and OS can
be exponential only if the post-progression hazard equals the
pre-progression death hazard. The proportion of PFS events that are
deaths is then fixed at $`\pi = \mathrm{mPFS}/\mathrm{mOS}`$, and the
model is the maximal independence model of Fleischer et al. (2009). For
a PFS median of 6 months and an OS median of 15 months this proportion
is 0.4.

`os.model = "expexp"` keeps both distributions exponential with a
smaller $`\pi`$ by letting the pre-progression death hazard decrease
over time, $`\lambda_O e^{-c t}`$. The post-progression hazard (in time
since randomization) is then determined by the exponential OS
distribution:
``` math
h_{12}(t) = \frac{\lambda_O e^{-\lambda_O t} - h_{02}(t) e^{-\lambda_P t}}
{e^{-\lambda_O t} - e^{-\lambda_P t}},
```
so it is not a free input and cannot depend on response.

``` r

ee <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
              os.model = "expexp", os.median = 15)
SurvEndpoint(ee, c(6, 12, 24), endpoint = "os")
#> [1] 0.7578583 0.5743492 0.3298770
exp(-log(2) / 15 * c(6, 12, 24))
#> [1] 0.7578583 0.5743492 0.3298770
```

## Decomposition of the correlations

In the illness-death model, $`\Delta`$ and $`D`$ are independent of PFS
given response. With $`m_r = (1 - \pi)/\gamma_r`$ and
$`\delta = m_1 - m_0`$,
``` math
\mathrm{Cov}(\mathrm{OS}, R) = \mathrm{Cov}(\mathrm{PFS}, R) + p(1-p)\delta, \quad
\mathrm{Cov}(\mathrm{PFS}, \mathrm{OS}) = \mathrm{Var}(\mathrm{PFS}) + \delta\,\mathrm{Cov}(\mathrm{PFS}, R).
```
The partial covariance of OS and $`R`$ given PFS is
$`\delta p (1-p) \{1 - \mathrm{Corr}(\mathrm{PFS}, R)^2\}`$, so its sign
is the sign of $`\delta`$. Without a response effect on post-progression
survival ($`\kappa = 1`$),
$`\mathrm{Corr}(\mathrm{OS}, R) = \mathrm{Corr}(\mathrm{PFS}, R)
\,\mathrm{Corr}(\mathrm{PFS}, \mathrm{OS})`$.

``` r

k1 <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
              os.median = 15)
k06 <- OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
               os.median = 15, pps.hr.resp = 0.6)
rbind(kappa_1 = CorEndpoints(k1), kappa_0.6 = CorEndpoints(k06))
#>           cor.pfs.resp cor.os.resp cor.pfs.os pcor.os.resp
#> kappa_1            0.4   0.2435478  0.6088695 3.817591e-17
#> kappa_0.6          0.4   0.3808329  0.6034093 1.908293e-01
```

With $`\kappa = 1`$, $`\mathrm{Corr}(\mathrm{OS}, R)`$ is 0.244, the
product of 0.4 and 0.609. When responders also live longer after
progression ($`\kappa = 0.6`$), it rises to 0.381. The ordering
$`\mathrm{Corr}(\mathrm{PFS}, R) \ge \mathrm{Corr}(\mathrm{OS}, R)`$ is
therefore a consequence of the inputs and not an assumption.

## Attainable correlation between PFS and response

The Gaussian copula approaches the Frechet bounds as the copula
correlation tends to -1 or 1.

``` r

CorBoundPFSResponse(orr = 0.3)
#>      lower      upper 
#> -0.5448300  0.7881852
CorBoundPFSResponse(orr = 0.3, pfs.median = 6, resp.tau = 1.5)
#>      lower      upper 
#> -0.4073680  0.7881852
```

## References

Fleischer, F., Gaschler-Markefski, B. and Bluhmki, E. (2009). A
statistical model for the dependence between progression-free survival
and overall survival. *Statistics in Medicine*, 28, 2669–2686.

Meller, M., Beyersmann, J. and Rufibach, K. (2019). Joint modeling of
progression-free and overall survival and computation of correlation
measures. *Statistics in Medicine*, 38, 4270–4289.
