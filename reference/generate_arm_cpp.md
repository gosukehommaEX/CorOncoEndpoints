# Generate patient-level data for one treatment group (internal)

Internal C++ kernel of
[`rOncoEndpoints`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md).
Random numbers are drawn from the `dqrng` generator in the following
order, each as one block of `nsim * n` values: latent PFS score,
independent normal for response, uniform for the type of the PFS event,
standard exponential for post-progression survival, uniform for the time
to response, uniform for accrual, standard exponential for dropout. The
draws are then transformed by `coronco::transform_one()` (see
`transform_core.h`).

## Usage

``` r
generate_arm_cpp(
  nsim,
  n,
  lam_p,
  theta,
  c_resp,
  z_tau,
  tau,
  os_model,
  pi_death,
  gam0,
  gam1,
  lam_o,
  c_dec,
  grid_g,
  grid_h,
  ttr_mode,
  ttr_rate,
  a_time,
  a_cum,
  d_hazard
)
```

## Arguments

- nsim:

  Number of simulated trials.

- n:

  Number of patients per trial in this group.

- lam_p, theta, c_resp, z_tau, tau:

  Parameters of PFS and response.

- os_model:

  0 for "idm", 1 for "expexp".

- pi_death, gam0, gam1:

  Parameters of the illness-death model.

- lam_o, c_dec, grid_g, grid_h:

  Parameters of the exponential-exponential model and its tabulated
  cumulative post-progression hazard.

- ttr_mode:

  1 to generate the time to response.

- ttr_rate:

  Rate of the untruncated time to response.

- a_time, a_cum:

  Accrual breakpoints and cumulative probabilities.

- d_hazard:

  Dropout hazard (0 for no dropout).

## Value

A list of vectors.
