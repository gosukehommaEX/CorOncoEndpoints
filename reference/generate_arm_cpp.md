# Generate patient-level data for one treatment group (internal)

Internal C++ kernel of
[`rOncoEndpoints`](https://gosukehommaEX.github.io/CorOncoEndpoints/reference/rOncoEndpoints.md).
Random numbers are drawn from the `dqrng` generator in the following
order, each as one block of `nsim * n` values: latent PFS score,
independent normal for response, uniform for the type of the PFS event,
standard exponential for post-progression survival, uniform for the time
to response, uniform for accrual, standard exponential for dropout, and,
only when `pc_hazard > 0`, standard exponential for the censoring of
PFS. The draws are then transformed by `coronco::transform_one()` (see
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
  a_counts,
  d_time,
  d_hazard,
  d_cum,
  pc_hazard,
  ep_mode
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

- a_time:

  Accrual breakpoints.

- a_counts:

  Number of patients of the group in each accrual interval (empty when
  all patients enter at `a_time[1]`).

- d_time, d_hazard, d_cum:

  Left ends, hazards and cumulative hazards of the pieces of the dropout
  distribution (one piece for exponential dropout; a zero hazard for no
  dropout).

- pc_hazard:

  Hazard of the censoring of PFS (0 for none).

- ep_mode:

  Endpoints of the group: 0 for all three, 1 for PFS and OS (response is
  `NA`), 2 for PFS and response (OS is `NA`), 3 for OS and response (the
  exponential time linked to response is OS, and PFS is `NA`). The
  random numbers drawn do not depend on `ep_mode`.

## Value

A list of vectors.
