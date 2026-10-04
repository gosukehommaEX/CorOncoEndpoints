// Pure C++ transformation core of CorOncoEndpoints (no R headers).
//
// Given the random inputs of one patient, these functions return the patient's
// progression-free survival (PFS), overall survival (OS), progression flag,
// response and time to response. The model is described in ?OncoArm:
//   PFS  = -log(Phi(-Z1)) / lam_p
//   R    = 1{Z1 > z_tau} 1{theta Z1 + sqrt(1 - theta^2) W > c_resp}
//   idm   : death before progression with probability pi_death, and
//           post-progression survival Exp(gam1) for responders, Exp(gam0) otherwise
//   expexp: death before progression with probability (lam_o / lam_p) exp(-c_dec PFS),
//           post-progression hazard h12(t) through its cumulative hazard G(t),
//           tabulated on a uniform grid and extended linearly with slope lam_o
//   OS   = PFS (death first) or the post-progression death time
//   TTR  = tau + Exp(ttr_rate) truncated to [0, PFS - tau) for responders

#ifndef CORONCO_TRANSFORM_CORE_H
#define CORONCO_TRANSFORM_CORE_H

#include <cmath>
#include <limits>

namespace coronco {

struct ArmParams {
  double lam_p;
  double theta;
  double c_resp;
  double z_tau;
  double tau;
  int os_model;          // 0 = "idm", 1 = "expexp"
  double pi_death;
  double gam0;
  double gam1;
  double lam_o;
  double c_dec;
  const double* grid_g;  // cumulative hazard G at t = 0, h, 2h, ...
  int grid_n;
  double grid_h;
  int ttr_mode;          // 1 when the time to response is generated
  double ttr_rate;
};

// log(Phi(-z)) for a standard normal distribution function Phi
inline double log_upper_norm(double z) {
  const double rt2 = 1.4142135623730951;
  if (z < 0.0) {
    return std::log1p(-0.5 * std::erfc(-z / rt2));
  }
  double v = 0.5 * std::erfc(z / rt2);
  if (v > 0.0) {
    return std::log(v);
  }
  // asymptotic expansion for extreme upper tails
  return -0.5 * z * z - std::log(z) - 0.5 * std::log(2.0 * 3.141592653589793);
}

// Cumulative post-progression hazard G(t) by linear interpolation on the grid
inline double cumhaz_interp(double t, const ArmParams& p) {
  const double tmax = (p.grid_n - 1) * p.grid_h;
  if (t >= tmax) {
    return p.grid_g[p.grid_n - 1] + p.lam_o * (t - tmax);
  }
  const double x = t / p.grid_h;
  int i = (int) std::floor(x);
  if (i < 0) i = 0;
  if (i > p.grid_n - 2) i = p.grid_n - 2;
  const double f = x - i;
  return p.grid_g[i] + f * (p.grid_g[i + 1] - p.grid_g[i]);
}

// Inverse of G by binary search and linear interpolation
inline double cumhaz_inverse(double y, const ArmParams& p) {
  const double gmax = p.grid_g[p.grid_n - 1];
  const double tmax = (p.grid_n - 1) * p.grid_h;
  if (y >= gmax) {
    return tmax + (y - gmax) / p.lam_o;
  }
  int lo = 0;
  int hi = p.grid_n - 1;
  while (hi - lo > 1) {
    const int mid = lo + (hi - lo) / 2;
    if (p.grid_g[mid] <= y) lo = mid; else hi = mid;
  }
  const double g0 = p.grid_g[lo];
  const double g1 = p.grid_g[hi];
  const double f = (g1 > g0) ? (y - g0) / (g1 - g0) : 0.0;
  return (lo + f) * p.grid_h;
}

// Piecewise-uniform accrual: breakpoints a_time[0..k] and cumulative
// probabilities a_cum[0..k] with a_cum[0] = 0 and a_cum[k] = 1
inline double accrual_inverse(double u, const double* a_time, const double* a_cum, int k) {
  if (k <= 0) return 0.0;
  for (int j = 0; j < k; ++j) {
    const double lo = a_cum[j];
    const double hi = a_cum[j + 1];
    if (hi > lo && (u < hi || j == k - 1)) {
      double f = (u - lo) / (hi - lo);
      if (f < 0.0) f = 0.0;
      if (f > 1.0) f = 1.0;
      return a_time[j] + f * (a_time[j + 1] - a_time[j]);
    }
  }
  return a_time[k];
}

// Transform the random inputs of one patient
inline void transform_one(double z1, double w, double u_type, double e_pps, double u_ttr,
                          const ArmParams& p,
                          double& pfs, double& os, int& prog, int& resp, double& ttr) {
  pfs = -log_upper_norm(z1) / p.lam_p;
  const double z2 = p.theta * z1 + std::sqrt(1.0 - p.theta * p.theta) * w;
  resp = (z1 > p.z_tau && z2 > p.c_resp) ? 1 : 0;

  bool death_first;
  if (p.os_model == 0) {
    death_first = u_type < p.pi_death;
  } else {
    death_first = u_type < (p.lam_o / p.lam_p) * std::exp(-p.c_dec * pfs);
  }
  prog = death_first ? 0 : 1;

  if (death_first) {
    os = pfs;
  } else if (p.os_model == 0) {
    const double rate = (resp == 1) ? p.gam1 : p.gam0;
    os = pfs + e_pps / rate;
  } else {
    os = cumhaz_inverse(cumhaz_interp(pfs, p) + e_pps, p);
    if (os < pfs) os = pfs;
  }

  ttr = std::numeric_limits<double>::quiet_NaN();
  if (p.ttr_mode == 1 && resp == 1) {
    const double width = pfs - p.tau;
    ttr = p.tau - std::log1p(-u_ttr * (-std::expm1(-p.ttr_rate * width))) / p.ttr_rate;
  }
}

}  // namespace coronco

#endif
