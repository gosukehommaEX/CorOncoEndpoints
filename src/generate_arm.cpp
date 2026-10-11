// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(dqrng)]]
#include <Rcpp.h>
#include <dqrng.h>
#include <cmath>
#include <vector>
#include "transform_core.h"
using namespace Rcpp;

//' Generate patient-level data for one treatment group (internal)
//'
//' Internal C++ kernel of \code{\link{rOncoEndpoints}}. Random numbers are
//' drawn from the \code{dqrng} generator in the following order, each as one
//' block of \code{nsim * n} values: latent PFS score, independent normal for
//' response, uniform for the type of the PFS event, standard exponential for
//' post-progression survival, uniform for the time to response, uniform for
//' accrual, standard exponential for dropout, and, only when
//' \code{pc_hazard > 0}, standard exponential for the censoring of PFS. The
//' draws are then transformed by \code{coronco::transform_one()} (see
//' \code{transform_core.h}).
//'
//' @param nsim Number of simulated trials.
//' @param n Number of patients per trial in this group.
//' @param lam_p,theta,c_resp,z_tau,tau Parameters of PFS and response.
//' @param os_model 0 for "idm", 1 for "expexp".
//' @param pi_death,gam0,gam1 Parameters of the illness-death model.
//' @param lam_o,c_dec,grid_g,grid_h Parameters of the exponential-exponential
//'   model and its tabulated cumulative post-progression hazard.
//' @param ttr_mode 1 to generate the time to response.
//' @param ttr_rate Rate of the untruncated time to response.
//' @param a_time Accrual breakpoints.
//' @param a_counts Number of patients of the group in each accrual interval
//'   (empty when all patients enter at \code{a_time[1]}).
//' @param d_time,d_hazard,d_cum Left ends, hazards and cumulative hazards of
//'   the pieces of the dropout distribution (one piece for exponential
//'   dropout; a zero hazard for no dropout).
//' @param pc_hazard Hazard of the censoring of PFS (0 for none).
//' @param ep_mode Endpoints of the group: 0 for all three, 1 for PFS and OS
//'   (response is \code{NA}), 2 for PFS and response (OS is \code{NA}), 3 for
//'   OS and response (the exponential time linked to response is OS, and PFS
//'   is \code{NA}). The random numbers drawn do not depend on \code{ep_mode}.
//' @return A list of vectors.
//' @keywords internal
// [[Rcpp::export]]
List generate_arm_cpp(int nsim, int n, double lam_p, double theta, double c_resp,
                      double z_tau, double tau, int os_model, double pi_death,
                      double gam0, double gam1, double lam_o, double c_dec,
                      NumericVector grid_g, double grid_h, int ttr_mode,
                      double ttr_rate, NumericVector a_time, IntegerVector a_counts,
                      NumericVector d_time, NumericVector d_hazard,
                      NumericVector d_cum, double pc_hazard, int ep_mode) {
  const int N = nsim * n;
  NumericVector z1 = dqrng::dqrnorm(N, 0.0, 1.0);
  NumericVector w = dqrng::dqrnorm(N, 0.0, 1.0);
  NumericVector u_type = dqrng::dqrunif(N, 0.0, 1.0);
  NumericVector e_pps = dqrng::dqrexp(N, 1.0);
  NumericVector u_ttr = dqrng::dqrunif(N, 0.0, 1.0);
  NumericVector u_acc = dqrng::dqrunif(N, 0.0, 1.0);
  NumericVector e_drop = dqrng::dqrexp(N, 1.0);
  NumericVector e_pc;
  if (pc_hazard > 0.0) e_pc = dqrng::dqrexp(N, 1.0);

  coronco::ArmParams par;
  par.lam_p = lam_p;
  par.theta = theta;
  par.c_resp = c_resp;
  par.z_tau = z_tau;
  par.tau = tau;
  par.os_model = os_model;
  par.pi_death = pi_death;
  par.gam0 = gam0;
  par.gam1 = gam1;
  par.lam_o = lam_o;
  par.c_dec = c_dec;
  par.grid_g = REAL(grid_g);
  par.grid_n = grid_g.size();
  par.grid_h = grid_h;
  par.ttr_mode = ttr_mode;
  par.ttr_rate = ttr_rate;
  // accrual interval of each patient of a simulated trial
  const int k_acc = a_counts.size();
  std::vector<int> interval(n, 0);
  if (k_acc > 0) {
    if (a_time.size() != k_acc + 1) stop("'a_time' must have one more element than 'a_counts'.");
    int total = 0;
    for (int j = 0; j < k_acc; ++j) total += a_counts[j];
    if (total != n) stop("The accrual counts must add up to the group size.");
    int r = 0;
    for (int j = 0; j < k_acc; ++j) {
      for (int c = 0; c < a_counts[j]; ++c) interval[r++] = j;
    }
  }
  const int k_drop = d_hazard.size();
  if (k_drop < 1 || d_time.size() != k_drop || d_cum.size() != k_drop) {
    stop("The dropout pieces are inconsistent.");
  }

  IntegerVector sim(N), progression(N), response(N);
  NumericVector accrual(N), pfs(N), os(N), ttr(N), dropout(N), pfs_cens(N);
  for (int i = 0; i < N; ++i) {
    double v_pfs, v_os, v_ttr;
    int v_prog, v_resp;
    coronco::transform_one(z1[i], w[i], u_type[i], e_pps[i], u_ttr[i], par,
                           v_pfs, v_os, v_prog, v_resp, v_ttr);
    if (ep_mode == 1) {
      v_resp = NA_INTEGER;
      v_ttr = NA_REAL;
    } else if (ep_mode == 2) {
      v_os = NA_REAL;
      v_prog = NA_INTEGER;
    } else if (ep_mode == 3) {
      v_os = v_pfs;
      v_pfs = NA_REAL;
      v_prog = NA_INTEGER;
    }
    sim[i] = i / n + 1;
    pfs[i] = v_pfs;
    os[i] = v_os;
    progression[i] = v_prog;
    response[i] = v_resp;
    ttr[i] = std::isnan(v_ttr) ? NA_REAL : v_ttr;
    accrual[i] = (k_acc > 0) ?
      coronco::accrual_position(u_acc[i], interval[i % n], REAL(a_time)) : a_time[0];
    dropout[i] = coronco::piecewise_exp_inverse(e_drop[i], REAL(d_time), REAL(d_hazard),
                                                REAL(d_cum), k_drop);
    pfs_cens[i] = (pc_hazard > 0.0) ? e_pc[i] / pc_hazard : R_PosInf;
  }
  return List::create(Named("sim") = sim, Named("accrual_time") = accrual,
                      Named("pfs_time") = pfs, Named("os_time") = os,
                      Named("progression") = progression,
                      Named("response") = response, Named("ttr") = ttr,
                      Named("dropout_time") = dropout,
                      Named("pfs_cens_time") = pfs_cens);
}
