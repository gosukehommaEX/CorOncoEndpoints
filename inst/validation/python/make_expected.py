# Expected values used in tests/testthat (run: python3 make_expected.py).
# Requires numpy and scipy; writes expected_values.json in the current folder.
# Expected values for tests/testthat of CorOncoEndpoints (independent Python computation)
import json, numpy as np
from scipy import special, optimize, integrate
import reference as R
out = {}
scen = {
 "A": dict(pfs_median=6, orr=0.3, death_prop=0.15, resp_cor=0.4, os_model="idm", os_median=15, kappa=0.6),
 "B": dict(pfs_median=6, orr=0.3, death_prop=0.15, resp_pfs_median=11, os_model="idm", pps_median=8, kappa=1.0, tau=1.5),
 "C": dict(pfs_median=6, orr=0.3, death_prop=0.15, resp_cor=0.4, os_model="expexp", os_median=15, tau=1.5, ttr_median=2.5),
 "D": dict(pfs_median=6, orr=0.3, death_prop=0.4, resp_cor=0.4, os_model="expexp", os_median=15),
}
for k, kw in scen.items():
    A = R.make_arm(**kw)
    d = dict(theta=A.theta, c_resp=A.c, z_tau=A.ztau if np.isfinite(A.ztau) else None)
    if A.os_model == "idm":
        d.update(gam0=A.gam0, gam1=A.gam1)
    else:
        d.update(c_dec=A.cdec, lam_o=A.lamO)
    d.update({kk: float(v) for kk, v in R.correlations(A).items()})
    for resp in ["all", "resp", "nonresp"]:
        d["os_surv12_" + resp] = float(R.surv_os(12.0, A, response=resp))
        d["os_median_" + resp] = float(R.quantile(0.5, A, "os", resp))
        d["pfs_median_" + resp] = float(R.quantile(0.5, A, "pfs", resp))
    d["pfs_surv12_resp"] = float(R.surv_pfs_resp(12.0, A))
    d["os_dens12_all"] = float(R.dens_os(12.0, A))
    d["os_haz12_all"] = d["os_dens12_all"] / d["os_surv12_all"]
    d["pfs_dens12_resp"] = float(A.lamP * np.exp(-A.lamP * 12) * R.q_of_s(12.0, A) / A.p)
    out[k] = d
# Fleischer 2009 and TrialSimulator documented examples (kappa = 1 illness-death)
def fl(l1, l2, l3):
    return dict(pfs_hazard=l1 + l2, death_prop=l2 / (l1 + l2), pps_hazard=l3)
def fl_arm(l1, l2, l3):
    return R.make_arm(pfs_median=np.log(2) / (l1 + l2), orr=0.3, death_prop=l2 / (l1 + l2), resp_cor=0.0,
                      os_model="idm", pps_hazard=l3, kappa=1.0)
out["fleischer"] = {
  "fig4_corr": float(R.correlations(fl_arm(2, 0.5, 1.5))["cor_pfs_os"]),
  "ex1_corr": float(R.correlations(fl_arm(0.284, 0.075, 0.128))["cor_pfs_os"]),
  "ex3_ctl_median": float(R.quantile(0.5, fl_arm(0.342, 0.0654, 0.0797), "os")),
  "ex3_ctl_l3swap_median": float(R.quantile(0.5, fl_arm(0.342, 0.0654, 0.0876), "os")),
  "ex3_trt_median": float(R.quantile(0.5, fl_arm(0.139, 0.0648, 0.0876), "os")),
}
ts = fl_arm(0.1, 0.05, 0.12)
out["trialsimulator"] = {"corr": float(R.correlations(ts)["cor_pfs_os"]), "median_pfs": float(np.log(2) / 0.15),
                         "median_os": float(R.quantile(0.5, ts, "os"))}
# Correlation bounds
p, lam, tau = 0.3, np.log(2) / 6, 1.5
up = -np.sqrt(p / (1 - p)) * np.log(p)
lo0 = np.sqrt((1 - p) / p) * np.log(1 - p)
a = lam * tau
eb = np.exp(-a) - p
lb = -np.log(eb)
lo_tau = ((a + 1) * np.exp(-a) - (lb + 1) * eb - p) / np.sqrt(p * (1 - p))
ztau = float(R.z_of_t(tau, lam))
out["bounds"] = dict(upper=up, lower_tau0=lo0, lower_tau=lo_tau,
                     gauss_lower_tau=float(R.corr_pfs_r_std(-0.9999, p, ztau)[0]),
                     gauss_upper_tau=float(R.corr_pfs_r_std(0.9999, p, ztau)[0]),
                     gauss_lower_tau0=float(R.corr_pfs_r_std(-0.9999, p, -np.inf)[0]))
json.dump(out, open("expected_values.json", "w"), indent=1)
for k, v in out.items():
    print(k)
    for kk, vv in v.items():
        print("   %-22s %.10g" % (kk, vv) if vv is not None else "   %-22s None" % kk)
